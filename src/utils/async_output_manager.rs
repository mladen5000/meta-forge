//! High-Performance Async Output Manager
//! =====================================
//!
//! Optimized output manager that eliminates file I/O bottlenecks through:
//! - Async/parallel file operations
//! - Batched writes with buffering
//! - Optional compression (disabled by default)
//! - Minimal metadata overhead
//! - Progress tracking for large operations

use anyhow::{Context, Result};
use chrono::Utc;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tokio::fs;
use tokio::io::{AsyncWriteExt, BufWriter};
use tokio::sync::{Mutex, Semaphore};
use tracing::{info, debug};

use crate::utils::intermediate_output::PipelineSection;

/// High-performance async output manager
#[derive(Debug, Clone)]
pub struct AsyncOutputManager {
    /// Base output directory
    pub base_output_dir: PathBuf,
    /// Current run directory
    pub run_dir: PathBuf,
    /// Run ID for this execution
    pub run_id: String,
    /// Configuration for output formats
    pub config: FastOutputConfig,
    /// Write semaphore to limit concurrent operations
    write_semaphore: Arc<Semaphore>,
    /// Buffer for batched operations
    write_buffer: Arc<Mutex<WriteBuffer>>,
}

/// Fast output configuration with performance optimizations
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastOutputConfig {
    /// Enable JSON output (primary format)
    pub enable_json: bool,
    /// Enable FASTA output for sequences (only when needed)
    pub enable_fasta: bool,
    /// Enable TSV output for tabular data (only when needed)
    pub enable_tsv: bool,
    /// Disable compression by default for speed
    pub compress_files: bool,
    /// Buffer size for writes (KB)
    pub buffer_size_kb: usize,
    /// Maximum concurrent writes
    pub max_concurrent_writes: usize,
    /// Skip checksum calculation for speed
    pub skip_checksums: bool,
    /// Batch write threshold
    pub batch_threshold: usize,
}

impl Default for FastOutputConfig {
    fn default() -> Self {
        Self {
            enable_json: true,
            enable_fasta: false, // Only enable when needed
            enable_tsv: false,   // Only enable when needed
            compress_files: false, // Disabled for speed
            buffer_size_kb: 64,    // 64KB buffer
            max_concurrent_writes: 8,
            skip_checksums: true, // Skip for speed
            batch_threshold: 10,  // Batch 10 operations
        }
    }
}

/// Write operation for batching
#[derive(Debug)]
struct WriteOperation {
    path: PathBuf,
    content: Vec<u8>,
    section: PipelineSection,
    format: String,
}

/// Write buffer for batching operations
#[derive(Debug, Default)]
struct WriteBuffer {
    operations: Vec<WriteOperation>,
    total_size: usize,
}

impl AsyncOutputManager {
    /// Create a new async output manager
    pub async fn new(base_output_dir: PathBuf, config: FastOutputConfig) -> Result<Self> {
        let now = Utc::now();
        let run_id = now.format("%d%m%y_%H%M%S").to_string();
        let run_dir = base_output_dir.join(format!("run_{run_id}"));

        // Create base directories
        fs::create_dir_all(&run_dir)
            .await
            .with_context(|| format!("Failed to create run directory: {}", run_dir.display()))?;

        info!("🗂️ Created run directory: {}", run_dir.display());

        let manager = Self {
            base_output_dir,
            run_dir,
            run_id,
            write_semaphore: Arc::new(Semaphore::new(config.max_concurrent_writes)),
            write_buffer: Arc::new(Mutex::new(WriteBuffer::default())),
            config,
        };

        // Create section subdirectories in parallel
        manager.create_section_directories().await?;

        Ok(manager)
    }

    /// Create all section subdirectories in parallel
    async fn create_section_directories(&self) -> Result<()> {
        let sections = [
            PipelineSection::Preprocessing,
            PipelineSection::Preprocessing,
            PipelineSection::Assembly,
            PipelineSection::Classification,
            PipelineSection::Classification,
            PipelineSection::Abundance,
            PipelineSection::Report,
        ];

        // Create directories in parallel
        let tasks: Vec<_> = sections
            .iter()
            .map(|section| async move {
                let section_dir = self.get_section_dir(section);
                fs::create_dir_all(&section_dir).await.with_context(|| {
                    format!("Failed to create section directory: {}", section_dir.display())
                })?;
                debug!("📁 Created section directory: {}", section_dir.display());
                Ok::<(), anyhow::Error>(())
            })
            .collect();

        // Wait for all directory creation tasks
        for task in tasks {
            task.await?;
        }

        Ok(())
    }

    /// Get the directory path for a specific section
    pub fn get_section_dir(&self, section: &PipelineSection) -> PathBuf {
        self.run_dir.join(section.as_str())
    }

    /// Save intermediate data asynchronously with batching
    pub async fn save_intermediate<T: Serialize>(
        &self,
        section: PipelineSection,
        filename: &str,
        data: &T,
        progress_callback: Option<Box<dyn Fn(f32) + Send + Sync>>,
    ) -> Result<()> {
        let section_dir = self.get_section_dir(&section);
        let json_path = section_dir.join(format!("{filename}.json"));

        // Serialize data
        let json_content = serde_json::to_string(data)
            .context("Failed to serialize data to JSON")?;

        if let Some(callback) = &progress_callback {
            callback(0.5);
        }

        // Add to batch or write immediately
        if self.config.batch_threshold > 1 {
            self.add_to_batch(WriteOperation {
                path: json_path,
                content: json_content.into_bytes(),
                section,
                format: "json".to_string(),
            }).await?;
        } else {
            self.write_file_direct(&json_path, json_content.as_bytes()).await?;
            debug!("💾 Saved JSON: {}", json_path.display());
        }

        if let Some(callback) = progress_callback {
            callback(1.0);
        }

        Ok(())
    }

    /// Save sequences asynchronously (only when FASTA is enabled)
    pub async fn save_sequences(
        &self,
        section: PipelineSection,
        filename: &str,
        sequences: &[(String, String)],
        progress_callback: Option<Box<dyn Fn(f32) + Send + Sync>>,
    ) -> Result<()> {
        if !self.config.enable_fasta {
            debug!("FASTA output disabled, skipping sequence save");
            return Ok(());
        }

        let section_dir = self.get_section_dir(&section);
        let fasta_path = section_dir.join(format!("{filename}.fasta"));

        // Build FASTA content efficiently
        let mut fasta_content = String::with_capacity(sequences.len() * 100); // Pre-allocate
        for (i, (id, sequence)) in sequences.iter().enumerate() {
            fasta_content.push_str(&format!(">{id}\n{sequence}\n"));

            if let Some(callback) = &progress_callback {
                if i % 1000 == 0 {
                    callback(i as f32 / sequences.len() as f32);
                }
            }
        }

        self.write_file_direct(&fasta_path, fasta_content.as_bytes()).await?;
        debug!("💾 Saved FASTA: {}", fasta_path.display());

        if let Some(callback) = progress_callback {
            callback(1.0);
        }

        Ok(())
    }

    /// Save TSV data asynchronously (only when TSV is enabled)
    pub async fn save_tsv(
        &self,
        section: PipelineSection,
        filename: &str,
        headers: &[String],
        rows: &[Vec<String>],
        progress_callback: Option<Box<dyn Fn(f32) + Send + Sync>>,
    ) -> Result<()> {
        if !self.config.enable_tsv {
            debug!("TSV output disabled, skipping TSV save");
            return Ok(());
        }

        let section_dir = self.get_section_dir(&section);
        let tsv_path = section_dir.join(format!("{filename}.tsv"));

        // Build TSV content efficiently
        let mut tsv_content = String::with_capacity(rows.len() * 50); // Pre-allocate
        tsv_content.push_str(&headers.join("\t"));
        tsv_content.push('\n');

        for (i, row) in rows.iter().enumerate() {
            tsv_content.push_str(&row.join("\t"));
            tsv_content.push('\n');

            if let Some(callback) = &progress_callback {
                if i % 1000 == 0 {
                    callback(i as f32 / rows.len() as f32);
                }
            }
        }

        self.write_file_direct(&tsv_path, tsv_content.as_bytes()).await?;
        debug!("💾 Saved TSV: {}", tsv_path.display());

        if let Some(callback) = progress_callback {
            callback(1.0);
        }

        Ok(())
    }

    /// Add operation to batch buffer
    async fn add_to_batch(&self, operation: WriteOperation) -> Result<()> {
        let mut buffer = self.write_buffer.lock().await;
        buffer.total_size += operation.content.len();
        buffer.operations.push(operation);

        // Flush if threshold reached
        if buffer.operations.len() >= self.config.batch_threshold {
            let operations = std::mem::take(&mut buffer.operations);
            buffer.total_size = 0;
            drop(buffer); // Release lock

            self.flush_batch(operations).await?;
        }

        Ok(())
    }

    /// Flush batched operations in parallel
    async fn flush_batch(&self, operations: Vec<WriteOperation>) -> Result<()> {
        if operations.is_empty() {
            return Ok(());
        }

        debug!("🔄 Flushing {} batched write operations", operations.len());

        // Write operations in parallel with semaphore limiting
        let tasks: Vec<_> = operations
            .into_iter()
            .map(|op| {
                let semaphore = self.write_semaphore.clone();
                async move {
                    let _permit = semaphore.acquire().await.unwrap();
                    self.write_file_direct(&op.path, &op.content).await
                        .with_context(|| format!("Failed to write file: {}", op.path.display()))
                }
            })
            .collect();

        // Wait for all writes to complete
        for task in tasks {
            task.await?;
        }

        debug!("✅ Batch flush completed");
        Ok(())
    }

    /// Write file directly with optimizations
    async fn write_file_direct(&self, path: &Path, content: &[u8]) -> Result<()> {
        // Create parent directory if needed
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).await?;
        }

        // Use buffered writer for better performance
        let file = fs::File::create(path).await?;
        let mut writer = BufWriter::with_capacity(
            self.config.buffer_size_kb * 1024,
            file
        );

        writer.write_all(content).await?;
        writer.flush().await?;

        Ok(())
    }

    /// Force flush all pending operations
    pub async fn flush_all(&self) -> Result<()> {
        let mut buffer = self.write_buffer.lock().await;
        let operations = std::mem::take(&mut buffer.operations);
        buffer.total_size = 0;
        drop(buffer);

        if !operations.is_empty() {
            self.flush_batch(operations).await?;
        }

        Ok(())
    }

    /// Get run statistics
    pub async fn get_stats(&self) -> OutputStats {
        let buffer = self.write_buffer.lock().await;
        OutputStats {
            pending_operations: buffer.operations.len(),
            pending_size_bytes: buffer.total_size,
            max_concurrent_writes: self.config.max_concurrent_writes,
            compression_enabled: self.config.compress_files,
        }
    }

    /// Convert to synchronous output manager for compatibility
    pub fn to_sync_manager(&self) -> crate::utils::intermediate_output::IntermediateOutputManager {
        let sync_config = crate::utils::intermediate_output::OutputConfig {
            enable_json: self.config.enable_json,
            enable_binary: false,
            enable_fasta: self.config.enable_fasta,
            enable_tsv: self.config.enable_tsv,
            compress_files: self.config.compress_files,
            max_file_size_mb: 100,
        };

        crate::utils::intermediate_output::IntermediateOutputManager {
            base_output_dir: self.base_output_dir.clone(),
            run_dir: self.run_dir.clone(),
            run_timestamp: Utc::now(),
            run_id: self.run_id.clone(),
            config: sync_config,
        }
    }
}

/// Output performance statistics
#[derive(Debug, Clone)]
pub struct OutputStats {
    pub pending_operations: usize,
    pub pending_size_bytes: usize,
    pub max_concurrent_writes: usize,
    pub compression_enabled: bool,
}

/// Convenience macro for progress tracking
#[macro_export]
macro_rules! track_progress {
    ($callback:expr, $progress:expr) => {
        if let Some(ref cb) = $callback {
            cb($progress);
        }
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;
    use serde_json::json;

    #[tokio::test]
    async fn test_async_output_manager_creation() {
        let temp_dir = tempdir().unwrap();
        let config = FastOutputConfig::default();

        let manager = AsyncOutputManager::new(temp_dir.path().to_path_buf(), config)
            .await
            .unwrap();

        // Check that run directory was created
        assert!(manager.run_dir.exists());

        // Check that section directories were created
        for section in [
            PipelineSection::Preprocessing,
            PipelineSection::Assembly,
            PipelineSection::Classification,
        ] {
            let section_dir = manager.get_section_dir(&section);
            assert!(section_dir.exists());
        }
    }

    #[tokio::test]
    async fn test_async_save_intermediate() {
        let temp_dir = tempdir().unwrap();
        let config = FastOutputConfig {
            batch_threshold: 1, // Disable batching for test
            ..FastOutputConfig::default()
        };

        let manager = AsyncOutputManager::new(temp_dir.path().to_path_buf(), config)
            .await
            .unwrap();

        let test_data = json!({"test": "data", "value": 42});

        // Save data
        manager
            .save_intermediate(PipelineSection::Assembly, "test_data", &test_data, None)
            .await
            .unwrap();

        // Check that file was created
        let json_path = manager
            .get_section_dir(&PipelineSection::Assembly)
            .join("test_data.json");
        assert!(json_path.exists());

        // Verify content
        let content = fs::read_to_string(&json_path).await.unwrap();
        let loaded: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(loaded, test_data);
    }

    #[tokio::test]
    async fn test_batched_writes() {
        let temp_dir = tempdir().unwrap();
        let config = FastOutputConfig {
            batch_threshold: 3, // Batch every 3 operations
            ..FastOutputConfig::default()
        };

        let manager = AsyncOutputManager::new(temp_dir.path().to_path_buf(), config)
            .await
            .unwrap();

        // Add multiple operations (should trigger batching)
        for i in 0..5 {
            let test_data = json!({"test": i});
            manager
                .save_intermediate(
                    PipelineSection::Assembly,
                    &format!("test_{}", i),
                    &test_data,
                    None,
                )
                .await
                .unwrap();
        }

        // Flush remaining operations
        manager.flush_all().await.unwrap();

        // Check that all files were created
        for i in 0..5 {
            let json_path = manager
                .get_section_dir(&PipelineSection::Assembly)
                .join(format!("test_{}.json", i));
            assert!(json_path.exists());
        }
    }

    #[tokio::test]
    async fn test_performance_config() {
        let temp_dir = tempdir().unwrap();
        let config = FastOutputConfig {
            enable_fasta: false,
            enable_tsv: false,
            compress_files: false,
            skip_checksums: true,
            max_concurrent_writes: 16,
            ..FastOutputConfig::default()
        };

        let manager = AsyncOutputManager::new(temp_dir.path().to_path_buf(), config)
            .await
            .unwrap();

        let stats = manager.get_stats().await;
        assert!(!stats.compression_enabled);
        assert_eq!(stats.max_concurrent_writes, 16);
    }
}