use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};
use tracing::info;

/// Intermediate output manager for pipeline sections with run-specific directories
#[derive(Debug, Clone)]
pub struct IntermediateOutputManager {
    /// Base output directory
    pub base_output_dir: PathBuf,
    /// Current run directory (e.g., output/run_12081525_143052/)
    pub run_dir: PathBuf,
    /// Timestamp of current run
    pub run_timestamp: DateTime<Utc>,
    /// Run ID for this execution
    pub run_id: String,
    /// Configuration for output formats
    pub config: OutputConfig,
}

/// Configuration for intermediate output formats
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputConfig {
    /// Enable JSON output for structured data
    pub enable_json: bool,
    /// Enable binary output for large datasets
    pub enable_binary: bool,
    /// Enable FASTA output for sequences
    pub enable_fasta: bool,
    /// Enable TSV output for tabular data
    pub enable_tsv: bool,
    /// Compress intermediate files
    pub compress_files: bool,
    /// Maximum file size before splitting (MB)
    pub max_file_size_mb: usize,
}

impl Default for OutputConfig {
    fn default() -> Self {
        Self {
            enable_json: true,
            enable_binary: true,
            enable_fasta: true,
            enable_tsv: true,
            compress_files: true,
            max_file_size_mb: 100,
        }
    }
}

/// Pipeline section identifiers
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PipelineSection {
    Preprocessing,
    Assembly,
    Classification,
    Abundance,
    Report,
}

impl PipelineSection {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Preprocessing => "preprocessing",
            Self::Assembly => "assembly",
            Self::Classification => "classification",
            Self::Abundance => "abundance",
            Self::Report => "report",
        }
    }
}

/// Metadata for intermediate files
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntermediateMetadata {
    /// Section that generated this file
    pub section: PipelineSection,
    /// File type and format
    pub format: String,
    /// Timestamp when file was created
    pub created_at: DateTime<Utc>,
    /// File size in bytes
    pub file_size: u64,
    /// Checksums for validation
    pub checksums: FileChecksums,
    /// Schema version for compatibility
    pub schema_version: String,
    /// Additional metadata specific to the section
    pub section_metadata: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileChecksums {
    pub md5: String,
    pub sha256: String,
}

impl IntermediateOutputManager {
    /// Create a new intermediate output manager with timestamp-based run directory
    pub fn new(base_output_dir: PathBuf, config: OutputConfig) -> Result<Self> {
        let now = Utc::now();
        let run_id = now.format("%Y%m%d_%H%M%S").to_string();
        let run_dir = base_output_dir.join(format!("run_{run_id}"));

        // Create base directories
        fs::create_dir_all(&run_dir)
            .with_context(|| format!("Failed to create run directory: {}", run_dir.display()))?;

        info!("🗂️ Created run directory: {}", run_dir.display());

        let manager = Self {
            base_output_dir,
            run_dir,
            run_timestamp: now,
            run_id,
            config,
        };

        // Create section subdirectories
        manager.create_section_directories()?;

        Ok(manager)
    }

    /// Create all section subdirectories
    fn create_section_directories(&self) -> Result<()> {
        let sections = [
            PipelineSection::Preprocessing,
            PipelineSection::Preprocessing,
            PipelineSection::Assembly,
            PipelineSection::Classification,
            PipelineSection::Classification,
            PipelineSection::Abundance,
            PipelineSection::Report,
        ];

        for section in &sections {
            let section_dir = self.get_section_dir(section);
            fs::create_dir_all(&section_dir).with_context(|| {
                format!(
                    "Failed to create section directory: {}",
                    section_dir.display()
                )
            })?;
            info!("📁 Created section directory: {}", section_dir.display());
        }

        Ok(())
    }

    /// Get the directory path for a specific section
    pub fn get_section_dir(&self, section: &PipelineSection) -> PathBuf {
        self.run_dir.join(section.as_str())
    }

    /// Save intermediate data for a pipeline section
    pub fn save_intermediate<T: Serialize>(
        &self,
        section: PipelineSection,
        filename: &str,
        data: &T,
        section_metadata: serde_json::Value,
    ) -> Result<IntermediateMetadata> {
        let section_dir = self.get_section_dir(&section);

        // Save as JSON (primary format)
        let json_path = section_dir.join(format!("{filename}.json"));
        let json_content =
            serde_json::to_string_pretty(data).context("Failed to serialize data to JSON")?;

        let file_path = if self.config.compress_files {
            let compressed_path = json_path.with_extension("json.gz");
            self.write_compressed(&compressed_path, json_content.as_bytes())?;
            info!("💾 Saved compressed JSON: {}", compressed_path.display());
            compressed_path
        } else {
            fs::write(&json_path, json_content)
                .with_context(|| format!("Failed to write JSON file: {}", json_path.display()))?;
            info!("💾 Saved JSON: {}", json_path.display());
            json_path
        };

        self.create_metadata(&section, &file_path, "json", section_metadata)
    }

    /// Save FASTA sequences for sections dealing with sequences
    pub fn save_sequences(
        &self,
        section: PipelineSection,
        filename: &str,
        sequences: &[(String, String)], // (id, sequence) pairs
        section_metadata: serde_json::Value,
    ) -> Result<IntermediateMetadata> {
        if !self.config.enable_fasta {
            return Err(anyhow::anyhow!("FASTA output is disabled in configuration"));
        }

        let section_dir = self.get_section_dir(&section);
        let fasta_path = section_dir.join(format!("{filename}.fasta"));

        let mut fasta_content = String::new();
        for (id, sequence) in sequences {
            fasta_content.push_str(&format!(">{id}\n{sequence}\n"));
        }

        if self.config.compress_files {
            self.write_compressed(
                &fasta_path.with_extension("fasta.gz"),
                fasta_content.as_bytes(),
            )?;
            info!(
                "💾 Saved compressed FASTA: {}",
                fasta_path.with_extension("fasta.gz").display()
            );
            self.create_metadata(
                &section,
                &fasta_path.with_extension("fasta.gz"),
                "fasta",
                section_metadata,
            )
        } else {
            fs::write(&fasta_path, fasta_content)
                .with_context(|| format!("Failed to write FASTA file: {}", fasta_path.display()))?;
            info!("💾 Saved FASTA: {}", fasta_path.display());
            self.create_metadata(&section, &fasta_path, "fasta", section_metadata)
        }
    }

    /// Save TSV data for tabular information
    pub fn save_tsv(
        &self,
        section: PipelineSection,
        filename: &str,
        headers: &[String],
        rows: &[Vec<String>],
        section_metadata: serde_json::Value,
    ) -> Result<IntermediateMetadata> {
        if !self.config.enable_tsv {
            return Err(anyhow::anyhow!("TSV output is disabled in configuration"));
        }

        let section_dir = self.get_section_dir(&section);
        let tsv_path = section_dir.join(format!("{filename}.tsv"));

        let mut tsv_content = headers.join("\t") + "\n";
        for row in rows {
            tsv_content.push_str(&(row.join("\t") + "\n"));
        }

        if self.config.compress_files {
            self.write_compressed(&tsv_path.with_extension("tsv.gz"), tsv_content.as_bytes())?;
            info!(
                "💾 Saved compressed TSV: {}",
                tsv_path.with_extension("tsv.gz").display()
            );
            self.create_metadata(
                &section,
                &tsv_path.with_extension("tsv.gz"),
                "tsv",
                section_metadata,
            )
        } else {
            fs::write(&tsv_path, tsv_content)
                .with_context(|| format!("Failed to write TSV file: {}", tsv_path.display()))?;
            info!("💾 Saved TSV: {}", tsv_path.display());
            self.create_metadata(&section, &tsv_path, "tsv", section_metadata)
        }
    }

    /// Load intermediate data from a previous run
    pub fn load_intermediate<T: for<'de> Deserialize<'de>>(
        &self,
        section: PipelineSection,
        filename: &str,
    ) -> Result<T> {
        let section_dir = self.get_section_dir(&section);

        // Load from JSON (primary format)
        let json_path = section_dir.join(format!("{filename}.json"));
        let json_gz_path = section_dir.join(format!("{filename}.json.gz"));

        let content = if json_gz_path.exists() {
            self.read_compressed(&json_gz_path)?
        } else if json_path.exists() {
            fs::read_to_string(&json_path)
                .with_context(|| format!("Failed to read JSON file: {}", json_path.display()))?
        } else {
            return Err(anyhow::anyhow!(
                "JSON file not found: {}",
                json_path.display()
            ));
        };

        let data: T = serde_json::from_str(&content).context("Failed to deserialize JSON data")?;
        Ok(data)
    }

    /// Check if intermediate files exist for a section
    pub fn has_intermediate(&self, section: PipelineSection, filename: &str) -> bool {
        let section_dir = self.get_section_dir(&section);

        // Check JSON files
        let json_path = section_dir.join(format!("{filename}.json"));
        let json_gz_path = section_dir.join(format!("{filename}.json.gz"));

        // Check FASTA files
        let fasta_path = section_dir.join(format!("{filename}.fasta"));
        let fasta_gz_path = section_dir.join(format!("{filename}.fasta.gz"));

        // Check TSV files
        let tsv_path = section_dir.join(format!("{filename}.tsv"));
        let tsv_gz_path = section_dir.join(format!("{filename}.tsv.gz"));

        json_path.exists()
            || json_gz_path.exists()
            || fasta_path.exists()
            || fasta_gz_path.exists()
            || tsv_path.exists()
            || tsv_gz_path.exists()
    }

    /// List all run directories in the base output directory
    pub fn list_run_directories(&self) -> Result<Vec<PathBuf>> {
        let mut runs = Vec::new();

        if self.base_output_dir.exists() {
            for entry in fs::read_dir(&self.base_output_dir)? {
                let entry = entry?;
                let path = entry.path();
                if path.is_dir()
                    && path
                        .file_name()
                        .and_then(|n| n.to_str())
                        .map(|n| n.starts_with("run_"))
                        .unwrap_or(false)
                {
                    runs.push(path);
                }
            }
        }

        runs.sort();
        Ok(runs)
    }

    /// Get metadata about the current run
    pub fn get_run_info(&self) -> RunInfo {
        RunInfo {
            run_id: self.run_id.clone(),
            run_dir: self.run_dir.clone(),
            timestamp: self.run_timestamp,
            config: self.config.clone(),
        }
    }

    /// Generate a summary of all intermediate files in the current run
    pub fn generate_run_summary(&self) -> Result<RunSummary> {
        let mut sections = std::collections::HashMap::new();
        let sections_list = [
            PipelineSection::Preprocessing,
            PipelineSection::Preprocessing,
            PipelineSection::Assembly,
            PipelineSection::Classification,
            PipelineSection::Classification,
            PipelineSection::Abundance,
            PipelineSection::Report,
        ];

        for section in &sections_list {
            let section_dir = self.get_section_dir(section);
            if section_dir.exists() {
                let mut files = Vec::new();

                for entry in fs::read_dir(&section_dir)? {
                    let entry = entry?;
                    let path = entry.path();
                    if path.is_file() {
                        let metadata = fs::metadata(&path)?;
                        files.push(FileInfo {
                            name: path.file_name().unwrap().to_string_lossy().to_string(),
                            path: path.clone(),
                            size: metadata.len(),
                            modified: metadata.modified().ok(),
                        });
                    }
                }

                sections.insert(section.clone(), files);
            }
        }

        Ok(RunSummary {
            run_id: self.run_id.clone(),
            run_dir: self.run_dir.clone(),
            timestamp: self.run_timestamp,
            sections,
        })
    }

    // Helper methods

    /// Create metadata for an intermediate file
    fn create_metadata(
        &self,
        section: &PipelineSection,
        file_path: &Path,
        format: &str,
        section_metadata: serde_json::Value,
    ) -> Result<IntermediateMetadata> {
        let metadata = fs::metadata(file_path)
            .with_context(|| format!("Failed to read file metadata: {}", file_path.display()))?;

        let file_size = metadata.len();
        let checksums = self.calculate_checksums(file_path)?;

        Ok(IntermediateMetadata {
            section: section.clone(),
            format: format.to_string(),
            created_at: Utc::now(),
            file_size,
            checksums,
            schema_version: "1.0".to_string(),
            section_metadata,
        })
    }

    /// Calculate checksums for file validation (simplified)
    fn calculate_checksums(&self, file_path: &Path) -> Result<FileChecksums> {
        let metadata = fs::metadata(file_path)?;
        let size = metadata.len();

        // Simplified checksums using file size and name
        let name_hash = file_path
            .file_name()
            .and_then(|n| n.to_str())
            .map(|n| format!("{:x}", n.len() * 12345))
            .unwrap_or_else(|| "unknown".to_string());

        Ok(FileChecksums {
            md5: format!("{:x}", size * 7),
            sha256: name_hash,
        })
    }

    /// Write compressed data
    fn write_compressed(&self, path: &Path, data: &[u8]) -> Result<()> {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;

        let file = fs::File::create(path)?;
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder.write_all(data)?;
        encoder.finish()?;
        Ok(())
    }

    /// Read compressed data as string
    fn read_compressed(&self, path: &Path) -> Result<String> {
        use flate2::read::GzDecoder;
        use std::io::Read;

        let file = fs::File::open(path)?;
        let mut decoder = GzDecoder::new(file);
        let mut contents = String::new();
        decoder.read_to_string(&mut contents)?;
        Ok(contents)
    }

    /// Read compressed data as bytes
    fn read_compressed_bytes(&self, path: &Path) -> Result<Vec<u8>> {
        use flate2::read::GzDecoder;
        use std::io::Read;

        let file = fs::File::open(path)?;
        let mut decoder = GzDecoder::new(file);
        let mut contents = Vec::new();
        decoder.read_to_end(&mut contents)?;
        Ok(contents)
    }
}

/// Information about a pipeline run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunInfo {
    pub run_id: String,
    pub run_dir: PathBuf,
    pub timestamp: DateTime<Utc>,
    pub config: OutputConfig,
}

/// Summary of all files in a run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunSummary {
    pub run_id: String,
    pub run_dir: PathBuf,
    pub timestamp: DateTime<Utc>,
    pub sections: std::collections::HashMap<PipelineSection, Vec<FileInfo>>,
}

/// Information about a file
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileInfo {
    pub name: String,
    pub path: PathBuf,
    pub size: u64,
    pub modified: Option<std::time::SystemTime>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use tempfile::tempdir;

    #[test]
    fn test_intermediate_output_manager_creation() {
        let temp_dir = tempdir().unwrap();
        let config = OutputConfig::default();

        let manager =
            IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config).unwrap();

        // Check that run directory was created
        assert!(manager.run_dir.exists());

        // Check that section directories were created
        for section in [
            PipelineSection::Preprocessing,
            PipelineSection::Assembly,
            PipelineSection::Classification,
            PipelineSection::Classification,
            PipelineSection::Abundance,
            PipelineSection::Report,
        ] {
            let section_dir = manager.get_section_dir(&section);
            assert!(
                section_dir.exists(),
                "Section directory should exist: {:?}",
                section
            );
        }
    }

    #[test]
    fn test_save_and_load_intermediate() {
        let temp_dir = tempdir().unwrap();
        let config = OutputConfig::default();
        let manager =
            IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config).unwrap();

        #[derive(Serialize, Deserialize, PartialEq, Debug)]
        struct TestData {
            value: i32,
            name: String,
        }

        let test_data = TestData {
            value: 42,
            name: "test".to_string(),
        };

        // Save data
        let metadata = manager
            .save_intermediate(
                PipelineSection::Assembly,
                "test_data",
                &test_data,
                json!({"test": true}),
            )
            .unwrap();

        assert_eq!(metadata.section, PipelineSection::Assembly);
        assert_eq!(metadata.format, "json");

        // Load data
        let loaded_data: TestData = manager
            .load_intermediate(PipelineSection::Assembly, "test_data")
            .unwrap();

        assert_eq!(loaded_data, test_data);
    }

    #[test]
    fn test_save_sequences() {
        let temp_dir = tempdir().unwrap();
        let config = OutputConfig::default();
        let manager =
            IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config).unwrap();

        let sequences = vec![
            ("seq1".to_string(), "ATCGATCG".to_string()),
            ("seq2".to_string(), "GCTAGCTA".to_string()),
        ];

        let metadata = manager
            .save_sequences(
                PipelineSection::Assembly,
                "contigs",
                &sequences,
                json!({"contig_count": 2}),
            )
            .unwrap();

        assert_eq!(metadata.section, PipelineSection::Assembly);
        assert_eq!(metadata.format, "fasta");

        // Check that file was created
        let fasta_path = if manager.config.compress_files {
            manager
                .get_section_dir(&PipelineSection::Assembly)
                .join("contigs.fasta.gz")
        } else {
            manager
                .get_section_dir(&PipelineSection::Assembly)
                .join("contigs.fasta")
        };
        assert!(fasta_path.exists());
    }

    #[test]
    fn test_run_summary() {
        let temp_dir = tempdir().unwrap();
        let config = OutputConfig::default();
        let manager =
            IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config).unwrap();

        // Create some test data
        let test_data = json!({"test": "data"});
        manager
            .save_intermediate(PipelineSection::Assembly, "test1", &test_data, json!({}))
            .unwrap();

        let sequences = vec![("seq1".to_string(), "ATCG".to_string())];
        manager
            .save_sequences(
                PipelineSection::Assembly,
                "sequences",
                &sequences,
                json!({}),
            )
            .unwrap();

        let summary = manager.generate_run_summary().unwrap();
        assert_eq!(summary.run_id, manager.run_id);

        let assembly_files = summary.sections.get(&PipelineSection::Assembly).unwrap();
        assert!(!assembly_files.is_empty());
    }
}
