//! Fast Mode CLI Interface
//! ======================
//!
//! Command-line interface for the high-performance pipeline
//! with optimized file I/O and minimal intermediate outputs.

use anyhow::Result;
use clap::{Parser, Subcommand};
use std::path::PathBuf;
use tracing::{info, level_filters::LevelFilter};
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

use crate::pipeline::fast_pipeline::{FastPipeline, FastPipelineConfig};
use crate::assembly::laptop_assembly::LaptopConfig;

/// Fast mode CLI arguments
#[derive(Parser)]
#[command(name = "meta-forge-fast")]
#[command(about = "High-performance metagenomic analysis pipeline")]
#[command(version = "1.0.0")]
pub struct FastCli {
    #[command(subcommand)]
    pub command: FastCommands,

    /// Output directory
    #[arg(short, long, default_value = "output")]
    pub output: PathBuf,

    /// Enable verbose logging
    #[arg(short, long)]
    pub verbose: bool,

    /// Maximum parallel file operations
    #[arg(long, default_value = "8")]
    pub max_parallel_writes: usize,

    /// Enable minimal output mode (JSON only)
    #[arg(long)]
    pub minimal: bool,

    /// Skip intermediate file saves during processing
    #[arg(long)]
    pub skip_intermediates: bool,
}

/// Fast mode commands
#[derive(Subcommand)]
pub enum FastCommands {
    /// Run the complete pipeline with optimized I/O
    Run {
        /// Input FASTQ files
        #[arg(required = true)]
        input_files: Vec<PathBuf>,

        /// Memory budget in MB for assembly
        #[arg(long, default_value = "2048")]
        memory_mb: usize,

        /// Number of CPU cores to use
        #[arg(long)]
        cpu_cores: Option<usize>,

        /// Assembly timeout in seconds
        #[arg(long, default_value = "300")]
        timeout: u64,
    },

    /// Benchmark file I/O performance
    Benchmark {
        /// Number of test files to create
        #[arg(long, default_value = "100")]
        test_files: usize,

        /// Size of each test file in KB
        #[arg(long, default_value = "100")]
        file_size_kb: usize,
    },

    /// Show performance statistics
    Stats {
        /// Run directory to analyze
        run_dir: Option<PathBuf>,
    },
}

impl FastCli {
    /// Initialize logging
    pub fn init_logging(&self) -> Result<()> {
        let level = if self.verbose {
            LevelFilter::DEBUG
        } else {
            LevelFilter::INFO
        };

        tracing_subscriber::registry()
            .with(fmt::layer().with_target(false))
            .with(
                EnvFilter::builder()
                    .with_default_directive(level.into())
                    .from_env_lossy(),
            )
            .init();

        Ok(())
    }

    /// Execute the CLI command
    pub async fn execute(self) -> Result<()> {
        self.init_logging()?;

        match self.command {
            FastCommands::Run {
                input_files,
                memory_mb,
                cpu_cores,
                timeout,
            } => {
                self.run_pipeline(input_files, memory_mb, cpu_cores, timeout).await
            }
            FastCommands::Benchmark {
                test_files,
                file_size_kb,
            } => {
                self.run_benchmark(test_files, file_size_kb).await
            }
            FastCommands::Stats { run_dir } => {
                self.show_stats(run_dir).await
            }
        }
    }

    /// Run the optimized pipeline
    async fn run_pipeline(
        &self,
        input_files: Vec<PathBuf>,
        memory_mb: usize,
        cpu_cores: Option<usize>,
        timeout: u64,
    ) -> Result<()> {
        info!("🚀 Starting high-performance metagenomic pipeline");
        info!("   📁 Input files: {}", input_files.len());
        info!("   📂 Output directory: {}", self.output.display());
        info!("   💾 Memory budget: {} MB", memory_mb);
        info!("   ⚡ Max parallel writes: {}", self.max_parallel_writes);
        info!("   🎯 Minimal output: {}", self.minimal);
        info!("   ⏭️  Skip intermediates: {}", self.skip_intermediates);

        // Validate input files
        for file in &input_files {
            if !file.exists() {
                return Err(anyhow::anyhow!("Input file not found: {}", file.display()));
            }
        }

        // Configure assembly settings
        let assembly_config = if let Some(cores) = cpu_cores {
            LaptopConfig::custom(memory_mb, cores, 31)?
        } else {
            let mut config = LaptopConfig::auto_detect();
            config.memory_budget_mb = memory_mb;
            config
        };

        // Configure pipeline
        let pipeline_config = FastPipelineConfig {
            output_dir: self.output.clone(),
            minimal_output: self.minimal,
            skip_intermediates: self.skip_intermediates,
            max_parallel_writes: self.max_parallel_writes,
            assembly_config,
            ..FastPipelineConfig::default()
        };

        // Create and run pipeline
        let mut pipeline = FastPipeline::new(pipeline_config).await?;
        let start_time = std::time::Instant::now();

        let report = pipeline.run(input_files).await?;

        let elapsed = start_time.elapsed();
        let stats = pipeline.get_output_stats().await;

        // Display results
        info!("✅ Pipeline completed successfully!");
        info!("   ⏱️  Total time: {:.2}s", elapsed.as_secs_f64());
        info!("   📊 Sample: {}", report.sample_name);
        info!("   📝 Report timestamp: {}", report.analysis_timestamp);
        info!("   💾 Pending writes: {}", stats.pending_operations);
        info!("   🗜️  Compression: {}", stats.compression_enabled);

        if self.verbose {
            info!("📋 Detailed Results:");
            info!("   Input summary: {}", serde_json::to_string_pretty(&report.input_summary)?);
            info!("   Results summary: {}", serde_json::to_string_pretty(&report.results_summary)?);
        }

        Ok(())
    }

    /// Run I/O performance benchmark
    async fn run_benchmark(&self, test_files: usize, file_size_kb: usize) -> Result<()> {
        info!("🏁 Running file I/O benchmark");
        info!("   📄 Test files: {}", test_files);
        info!("   📏 File size: {} KB", file_size_kb);

        let benchmark_dir = self.output.join("benchmark");
        tokio::fs::create_dir_all(&benchmark_dir).await?;

        // Generate test data
        let test_data = vec![b'A'; file_size_kb * 1024];

        // Benchmark synchronous writes
        let sync_start = std::time::Instant::now();
        for i in 0..test_files {
            let path = benchmark_dir.join(format!("sync_test_{}.txt", i));
            std::fs::write(&path, &test_data)?;
        }
        let sync_duration = sync_start.elapsed();

        // Benchmark async writes
        let async_start = std::time::Instant::now();
        let mut async_tasks = Vec::new();

        for i in 0..test_files {
            let path = benchmark_dir.join(format!("async_test_{}.txt", i));
            let data = test_data.clone();
            async_tasks.push(async move {
                tokio::fs::write(&path, &data).await
            });
        }

        // Wait for all async writes
        for task in async_tasks {
            task.await?;
        }
        let async_duration = async_start.elapsed();

        // Display results
        info!("📊 Benchmark Results:");
        info!("   📝 Synchronous writes: {:.2}s ({:.1} MB/s)",
              sync_duration.as_secs_f64(),
              (test_files * file_size_kb) as f64 / 1024.0 / sync_duration.as_secs_f64());
        info!("   ⚡ Asynchronous writes: {:.2}s ({:.1} MB/s)",
              async_duration.as_secs_f64(),
              (test_files * file_size_kb) as f64 / 1024.0 / async_duration.as_secs_f64());
        info!("   🚀 Speedup: {:.2}x", sync_duration.as_secs_f64() / async_duration.as_secs_f64());

        // Cleanup
        tokio::fs::remove_dir_all(&benchmark_dir).await?;

        Ok(())
    }

    /// Show performance statistics
    async fn show_stats(&self, run_dir: Option<PathBuf>) -> Result<()> {
        let target_dir = run_dir.unwrap_or_else(|| self.output.clone());

        info!("📊 Performance Statistics");
        info!("   📂 Directory: {}", target_dir.display());

        if !target_dir.exists() {
            return Err(anyhow::anyhow!("Directory not found: {}", target_dir.display()));
        }

        // Analyze directory structure
        let mut total_files = 0;
        let mut total_size = 0u64;
        let mut file_types = std::collections::HashMap::new();

        let mut entries = tokio::fs::read_dir(&target_dir).await?;
        while let Some(entry) = entries.next_entry().await? {
            let path = entry.path();
            if path.is_file() {
                total_files += 1;
                let metadata = tokio::fs::metadata(&path).await?;
                total_size += metadata.len();

                if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
                    *file_types.entry(ext.to_string()).or_insert(0) += 1;
                }
            }
        }

        info!("   📄 Total files: {}", total_files);
        info!("   💾 Total size: {:.2} MB", total_size as f64 / 1024.0 / 1024.0);

        if !file_types.is_empty() {
            info!("   📋 File types:");
            for (ext, count) in file_types {
                info!("     .{}: {} files", ext, count);
            }
        }

        Ok(())
    }
}

/// Main entry point for fast mode CLI
pub async fn run_fast_cli() -> Result<()> {
    let cli = FastCli::parse();
    cli.execute().await
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[tokio::test]
    async fn test_benchmark() {
        let temp_dir = tempdir().unwrap();
        let cli = FastCli::parse_from([
            "meta-forge-fast",
            "benchmark",
            "--test-files", "10",
            "--file-size-kb", "10"
        ]);

        // This would normally run the benchmark
        // For testing, we just verify the CLI parsing works
        match cli.command {
            FastCommands::Benchmark { test_files, file_size_kb } => {
                assert_eq!(test_files, 10);
                assert_eq!(file_size_kb, 10);
            }
            _ => panic!("Expected benchmark command"),
        }
    }

    #[test]
    fn test_cli_parsing() {
        let cli = FastCli::parse_from([
            "meta-forge-fast",
            "--output", "/tmp/test",
            "--max-parallel-writes", "16",
            "--minimal",
            "run",
            "test1.fastq",
            "test2.fastq"
        ]);

        assert_eq!(cli.output, PathBuf::from("/tmp/test"));
        assert_eq!(cli.max_parallel_writes, 16);
        assert!(cli.minimal);

        match cli.command {
            FastCommands::Run { input_files, .. } => {
                assert_eq!(input_files.len(), 2);
                assert_eq!(input_files[0], PathBuf::from("test1.fastq"));
                assert_eq!(input_files[1], PathBuf::from("test2.fastq"));
            }
            _ => panic!("Expected run command"),
        }
    }
}