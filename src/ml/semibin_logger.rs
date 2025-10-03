//! SemiBin2 execution logging and progress tracking
//!
//! This module provides detailed monitoring and logging for SemiBin2 integration,
//! including subprocess output parsing, progress tracking, and result validation.

use anyhow::{Context, Result};
use std::io::{BufRead, BufReader};
use std::process::{Child, Command, Stdio};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::thread;
use tracing::{debug, info, warn, error};

use crate::ml::classification_reporter::{ClassificationReporter, ClassificationStage};

/// SemiBin2 execution configuration
#[derive(Debug, Clone)]
pub struct SemiBinConfig {
    /// Path to semibin2 executable
    pub executable: PathBuf,
    /// Input contigs FASTA file
    pub contigs_file: PathBuf,
    /// BAM alignment files
    pub bam_files: Vec<PathBuf>,
    /// Output directory
    pub output_dir: PathBuf,
    /// Number of threads to use
    pub threads: usize,
    /// Minimum contig length
    pub min_contig_length: usize,
    /// Environment type (e.g., "human_gut", "ocean", "soil")
    pub environment: Option<String>,
    /// Additional command-line arguments
    pub extra_args: Vec<String>,
}

/// SemiBin2 execution result
#[derive(Debug, Clone)]
pub struct SemiBinResult {
    /// Exit code
    pub exit_code: i32,
    /// Number of bins created
    pub num_bins: usize,
    /// Output bin directory
    pub bin_dir: PathBuf,
    /// Path to binning result file
    pub result_file: PathBuf,
    /// Execution time in seconds
    pub execution_time_seconds: f64,
    /// Total contigs processed
    pub contigs_processed: usize,
    /// Total contigs binned
    pub contigs_binned: usize,
}

/// SemiBin2 progress event
#[derive(Debug, Clone)]
pub enum SemiBinEvent {
    /// Process started
    Started,
    /// Reading input files
    ReadingInput { file_type: String },
    /// Feature extraction progress
    FeatureExtraction { processed: usize, total: usize },
    /// Clustering/binning stage
    Clustering { iteration: usize },
    /// Writing output
    WritingOutput { bin_count: usize },
    /// Process completed
    Completed { bins_created: usize },
    /// Error occurred
    Error { message: String },
    /// Generic log message
    Log { level: LogLevel, message: String },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogLevel {
    Debug,
    Info,
    Warning,
    Error,
}

/// SemiBin2 executor with progress tracking
pub struct SemiBinExecutor {
    config: SemiBinConfig,
    reporter: Option<Arc<Mutex<ClassificationReporter>>>,
}

impl SemiBinExecutor {
    /// Create a new SemiBin executor
    pub fn new(config: SemiBinConfig) -> Self {
        Self {
            config,
            reporter: None,
        }
    }

    /// Attach a classification reporter for integrated logging
    pub fn with_reporter(mut self, reporter: Arc<Mutex<ClassificationReporter>>) -> Self {
        self.reporter = Some(reporter);
        self
    }

    /// Execute SemiBin2 and track progress
    pub fn execute(&self) -> Result<SemiBinResult> {
        info!("🔬 Starting SemiBin2 execution");
        info!("   Executable: {}", self.config.executable.display());
        info!("   Contigs: {}", self.config.contigs_file.display());
        info!("   BAM files: {}", self.config.bam_files.len());
        info!("   Output: {}", self.config.output_dir.display());
        info!("   Threads: {}", self.config.threads);
        info!("   Min contig length: {}", self.config.min_contig_length);

        if let Some(env) = &self.config.environment {
            info!("   Environment: {}", env);
        }

        debug!("Configuration details:");
        for (i, bam) in self.config.bam_files.iter().enumerate() {
            debug!("   BAM[{}]: {}", i, bam.display());
        }
        if !self.config.extra_args.is_empty() {
            debug!("   Extra args: {:?}", self.config.extra_args);
        }

        let start_time = std::time::Instant::now();

        // Update reporter if available
        if let Some(reporter) = &self.reporter {
            if let Ok(mut r) = reporter.lock() {
                info!("   Reporter attached - tracking progress");
                r.start_stage(ClassificationStage::Clustering);
            }
        } else {
            debug!("   No reporter attached - minimal tracking");
        }

        // Build command
        let mut cmd = Command::new(&self.config.executable);
        cmd.arg("single_easy_bin")
            .arg("-i")
            .arg(&self.config.contigs_file)
            .arg("-o")
            .arg(&self.config.output_dir)
            .arg("-t")
            .arg(self.config.threads.to_string())
            .arg("--min-len")
            .arg(self.config.min_contig_length.to_string());

        // Add BAM files
        for bam_file in &self.config.bam_files {
            cmd.arg("-b").arg(bam_file);
        }

        // Add environment if specified
        if let Some(env) = &self.config.environment {
            cmd.arg("--environment").arg(env);
        }

        // Add extra arguments
        for arg in &self.config.extra_args {
            cmd.arg(arg);
        }

        // Set up stdout/stderr capture
        cmd.stdout(Stdio::piped()).stderr(Stdio::piped());

        info!("   Launching SemiBin2 subprocess...");
        debug!("   Command: {:?}", cmd);

        // Spawn process
        let mut child = cmd
            .spawn()
            .context("Failed to spawn SemiBin2 process")?;

        // Capture and parse output in separate threads
        let stdout_handle = child.stdout.take();
        let stderr_handle = child.stderr.take();

        let stdout_thread = if let Some(stdout) = stdout_handle {
            let reporter = self.reporter.clone();
            Some(thread::spawn(move || {
                Self::process_stdout(stdout, reporter);
            }))
        } else {
            None
        };

        let stderr_thread = if let Some(stderr) = stderr_handle {
            let reporter = self.reporter.clone();
            Some(thread::spawn(move || {
                Self::process_stderr(stderr, reporter);
            }))
        } else {
            None
        };

        // Wait for process completion
        debug!("   Waiting for SemiBin2 process to complete...");
        let exit_status = child
            .wait()
            .context("Failed to wait for SemiBin2 process")?;

        // Wait for output threads
        debug!("   Collecting output from threads...");
        if let Some(handle) = stdout_thread {
            let _ = handle.join();
        }
        if let Some(handle) = stderr_thread {
            let _ = handle.join();
        }

        let execution_time = start_time.elapsed().as_secs_f64();
        let exit_code = exit_status.code().unwrap_or(-1);

        info!("   SemiBin2 completed with exit code: {}", exit_code);
        info!("   Execution time: {:.2}s", execution_time);

        if !exit_status.success() {
            error!("❌ SemiBin2 execution failed with exit code: {}", exit_code);
            error!("   Check stderr output above for details");
            return Err(anyhow::anyhow!(
                "SemiBin2 failed with exit code: {}",
                exit_code
            ));
        }

        // Parse results
        info!("   Parsing SemiBin2 output files...");
        let result = self.parse_results(execution_time, exit_code)?;

        // Update reporter
        if let Some(reporter) = &self.reporter {
            if let Ok(mut r) = reporter.lock() {
                r.complete_stage(ClassificationStage::Clustering, result.contigs_processed);
            }
        }

        info!("✅ SemiBin2 execution completed successfully");
        info!("   Bins created: {}", result.num_bins);
        info!("   Contigs binned: {}/{}", result.contigs_binned, result.contigs_processed);

        Ok(result)
    }

    /// Process stdout and emit events
    fn process_stdout(
        stdout: impl std::io::Read,
        reporter: Option<Arc<Mutex<ClassificationReporter>>>,
    ) {
        let reader = BufReader::new(stdout);

        for line in reader.lines() {
            if let Ok(line) = line {
                let trimmed = line.trim();

                // Parse SemiBin2 output patterns
                if trimmed.contains("Reading") {
                    debug!("📖 {}", trimmed);
                } else if trimmed.contains("Extracting") || trimmed.contains("features") {
                    debug!("🔍 {}", trimmed);

                    // Try to parse progress
                    if let Some(progress) = Self::parse_progress(&trimmed) {
                        if let Some(ref reporter) = reporter {
                            if let Ok(r) = reporter.lock() {
                                r.log_feature_extraction_progress(progress.0, progress.1);
                            }
                        }
                    }
                } else if trimmed.contains("Clustering") || trimmed.contains("iteration") {
                    debug!("🎯 {}", trimmed);
                } else if trimmed.contains("Writing") || trimmed.contains("bins") {
                    info!("💾 {}", trimmed);
                } else if !trimmed.is_empty() {
                    debug!("   {}", trimmed);
                }
            }
        }
    }

    /// Process stderr and emit warnings/errors
    fn process_stderr(
        stderr: impl std::io::Read,
        _reporter: Option<Arc<Mutex<ClassificationReporter>>>,
    ) {
        let reader = BufReader::new(stderr);

        for line in reader.lines() {
            if let Ok(line) = line {
                let trimmed = line.trim();
                if !trimmed.is_empty() {
                    // Check severity
                    if trimmed.to_lowercase().contains("error") {
                        error!("❌ SemiBin2: {}", trimmed);
                    } else if trimmed.to_lowercase().contains("warning") {
                        warn!("⚠️  SemiBin2: {}", trimmed);
                    } else {
                        debug!("   SemiBin2: {}", trimmed);
                    }
                }
            }
        }
    }

    /// Parse progress information from output line
    fn parse_progress(line: &str) -> Option<(usize, usize)> {
        // Look for patterns like "Processing 500/1000" or "500 of 1000"
        let parts: Vec<&str> = line.split_whitespace().collect();

        for (i, part) in parts.iter().enumerate() {
            if part.contains('/') {
                // Format: "500/1000"
                let nums: Vec<&str> = part.split('/').collect();
                if nums.len() == 2 {
                    if let (Ok(current), Ok(total)) = (nums[0].parse(), nums[1].parse()) {
                        return Some((current, total));
                    }
                }
            } else if *part == "of" && i > 0 && i < parts.len() - 1 {
                // Format: "500 of 1000"
                if let (Ok(current), Ok(total)) = (parts[i - 1].parse(), parts[i + 1].parse()) {
                    return Some((current, total));
                }
            }
        }

        None
    }

    /// Parse SemiBin2 output files to get results
    fn parse_results(&self, execution_time: f64, exit_code: i32) -> Result<SemiBinResult> {
        debug!("   Locating output bin directory...");
        let bin_dir = self.config.output_dir.join("output_bins");

        // Try multiple possible result file locations
        debug!("   Searching for result files...");
        let result_file = if self.config.output_dir.join("output_recluster_bins.tsv").exists() {
            debug!("   Found: output_recluster_bins.tsv");
            self.config.output_dir.join("output_recluster_bins.tsv")
        } else if self.config.output_dir.join("output_bins.tsv").exists() {
            debug!("   Found: output_bins.tsv");
            self.config.output_dir.join("output_bins.tsv")
        } else {
            debug!("   Using fallback: binning_results.tsv");
            self.config.output_dir.join("binning_results.tsv")
        };

        info!("   Result file: {}", result_file.display());

        // Count bins
        debug!("   Counting bins in: {}", bin_dir.display());
        let num_bins = if bin_dir.exists() {
            let count = std::fs::read_dir(&bin_dir)
                .map(|entries| entries.filter_map(Result::ok).count())
                .unwrap_or(0);
            debug!("   Found {} bins", count);
            count
        } else {
            warn!("   Bin directory not found: {}", bin_dir.display());
            0
        };

        // Parse result file for contig counts
        debug!("   Parsing result file for contig statistics...");
        let (contigs_processed, contigs_binned) = if result_file.exists() {
            Self::parse_result_file(&result_file)?
        } else {
            warn!("   Result file not found: {}", result_file.display());
            (0, 0)
        };

        Ok(SemiBinResult {
            exit_code,
            num_bins,
            bin_dir,
            result_file,
            execution_time_seconds: execution_time,
            contigs_processed,
            contigs_binned,
        })
    }

    /// Parse TSV result file to count contigs
    fn parse_result_file(path: &PathBuf) -> Result<(usize, usize)> {
        use std::collections::HashSet;
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        debug!("   Reading result file: {}", path.display());
        let file = File::open(path)
            .context(format!("Failed to open result file: {}", path.display()))?;
        let reader = BufReader::new(file);

        let mut all_contigs = HashSet::new();
        let mut binned_contigs = HashSet::new();
        let mut line_count = 0;

        for (idx, line) in reader.lines().enumerate() {
            if idx == 0 {
                debug!("   Skipping header line");
                continue; // Skip header
            }

            if let Ok(line) = line {
                line_count += 1;
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() >= 2 {
                    let contig_id = parts[0];
                    let bin_id = parts[1];

                    all_contigs.insert(contig_id.to_string());

                    if bin_id != "unbinned" && !bin_id.is_empty() {
                        binned_contigs.insert(contig_id.to_string());
                    }
                } else if !line.trim().is_empty() {
                    warn!("   Malformed line {} in result file", idx + 1);
                }
            }
        }

        let total = all_contigs.len();
        let binned = binned_contigs.len();
        debug!("   Parsed {} lines from result file", line_count);
        debug!("   Total contigs: {}, Binned: {}", total, binned);

        if total == 0 {
            warn!("   No contigs found in result file");
        }

        Ok((total, binned))
    }
}

/// Validate SemiBin2 installation and get version
pub fn validate_semibin_installation(executable: &PathBuf) -> Result<String> {
    info!("🔍 Validating SemiBin2 installation");

    let output = Command::new(executable)
        .arg("--version")
        .output()
        .context("Failed to execute semibin2 --version")?;

    if !output.status.success() {
        return Err(anyhow::anyhow!(
            "SemiBin2 not found or not executable: {}",
            executable.display()
        ));
    }

    let version = String::from_utf8_lossy(&output.stdout)
        .trim()
        .to_string();

    info!("✅ SemiBin2 found: {}", version);
    Ok(version)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_progress() {
        assert_eq!(
            SemiBinExecutor::parse_progress("Processing 500/1000 contigs"),
            Some((500, 1000))
        );

        assert_eq!(
            SemiBinExecutor::parse_progress("Extracted features for 250 of 500 contigs"),
            Some((250, 500))
        );

        assert_eq!(
            SemiBinExecutor::parse_progress("No progress here"),
            None
        );
    }

    #[test]
    fn test_log_level() {
        assert_eq!(LogLevel::Info, LogLevel::Info);
        assert_ne!(LogLevel::Error, LogLevel::Warning);
    }
}
