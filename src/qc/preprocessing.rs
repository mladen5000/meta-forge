//! Preprocessing module for reading and processing input files
//!
//! Handles:
//! - File format detection (FASTQ, FASTA, compressed)
//! - Reading FASTQ/FASTA files
//! - Applying QC pipeline to raw reads
//! - Generating preprocessing statistics

use anyhow::{Context, Result};
use bio::io::{fasta, fastq};
use colored::Colorize;
use std::path::Path;
use tracing::{debug, info};

use super::{QCPipeline, QCPipelineConfig};
use crate::core::data_structures::{CorrectedRead, CorrectionMetadata};

/// File format types supported by the preprocessor
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fastq,
    FastqGz,
    Fasta,
    FastaGz,
    Unknown,
}

/// Preprocessor for handling raw sequencing data
pub struct Preprocessor {
    qc_config: QCPipelineConfig,
}

impl Preprocessor {
    /// Create a new preprocessor with default QC configuration
    pub fn new() -> Self {
        Self {
            qc_config: QCPipelineConfig::default(),
        }
    }

    /// Create a new preprocessor with custom QC configuration
    pub fn with_config(qc_config: QCPipelineConfig) -> Self {
        Self { qc_config }
    }

    /// Detect file format from extension
    pub fn detect_file_format(&self, path: &Path) -> Result<FileFormat> {
        let extension = path.extension().and_then(|e| e.to_str()).unwrap_or("");

        let format = match extension.to_lowercase().as_str() {
            "fastq" | "fq" => FileFormat::Fastq,
            "gz" if path.to_str().unwrap_or("").contains(".fastq.")
                || path.to_str().unwrap_or("").contains(".fq.") =>
            {
                FileFormat::FastqGz
            }
            "fasta" | "fa" | "fna" => FileFormat::Fasta,
            "gz" if path.to_str().unwrap_or("").contains(".fasta.")
                || path.to_str().unwrap_or("").contains(".fa.") =>
            {
                FileFormat::FastaGz
            }
            _ => FileFormat::Unknown,
        };

        Ok(format)
    }

    /// Process a FASTQ file with QC pipeline
    pub fn process_fastq_file(&self, file_path: &Path) -> Result<Vec<CorrectedRead>> {
        self.process_fastq_file_with_kmer_cache(file_path, None)
    }

    /// Process a FASTQ file with QC pipeline and optionally pre-populate k-mer cache
    ///
    /// # Arguments
    /// * `file_path` - Path to the FASTQ file
    /// * `k_size` - Optional k-mer size to pre-populate cache. If None, cache is not populated.
    pub fn process_fastq_file_with_kmer_cache(&self, file_path: &Path, k_size: Option<usize>) -> Result<Vec<CorrectedRead>> {
        info!("{}", "📖 Reading FASTQ file...".bright_cyan());
        let reader = fastq::Reader::from_file(file_path).context("Failed to open FASTQ file")?;

        // Create QC pipeline
        let mut qc_pipeline = QCPipeline::new(self.qc_config.clone());

        let mut raw_reads = Vec::new();

        // First pass: read all sequences
        for (read_id, record_result) in reader.records().enumerate() {
            let record = record_result?;
            let sequence = std::str::from_utf8(record.seq())?;
            let quality_scores = record.qual().to_vec();

            let read = CorrectedRead {
                id: read_id,
                original: sequence.to_string(),
                corrected: sequence.to_string(),
                corrections: Vec::new(),
                quality_scores,
                correction_metadata: CorrectionMetadata {
                    algorithm: "qc".to_string(),
                    confidence_threshold: 0.0,
                    context_window: 0,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: Vec::new(),
            };

            raw_reads.push(read);

            if read_id % 50000 == 0 && read_id > 0 {
                info!(
                    "  {} {}",
                    "📊 Loaded".bright_blue(),
                    format!("{} reads", read_id).white()
                );
            }
        }

        info!(
            "{} {}",
            "✅ Loaded".bright_green(),
            format!("{} total reads", raw_reads.len()).white().bold()
        );

        // Second pass: QC filtering and trimming
        info!(
            "{}",
            "🔬 Applying quality control filters...".bright_yellow()
        );

        // Debug: Check first few reads
        if let Some(first_read) = raw_reads.first() {
            let avg_qual = first_read
                .quality_scores
                .iter()
                .map(|&q| q as f64)
                .sum::<f64>()
                / first_read.quality_scores.len() as f64;
            debug!(
                "  First read: len={}, avg_qual_raw={:.1}",
                first_read.original.len(),
                avg_qual
            );
        }

        let mut corrected_reads = qc_pipeline.process_reads(&raw_reads);

        // Get and display QC stats
        let qc_stats = qc_pipeline.stats();

        // Print detailed QC statistics with colors
        self.print_qc_stats(&qc_stats);

        if qc_stats.reads_passed == 0 && qc_stats.reads_input > 0 {
            tracing::warn!(
                "{}",
                "⚠️  Warning: All reads filtered out! Check quality settings."
                    .bright_yellow()
                    .bold()
            );
        }

        // Pre-populate k-mer hash cache if k_size is specified
        if let Some(k) = k_size {
            info!("{}", format!("🧬 Pre-populating k-mer cache (k={})...", k).bright_cyan());
            let start = std::time::Instant::now();

            use rayon::prelude::*;
            corrected_reads.par_iter_mut().for_each(|read| {
                read.populate_kmer_hash_cache(k);
            });

            let elapsed = start.elapsed();
            info!(
                "{}",
                format!(
                    "✅ K-mer cache populated for {} reads in {:.2}s",
                    corrected_reads.len(),
                    elapsed.as_secs_f64()
                )
                .bright_green()
            );
        }

        Ok(corrected_reads)
    }

    /// Process a FASTA file (no QC, just load)
    pub fn process_fasta_file(&self, file_path: &Path) -> Result<Vec<CorrectedRead>> {
        self.process_fasta_file_with_kmer_cache(file_path, None)
    }

    /// Process a FASTA file and optionally pre-populate k-mer cache
    ///
    /// # Arguments
    /// * `file_path` - Path to the FASTA file
    /// * `k_size` - Optional k-mer size to pre-populate cache. If None, cache is not populated.
    pub fn process_fasta_file_with_kmer_cache(&self, file_path: &Path, k_size: Option<usize>) -> Result<Vec<CorrectedRead>> {
        info!("{}", "📖 Reading FASTA file...".bright_cyan());
        let reader = fasta::Reader::from_file(file_path).context("Failed to open FASTA file")?;

        let mut corrected_reads = Vec::new();

        for (read_id, record_result) in reader.records().enumerate() {
            let record = record_result?;
            let sequence = std::str::from_utf8(record.seq())?;

            let corrected_read = CorrectedRead {
                id: read_id,
                original: sequence.to_string(),
                corrected: sequence.to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30; sequence.len()], // Default quality
                correction_metadata: CorrectionMetadata {
                    algorithm: "none".to_string(),
                    confidence_threshold: 0.0,
                    context_window: 0,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: Vec::new(),
            };

            corrected_reads.push(corrected_read);

            if read_id % 50000 == 0 && read_id > 0 {
                info!(
                    "  {} {}",
                    "📊 Loaded".bright_blue(),
                    format!("{} reads", read_id).white()
                );
            }
        }

        info!(
            "{} {}",
            "✅ Loaded".bright_green(),
            format!("{} total reads", corrected_reads.len())
                .white()
                .bold()
        );

        // Pre-populate k-mer hash cache if k_size is specified
        if let Some(k) = k_size {
            info!("{}", format!("🧬 Pre-populating k-mer cache (k={})...", k).bright_cyan());
            let start = std::time::Instant::now();

            use rayon::prelude::*;
            corrected_reads.par_iter_mut().for_each(|read| {
                read.populate_kmer_hash_cache(k);
            });

            let elapsed = start.elapsed();
            info!(
                "{}",
                format!(
                    "✅ K-mer cache populated for {} reads in {:.2}s",
                    corrected_reads.len(),
                    elapsed.as_secs_f64()
                )
                .bright_green()
            );
        }

        Ok(corrected_reads)
    }

    /// Print QC statistics with colored output
    fn print_qc_stats(&self, qc_stats: &crate::qc::QCStats) {
        use colored::Colorize;

        println!(
            "\n{}",
            "═══════════════════════════════════════════".bright_cyan()
        );
        println!("{}", "   QUALITY CONTROL STATISTICS".bright_cyan().bold());
        println!(
            "{}",
            "═══════════════════════════════════════════".bright_cyan()
        );

        let pass_rate = if qc_stats.reads_input > 0 {
            (qc_stats.reads_passed as f64 / qc_stats.reads_input as f64) * 100.0
        } else {
            0.0
        };

        println!("\n{}", "📊 Input/Output Summary:".bright_blue().bold());
        println!(
            "  {} {}",
            "Input reads:".white(),
            format!("{:>10}", qc_stats.reads_input).yellow()
        );
        println!(
            "  {} {}",
            "Passed:     ".white(),
            format!("{:>10} ({:.1}%)", qc_stats.reads_passed, pass_rate)
                .bright_green()
                .bold()
        );
        println!(
            "  {} {}",
            "Failed:     ".white(),
            format!("{:>10} ({:.1}%)", qc_stats.reads_failed, 100.0 - pass_rate).bright_red()
        );

        if qc_stats.reads_failed > 0 {
            println!("\n{}", "❌ Failure Reasons:".bright_red().bold());
            println!(
                "  {} {:>10}",
                "Quality:".white(),
                qc_stats.reads_failed_quality.to_string().red()
            );
            println!(
                "  {} {:>10}",
                "Length: ".white(),
                qc_stats.reads_failed_length.to_string().red()
            );
        }

        if qc_stats.adapters_detected > 0 {
            println!("\n{}", "✂️  Adapter Trimming:".bright_magenta().bold());
            println!(
                "  {} {:>10}",
                "Adapters detected:".white(),
                qc_stats.adapters_detected.to_string().magenta()
            );
            println!(
                "  {} {:>10}",
                "Bases removed:   ".white(),
                qc_stats.bases_trimmed_adapter.to_string().magenta()
            );
        }

        println!("\n{}", "📏 Length Statistics:".bright_blue().bold());
        println!(
            "  {} {:>10.1} bp",
            "Before:".white(),
            qc_stats.mean_length_before
        );
        println!(
            "  {} {:>10.1} bp",
            "After: ".white(),
            qc_stats.mean_length_after
        );

        println!("\n{}", "⭐ Quality Metrics:".bright_blue().bold());
        println!(
            "  {} Q{:.1}  →  Q{:.1}",
            "Average:".white(),
            qc_stats.mean_quality_before,
            qc_stats.mean_quality_after
        );
        println!(
            "  {} {:.1}%  →  {:.1}%",
            "Q20:    ".white(),
            qc_stats.q20_percentage_before,
            qc_stats.q20_percentage_after
        );
        println!(
            "  {} {:.1}%  →  {:.1}%",
            "Q30:    ".white(),
            qc_stats.q30_percentage_before,
            qc_stats.q30_percentage_after
        );

        println!(
            "\n{}",
            "═══════════════════════════════════════════\n".bright_cyan()
        );
    }
}

impl Default for Preprocessor {
    fn default() -> Self {
        Self::new()
    }
}
