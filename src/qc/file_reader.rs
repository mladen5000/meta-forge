//! QC File Reader Module
//!
//! Handles preprocessing and reading of input sequencing files.
//! Supports FASTQ and FASTA formats (both plain and gzipped).
//!
//! # Features
//! - Automatic file format detection
//! - FASTQ quality filtering with QC pipeline
//! - FASTA parsing with default quality scores
//! - Detailed QC statistics and logging
//! - Intermediate output checkpoint management

use anyhow::{Context, Result};
use colored::Colorize;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tracing::{debug, info, warn};

use crate::core::{CorrectedRead, CorrectionMetadata, FileFormat};
use crate::qc::{QCPipeline, QCPipelineConfig};
use crate::qc::qc_stats::QCStats;
use crate::utils::intermediate_output::{IntermediateOutputManager, PipelineSection};

/// File reader for preprocessing input sequencing data
///
/// Handles file format detection, QC filtering, and read processing.
pub struct FileReader {
    qc_pipeline: Arc<QCPipeline>,
    output_manager: IntermediateOutputManager,
}

impl FileReader {
    /// Create a new FileReader with QC pipeline and output manager
    pub fn new(
        qc_pipeline: Arc<QCPipeline>,
        output_manager: IntermediateOutputManager,
    ) -> Self {
        Self {
            qc_pipeline,
            output_manager,
        }
    }

    /// Preprocess input files with error correction
    pub async fn preprocess_inputs(
        &self,
        inputs: &[PathBuf],
    ) -> Result<(Vec<CorrectedRead>, QCStats)> {
        info!("📋 Starting preprocessing of {} input files", inputs.len());
        let mut all_reads = Vec::new();
        let mut combined_stats = QCStats::default();

        // Save preprocessing initialization status immediately
        let preprocessing_init = serde_json::json!({
            "status": "started",
            "input_files": inputs.len(),
            "file_paths": inputs.iter().map(|p| p.to_string_lossy()).collect::<Vec<_>>(),
            "timestamp": chrono::Utc::now().to_rfc3339()
        });

        self.output_manager.save_intermediate(
            PipelineSection::Preprocessing,
            "initialization_status",
            &preprocessing_init,
            serde_json::json!({"phase": "initialization"}),
        )?;
        info!("💾 Saved preprocessing initialization status");

        for (file_idx, input_file) in inputs.iter().enumerate() {
            info!(
                "📖 Processing input file {}/{}: {}",
                file_idx + 1,
                inputs.len(),
                input_file.display()
            );

            // Save per-file processing status
            let file_status = serde_json::json!({
                "file_path": input_file.to_string_lossy(),
                "file_index": file_idx,
                "status": "processing",
                "timestamp": chrono::Utc::now().to_rfc3339()
            });

            self.output_manager.save_intermediate(
                PipelineSection::Preprocessing,
                &format!("file_{file_idx}_processing"),
                &file_status,
                serde_json::json!({"file_processing": true}),
            )?;

            // Determine file format
            let format = self.detect_file_format(input_file)?;

            // Read and process based on format
            let (reads, stats) = match format {
                FileFormat::Fastq | FileFormat::FastqGz => {
                    self.process_fastq_file(input_file).await?
                }
                FileFormat::Fasta | FileFormat::FastaGz => {
                    let reads = self.process_fasta_file(input_file).await?;
                    (reads, QCStats::default())
                }
                FileFormat::Unknown => {
                    return Err(anyhow::anyhow!(
                        "Unsupported file format: {}",
                        input_file.display()
                    ));
                }
            };

            // Aggregate stats (same as in preprocess_inputs_with_progress)
            combined_stats.reads_input += stats.reads_input;
            combined_stats.reads_passed += stats.reads_passed;
            combined_stats.reads_failed += stats.reads_failed;
            combined_stats.reads_failed_quality += stats.reads_failed_quality;
            combined_stats.reads_failed_length += stats.reads_failed_length;
            combined_stats.reads_failed_adapter += stats.reads_failed_adapter;
            combined_stats.bases_trimmed_quality += stats.bases_trimmed_quality;
            combined_stats.bases_trimmed_adapter += stats.bases_trimmed_adapter;
            combined_stats.total_bases_before += stats.total_bases_before;
            combined_stats.total_bases_after += stats.total_bases_after;
            combined_stats.adapters_detected += stats.adapters_detected;

            // Merge adapter types
            for (adapter, count) in stats.adapter_types {
                *combined_stats.adapter_types.entry(adapter).or_insert(0) += count;
            }

            // Accumulate quality/length sums (will average later)
            combined_stats.mean_quality_before +=
                stats.mean_quality_before * stats.reads_input as f64;
            combined_stats.mean_quality_after +=
                stats.mean_quality_after * stats.reads_passed as f64;
            combined_stats.mean_length_before +=
                stats.mean_length_before * stats.reads_input as f64;
            combined_stats.mean_length_after += stats.mean_length_after * stats.reads_passed as f64;
            combined_stats.q20_percentage_before +=
                stats.q20_percentage_before * stats.reads_input as f64;
            combined_stats.q30_percentage_before +=
                stats.q30_percentage_before * stats.reads_input as f64;
            combined_stats.q20_percentage_after +=
                stats.q20_percentage_after * stats.reads_passed as f64;
            combined_stats.q30_percentage_after +=
                stats.q30_percentage_after * stats.reads_passed as f64;

            all_reads.extend(reads.clone());

            // Save completed file processing status
            let file_complete = serde_json::json!({
                "file_path": input_file.to_string_lossy(),
                "file_index": file_idx,
                "status": "completed",
                "reads_processed": reads.len(),
                "total_reads_so_far": all_reads.len(),
                "timestamp": chrono::Utc::now().to_rfc3339()
            });

            self.output_manager.save_intermediate(
                PipelineSection::Preprocessing,
                &format!("file_{file_idx}_completed"),
                &file_complete,
                serde_json::json!({"file_completed": true}),
            )?;

            info!(
                "✅ Completed processing file {}: {} reads",
                input_file.display(),
                reads.len()
            );
        }

        // Save preprocessing completion summary
        let preprocessing_summary = serde_json::json!({
            "status": "completed",
            "total_files_processed": inputs.len(),
            "total_reads_processed": all_reads.len(),
            "average_read_length": all_reads.iter().map(|r| r.corrected.len()).sum::<usize>() as f64 / all_reads.len() as f64,
            "total_corrections": all_reads.iter().map(|r| r.corrections.len()).sum::<usize>(),
            "timestamp": chrono::Utc::now().to_rfc3339()
        });

        self.output_manager.save_intermediate(
            PipelineSection::Preprocessing,
            "final_summary",
            &preprocessing_summary,
            serde_json::json!({"final_summary": true}),
        )?;
        info!("💾 Saved preprocessing final summary");

        // Calculate final averages
        if combined_stats.reads_input > 0 {
            combined_stats.mean_quality_before /= combined_stats.reads_input as f64;
            combined_stats.mean_length_before /= combined_stats.reads_input as f64;
            combined_stats.q20_percentage_before /= combined_stats.reads_input as f64;
            combined_stats.q30_percentage_before /= combined_stats.reads_input as f64;
        }

        if combined_stats.reads_passed > 0 {
            combined_stats.mean_quality_after /= combined_stats.reads_passed as f64;
            combined_stats.mean_length_after /= combined_stats.reads_passed as f64;
            combined_stats.q20_percentage_after /= combined_stats.reads_passed as f64;
            combined_stats.q30_percentage_after /= combined_stats.reads_passed as f64;
        }

        info!(
            "📊 Preprocessing completed: {} reads from {} files",
            all_reads.len(),
            inputs.len()
        );
        Ok((all_reads, combined_stats))
    }

    /// Process FASTQ file with quality filtering
    async fn process_fastq_file(
        &self,
        file_path: &Path,
    ) -> Result<(Vec<CorrectedRead>, QCStats)> {
        use bio::io::fastq;

        info!("{}", "📖 Reading FASTQ file...".bright_cyan());
        let reader = fastq::Reader::from_file(file_path).context("Failed to open FASTQ file")?;

        // Create QC pipeline
        let qc_config = QCPipelineConfig::default();
        let mut qc_pipeline = QCPipeline::new(qc_config);

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

        let corrected_reads = qc_pipeline.process_reads(&raw_reads);

        // Get and display QC stats
        let qc_stats = qc_pipeline.stats();

        // Print detailed QC statistics with colors
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

        if qc_stats.reads_passed == 0 && qc_stats.reads_input > 0 {
            warn!(
                "{}",
                "⚠️  Warning: All reads filtered out! Check quality settings."
                    .bright_yellow()
                    .bold()
            );
        }

        Ok((corrected_reads, qc_stats))
    }

    /// Process FASTA file (no quality filtering)
    async fn process_fasta_file(&self, file_path: &Path) -> Result<Vec<CorrectedRead>> {
        use bio::io::fasta;

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
        }

        Ok(corrected_reads)
    }

    /// Detect file format from extension
    fn detect_file_format(&self, file_path: &Path) -> Result<FileFormat> {
        // First check if this looks like a CLI flag (starts with -)
        let path_str = file_path.to_string_lossy();
        if path_str.starts_with('-') {
            return Err(anyhow::anyhow!(
                "Invalid input: '{}' looks like a command-line flag. \
                Global flags like -m, -j, -o must come BEFORE the subcommand.\n\
                Correct: meta-forge -m 4096 -j 8 analyze input.fastq\n\
                Wrong: meta-forge analyze input.fastq -m 4096",
                path_str
            ));
        }

        let extension = file_path
            .extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("");

        match extension {
            "fastq" | "fq" => Ok(FileFormat::Fastq),
            "gz" => {
                // Check the extension before .gz
                let stem = file_path
                    .file_stem()
                    .and_then(|stem| Path::new(stem).extension())
                    .and_then(|ext| ext.to_str())
                    .unwrap_or("");

                match stem {
                    "fastq" | "fq" => Ok(FileFormat::FastqGz),
                    "fasta" | "fa" | "fas" => Ok(FileFormat::FastaGz),
                    _ => Ok(FileFormat::Unknown),
                }
            }
            "fasta" | "fa" | "fas" => Ok(FileFormat::Fasta),
            _ => Ok(FileFormat::Unknown),
        }
    }
}
