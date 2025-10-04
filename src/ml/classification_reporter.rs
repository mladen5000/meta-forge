//! Comprehensive classification progress and quality reporting
//!
//! This module provides detailed logging, metrics, and intermediate outputs
//! for the classification/binning stage of the metagenomics pipeline.

use anyhow::{Context, Result};
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use tracing::{debug, info};

use crate::utils::intermediate_output::{IntermediateOutputManager, PipelineSection};

/// Detailed classification progress tracker
#[derive(Debug, Clone)]
pub struct ClassificationReporter {
    /// Output manager for saving results
    output_manager: IntermediateOutputManager,
    /// Current classification run metadata
    run_metadata: ClassificationRunMetadata,
    /// Timing information for each stage
    stage_timings: HashMap<ClassificationStage, StageTiming>,
    /// Quality metrics for bins
    bin_metrics: Vec<BinQualityMetrics>,
}

/// Metadata about the classification run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassificationRunMetadata {
    /// Total number of contigs to classify
    pub total_contigs: usize,
    /// Number of contigs passing minimum length filter
    pub valid_contigs: usize,
    /// Classification method used
    pub method: ClassificationMethod,
    /// Configuration parameters
    pub config: serde_json::Value,
    /// Start timestamp
    pub started_at: DateTime<Utc>,
    /// End timestamp (if completed)
    pub completed_at: Option<DateTime<Utc>>,
}

/// Classification method types
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum ClassificationMethod {
    /// ML-based k-mer classification
    MLKmer,
    /// External SemiBin2 integration
    SemiBin2,
    /// Hybrid approach
    Hybrid,
}

/// Classification pipeline stages
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ClassificationStage {
    /// Initial contig filtering
    Filtering,
    /// Feature extraction (k-mer frequencies, coverage)
    FeatureExtraction,
    /// Feature normalization
    Normalization,
    /// Clustering/binning algorithm
    Clustering,
    /// Bin quality assessment
    QualityAssessment,
    /// Writing output files
    OutputWriting,
}

impl ClassificationStage {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Filtering => "Filtering",
            Self::FeatureExtraction => "Feature Extraction",
            Self::Normalization => "Normalization",
            Self::Clustering => "Clustering",
            Self::QualityAssessment => "Quality Assessment",
            Self::OutputWriting => "Output Writing",
        }
    }

    pub fn description(&self) -> &'static str {
        match self {
            Self::Filtering => "Filtering contigs by minimum length and quality criteria",
            Self::FeatureExtraction => "Extracting k-mer frequencies and coverage features",
            Self::Normalization => "Normalizing feature vectors (Z-score or Min-Max)",
            Self::Clustering => "Performing k-means clustering to create bins",
            Self::QualityAssessment => "Computing bin quality metrics (completeness, contamination)",
            Self::OutputWriting => "Writing binning results to output files",
        }
    }
}

/// Timing information for a pipeline stage
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageTiming {
    pub stage: ClassificationStage,
    pub started_at: DateTime<Utc>,
    pub completed_at: Option<DateTime<Utc>>,
    pub duration_seconds: Option<f64>,
    pub items_processed: usize,
    pub throughput_per_second: Option<f64>,
}

/// Quality metrics for a single bin
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinQualityMetrics {
    /// Bin identifier
    pub bin_id: usize,
    /// Number of contigs in bin
    pub contig_count: usize,
    /// Total base pairs in bin
    pub total_bp: usize,
    /// Average contig length
    pub avg_contig_length: f64,
    /// Median contig length
    pub median_contig_length: usize,
    /// Average coverage of contigs in bin
    pub avg_coverage: f64,
    /// Coverage standard deviation
    pub coverage_std: f64,
    /// Average GC content
    pub avg_gc_content: f64,
    /// GC content standard deviation
    pub gc_std: f64,
    /// Estimated completeness (0-100%)
    pub completeness: Option<f64>,
    /// Estimated contamination (0-100%)
    pub contamination: Option<f64>,
    /// Average confidence score for bin assignments
    pub avg_confidence: f64,
    /// N50 statistic
    pub n50: usize,
    /// Longest contig in bin
    pub longest_contig: usize,
}

/// Classification result with detailed statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassificationResults {
    /// Run metadata
    pub metadata: ClassificationRunMetadata,
    /// Per-bin quality metrics
    pub bin_metrics: Vec<BinQualityMetrics>,
    /// Stage timings
    pub stage_timings: Vec<StageTiming>,
    /// Overall statistics
    pub overall_stats: OverallStatistics,
}

/// Overall classification statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverallStatistics {
    /// Total contigs classified
    pub total_classified: usize,
    /// Total contigs filtered out
    pub total_filtered: usize,
    /// Number of bins created
    pub num_bins: usize,
    /// Total base pairs classified
    pub total_bp_classified: usize,
    /// Average bin size (contigs)
    pub avg_bin_size: f64,
    /// Average bin size (bp)
    pub avg_bin_size_bp: f64,
    /// Total execution time (seconds)
    pub total_time_seconds: f64,
    /// Overall throughput (contigs/second)
    pub throughput_contigs_per_sec: f64,
}

impl ClassificationReporter {
    /// Create a new classification reporter
    pub fn new(
        output_manager: IntermediateOutputManager,
        total_contigs: usize,
        method: ClassificationMethod,
        config: serde_json::Value,
    ) -> Self {
        let run_metadata = ClassificationRunMetadata {
            total_contigs,
            valid_contigs: 0,
            method,
            config,
            started_at: Utc::now(),
            completed_at: None,
        };

        info!("🔬 Classification Run Started");
        info!("   Total contigs: {}", total_contigs);
        info!("   Method: {:?}", run_metadata.method);
        info!("   Started at: {}", run_metadata.started_at.format("%Y-%m-%d %H:%M:%S"));

        Self {
            output_manager,
            run_metadata,
            stage_timings: HashMap::new(),
            bin_metrics: Vec::new(),
        }
    }

    /// Start tracking a pipeline stage
    pub fn start_stage(&mut self, stage: ClassificationStage) {
        let timing = StageTiming {
            stage,
            started_at: Utc::now(),
            completed_at: None,
            duration_seconds: None,
            items_processed: 0,
            throughput_per_second: None,
        };

        info!("▶️  {} - {}", stage.as_str(), stage.description());
        self.stage_timings.insert(stage, timing);
    }

    /// Complete a pipeline stage with item count
    pub fn complete_stage(&mut self, stage: ClassificationStage, items_processed: usize) {
        if let Some(timing) = self.stage_timings.get_mut(&stage) {
            let completed_at = Utc::now();
            let duration = (completed_at - timing.started_at).num_milliseconds() as f64 / 1000.0;
            let throughput = if duration > 0.0 {
                items_processed as f64 / duration
            } else {
                0.0
            };

            timing.completed_at = Some(completed_at);
            timing.duration_seconds = Some(duration);
            timing.items_processed = items_processed;
            timing.throughput_per_second = Some(throughput);

            info!(
                "✅ {} completed in {:.2}s ({} items, {:.1} items/sec)",
                stage.as_str(),
                duration,
                items_processed,
                throughput
            );
        }
    }

    /// Update valid contig count after filtering
    pub fn set_valid_contigs(&mut self, count: usize) {
        self.run_metadata.valid_contigs = count;
        let filtered = self.run_metadata.total_contigs - count;
        info!(
            "📊 Filtering complete: {} valid contigs, {} filtered out ({:.1}% retained)",
            count,
            filtered,
            (count as f64 / self.run_metadata.total_contigs as f64) * 100.0
        );
    }

    /// Log feature extraction progress
    pub fn log_feature_extraction_progress(&self, processed: usize, total: usize) {
        if processed % 100 == 0 || processed == total {
            let percent = (processed as f64 / total as f64) * 100.0;
            debug!(
                "   Extracting features: {}/{} ({:.1}%)",
                processed, total, percent
            );
        }
    }

    /// Log normalization details
    pub fn log_normalization(&self, method: &str, n_features: usize, n_samples: usize) {
        info!(
            "   Normalization: {} method, {} features, {} samples",
            method, n_features, n_samples
        );
    }

    /// Log clustering progress
    pub fn log_clustering_iteration(&self, iteration: usize, changed_count: usize) {
        debug!(
            "   Clustering iteration {}: {} assignments changed",
            iteration, changed_count
        );
    }

    /// Log clustering completion
    pub fn log_clustering_complete(&self, iterations: usize, num_bins: usize) {
        info!(
            "   Clustering converged after {} iterations, created {} bins",
            iterations, num_bins
        );
    }

    /// Add bin quality metrics
    pub fn add_bin_metrics(&mut self, metrics: BinQualityMetrics) {
        debug!(
            "   Bin {}: {} contigs, {:.1} MB, N50={}, completeness={:.1}%",
            metrics.bin_id,
            metrics.contig_count,
            metrics.total_bp as f64 / 1_000_000.0,
            metrics.n50,
            metrics.completeness.unwrap_or(0.0)
        );
        self.bin_metrics.push(metrics);
    }

    /// Complete the classification run
    pub fn complete_run(&mut self) -> Result<ClassificationResults> {
        self.run_metadata.completed_at = Some(Utc::now());

        let total_time = (self.run_metadata.completed_at.unwrap() - self.run_metadata.started_at)
            .num_milliseconds() as f64
            / 1000.0;

        // Calculate overall statistics
        let num_bins = self.bin_metrics.len();
        let total_classified = self
            .bin_metrics
            .iter()
            .map(|b| b.contig_count)
            .sum::<usize>();
        let total_bp_classified = self.bin_metrics.iter().map(|b| b.total_bp).sum::<usize>();
        let avg_bin_size = if num_bins > 0 {
            total_classified as f64 / num_bins as f64
        } else {
            0.0
        };
        let avg_bin_size_bp = if num_bins > 0 {
            total_bp_classified as f64 / num_bins as f64
        } else {
            0.0
        };

        let overall_stats = OverallStatistics {
            total_classified,
            total_filtered: self.run_metadata.total_contigs - total_classified,
            num_bins,
            total_bp_classified,
            avg_bin_size,
            avg_bin_size_bp,
            total_time_seconds: total_time,
            throughput_contigs_per_sec: total_classified as f64 / total_time,
        };

        // Sort bin metrics by bin ID
        self.bin_metrics.sort_by_key(|m| m.bin_id);

        let results = ClassificationResults {
            metadata: self.run_metadata.clone(),
            bin_metrics: self.bin_metrics.clone(),
            stage_timings: self.stage_timings.values().cloned().collect(),
            overall_stats: overall_stats.clone(),
        };

        // Print summary
        info!("🎯 Classification Complete");
        info!("   Total time: {:.2}s", total_time);
        info!(
            "   Bins created: {} ({:.1} contigs/bin avg)",
            num_bins, avg_bin_size
        );
        info!(
            "   Total classified: {} contigs ({} bp)",
            total_classified, total_bp_classified
        );
        info!(
            "   Throughput: {:.1} contigs/sec",
            overall_stats.throughput_contigs_per_sec
        );

        // Save detailed results
        self.save_results(&results)?;

        Ok(results)
    }

    /// Save comprehensive results to output directory
    fn save_results(&self, results: &ClassificationResults) -> Result<()> {
        // Save full JSON results
        self.output_manager
            .save_intermediate(
                PipelineSection::Classification,
                "classification_results",
                results,
                serde_json::json!({
                    "method": results.metadata.method,
                    "num_bins": results.bin_metrics.len(),
                }),
            )
            .context("Failed to save classification results JSON")?;

        // Save bin metrics as TSV
        self.save_bin_metrics_tsv(&results.bin_metrics)?;

        // Save stage timings as TSV
        self.save_stage_timings_tsv(&results.stage_timings)?;

        // Save summary report
        self.save_summary_report(results)?;

        Ok(())
    }

    /// Save bin metrics as TSV
    fn save_bin_metrics_tsv(&self, metrics: &[BinQualityMetrics]) -> Result<()> {
        let headers = vec![
            "bin_id".to_string(),
            "contig_count".to_string(),
            "total_bp".to_string(),
            "avg_contig_length".to_string(),
            "median_contig_length".to_string(),
            "avg_coverage".to_string(),
            "coverage_std".to_string(),
            "avg_gc_content".to_string(),
            "gc_std".to_string(),
            "completeness".to_string(),
            "contamination".to_string(),
            "avg_confidence".to_string(),
            "n50".to_string(),
            "longest_contig".to_string(),
        ];

        let rows: Vec<Vec<String>> = metrics
            .iter()
            .map(|m| {
                vec![
                    m.bin_id.to_string(),
                    m.contig_count.to_string(),
                    m.total_bp.to_string(),
                    format!("{:.2}", m.avg_contig_length),
                    m.median_contig_length.to_string(),
                    format!("{:.2}", m.avg_coverage),
                    format!("{:.2}", m.coverage_std),
                    format!("{:.4}", m.avg_gc_content),
                    format!("{:.4}", m.gc_std),
                    m.completeness
                        .map(|v| format!("{:.2}", v))
                        .unwrap_or_else(|| "NA".to_string()),
                    m.contamination
                        .map(|v| format!("{:.2}", v))
                        .unwrap_or_else(|| "NA".to_string()),
                    format!("{:.4}", m.avg_confidence),
                    m.n50.to_string(),
                    m.longest_contig.to_string(),
                ]
            })
            .collect();

        self.output_manager
            .save_tsv(
                PipelineSection::Classification,
                "bin_quality_metrics",
                &headers,
                &rows,
                serde_json::json!({ "num_bins": metrics.len() }),
            )
            .context("Failed to save bin metrics TSV")?;

        info!(
            "💾 Saved bin quality metrics: bin_quality_metrics.tsv ({} bins)",
            metrics.len()
        );

        Ok(())
    }

    /// Save stage timings as TSV
    fn save_stage_timings_tsv(&self, timings: &[StageTiming]) -> Result<()> {
        let headers = vec![
            "stage".to_string(),
            "duration_seconds".to_string(),
            "items_processed".to_string(),
            "throughput_per_second".to_string(),
        ];

        let rows: Vec<Vec<String>> = timings
            .iter()
            .map(|t| {
                vec![
                    t.stage.as_str().to_string(),
                    t.duration_seconds
                        .map(|d| format!("{:.3}", d))
                        .unwrap_or_else(|| "NA".to_string()),
                    t.items_processed.to_string(),
                    t.throughput_per_second
                        .map(|tp| format!("{:.2}", tp))
                        .unwrap_or_else(|| "NA".to_string()),
                ]
            })
            .collect();

        self.output_manager
            .save_tsv(
                PipelineSection::Classification,
                "stage_timings",
                &headers,
                &rows,
                serde_json::json!({ "num_stages": timings.len() }),
            )
            .context("Failed to save stage timings TSV")?;

        Ok(())
    }

    /// Save human-readable summary report
    fn save_summary_report(&self, results: &ClassificationResults) -> Result<()> {
        let mut report = String::new();

        report.push_str("# Classification Run Summary\n\n");
        report.push_str(&format!(
            "**Method:** {:?}\n",
            results.metadata.method
        ));
        report.push_str(&format!(
            "**Started:** {}\n",
            results.metadata.started_at.format("%Y-%m-%d %H:%M:%S")
        ));
        report.push_str(&format!(
            "**Completed:** {}\n",
            results
                .metadata
                .completed_at
                .map(|t| t.format("%Y-%m-%d %H:%M:%S").to_string())
                .unwrap_or_else(|| "In Progress".to_string())
        ));
        report.push_str(&format!(
            "**Duration:** {:.2} seconds\n\n",
            results.overall_stats.total_time_seconds
        ));

        report.push_str("## Overall Statistics\n\n");
        report.push_str(&format!(
            "- **Total Contigs:** {}\n",
            results.metadata.total_contigs
        ));
        report.push_str(&format!(
            "- **Valid Contigs:** {}\n",
            results.metadata.valid_contigs
        ));
        report.push_str(&format!(
            "- **Bins Created:** {}\n",
            results.overall_stats.num_bins
        ));
        report.push_str(&format!(
            "- **Total Classified:** {} contigs ({:.2} MB)\n",
            results.overall_stats.total_classified,
            results.overall_stats.total_bp_classified as f64 / 1_000_000.0
        ));
        report.push_str(&format!(
            "- **Average Bin Size:** {:.1} contigs ({:.2} MB)\n",
            results.overall_stats.avg_bin_size,
            results.overall_stats.avg_bin_size_bp / 1_000_000.0
        ));
        report.push_str(&format!(
            "- **Throughput:** {:.1} contigs/sec\n\n",
            results.overall_stats.throughput_contigs_per_sec
        ));

        report.push_str("## Bin Summary\n\n");
        report.push_str("| Bin | Contigs | Size (MB) | N50 | Avg Coverage | GC% | Completeness | Contamination |\n");
        report.push_str("|-----|---------|-----------|-----|--------------|-----|--------------|---------------|\n");

        for bin in &results.bin_metrics {
            report.push_str(&format!(
                "| {} | {} | {:.2} | {} | {:.1}x | {:.1}% | {:.1}% | {:.1}% |\n",
                bin.bin_id,
                bin.contig_count,
                bin.total_bp as f64 / 1_000_000.0,
                bin.n50,
                bin.avg_coverage,
                bin.avg_gc_content * 100.0,
                bin.completeness.unwrap_or(0.0),
                bin.contamination.unwrap_or(0.0)
            ));
        }

        report.push_str("\n## Stage Timings\n\n");
        for timing in &results.stage_timings {
            if let Some(duration) = timing.duration_seconds {
                report.push_str(&format!(
                    "- **{}:** {:.2}s ({} items, {:.1} items/sec)\n",
                    timing.stage.as_str(),
                    duration,
                    timing.items_processed,
                    timing.throughput_per_second.unwrap_or(0.0)
                ));
            }
        }

        // Save as markdown in docs directory
        let summary_path = self
            .output_manager
            .get_section_dir(&PipelineSection::Classification)
            .join("SUMMARY.md");

        std::fs::write(&summary_path, report)
            .with_context(|| format!("Failed to write summary report: {}", summary_path.display()))?;

        info!("📄 Saved summary report: SUMMARY.md");

        Ok(())
    }
}

/// Calculate bin quality metrics from contig assignments
pub fn calculate_bin_metrics(
    bin_id: usize,
    contigs: &[(String, usize, f64, f64)], // (sequence, length, coverage, gc_content)
    confidences: &[f64],
) -> BinQualityMetrics {
    let contig_count = contigs.len();
    let total_bp: usize = contigs.iter().map(|(_, len, _, _)| len).sum();

    let avg_contig_length = if contig_count > 0 {
        total_bp as f64 / contig_count as f64
    } else {
        0.0
    };

    // Calculate median length
    let mut lengths: Vec<usize> = contigs.iter().map(|(_, len, _, _)| *len).collect();
    lengths.sort_unstable();
    let median_contig_length = if !lengths.is_empty() {
        lengths[lengths.len() / 2]
    } else {
        0
    };

    // Calculate coverage statistics
    let coverages: Vec<f64> = contigs.iter().map(|(_, _, cov, _)| *cov).collect();
    let avg_coverage = if !coverages.is_empty() {
        coverages.iter().sum::<f64>() / coverages.len() as f64
    } else {
        0.0
    };

    let coverage_std = if coverages.len() > 1 {
        let variance = coverages
            .iter()
            .map(|c| (c - avg_coverage).powi(2))
            .sum::<f64>()
            / coverages.len() as f64;
        variance.sqrt()
    } else {
        0.0
    };

    // Calculate GC content statistics
    let gc_contents: Vec<f64> = contigs.iter().map(|(_, _, _, gc)| *gc).collect();
    let avg_gc_content = if !gc_contents.is_empty() {
        gc_contents.iter().sum::<f64>() / gc_contents.len() as f64
    } else {
        0.0
    };

    let gc_std = if gc_contents.len() > 1 {
        let variance = gc_contents
            .iter()
            .map(|gc| (gc - avg_gc_content).powi(2))
            .sum::<f64>()
            / gc_contents.len() as f64;
        variance.sqrt()
    } else {
        0.0
    };

    // Calculate N50
    let n50 = calculate_n50(&lengths);

    // Find longest contig
    let longest_contig = lengths.iter().max().copied().unwrap_or(0);

    // Average confidence
    let avg_confidence = if !confidences.is_empty() {
        confidences.iter().sum::<f64>() / confidences.len() as f64
    } else {
        0.0
    };

    BinQualityMetrics {
        bin_id,
        contig_count,
        total_bp,
        avg_contig_length,
        median_contig_length,
        avg_coverage,
        coverage_std,
        avg_gc_content,
        gc_std,
        completeness: None, // Requires external tools like CheckM
        contamination: None, // Requires external tools like CheckM
        avg_confidence,
        n50,
        longest_contig,
    }
}

/// Calculate N50 statistic
fn calculate_n50(sorted_lengths: &[usize]) -> usize {
    let total_length: usize = sorted_lengths.iter().sum();
    let half_length = total_length / 2;

    let mut cumulative = 0;
    for &length in sorted_lengths.iter().rev() {
        cumulative += length;
        if cumulative >= half_length {
            return length;
        }
    }

    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_calculate_n50() {
        let lengths = vec![100, 200, 300, 400, 500];
        let n50 = calculate_n50(&lengths);
        // Total = 1500, half = 750
        // 500 + 400 = 900 >= 750, so N50 = 400
        assert_eq!(n50, 400);
    }

    #[test]
    fn test_bin_metrics_calculation() {
        let contigs = vec![
            ("seq1".to_string(), 1000, 10.0, 0.5),
            ("seq2".to_string(), 2000, 15.0, 0.55),
            ("seq3".to_string(), 1500, 12.0, 0.52),
        ];
        let confidences = vec![0.9, 0.85, 0.88];

        let metrics = calculate_bin_metrics(1, &contigs, &confidences);

        assert_eq!(metrics.bin_id, 1);
        assert_eq!(metrics.contig_count, 3);
        assert_eq!(metrics.total_bp, 4500);
        assert_eq!(metrics.avg_contig_length, 1500.0);
        assert_eq!(metrics.median_contig_length, 1500);
        assert!((metrics.avg_coverage - 12.333).abs() < 0.01);
        assert!((metrics.avg_gc_content - 0.523).abs() < 0.01);
        assert_eq!(metrics.longest_contig, 2000);
    }
}
