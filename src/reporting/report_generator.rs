//! Report generation module for creating comprehensive analysis reports

use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use tracing::info;

use crate::core::data_structures::Contig;
use crate::pipeline::complete_integration::{AssemblyResults, TaxonomicClassification, AbundanceProfile};
use crate::utils::kraken_reporter::{KrakenClassification, KrakenReporter};

/// Complete analysis report structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisReport {
    pub sample_name: String,
    pub timestamp: chrono::DateTime<chrono::Utc>,
    pub summary: ReportSummary,
    pub quality_metrics: QualityMetrics,
    pub taxonomic_composition: Vec<TaxonomicClassification>,
    pub abundance_data: AbundanceProfile,
    pub performance_metrics: PerformanceMetrics,
}

/// Summary statistics for the report
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportSummary {
    pub total_contigs: usize,
    pub total_length: usize,
    pub n50: usize,
    pub mean_coverage: f64,
    pub unique_species: usize,
    pub diversity_index: f64,
}

/// Quality metrics for assembly and classification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
    pub assembly_completeness: f64,
    pub classification_confidence: f64,
    pub coverage_uniformity: f64,
}

/// Performance metrics for the analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
    pub total_processing_time: std::time::Duration,
    pub peak_memory_usage: usize,
    pub reads_processed: u64,
    pub errors_corrected: u64,
    pub repeats_resolved: u64,
}

/// Report generator for creating various report formats
pub struct ReportGenerator {
    output_formats: OutputFormats,
}

#[derive(Debug, Clone)]
pub struct OutputFormats {
    pub json: bool,
    pub html: bool,
    pub tsv: bool,
    pub markdown: bool,
    pub kraken: bool,
}

impl Default for OutputFormats {
    fn default() -> Self {
        Self {
            json: true,
            html: true,
            tsv: true,
            markdown: true,
            kraken: true,
        }
    }
}

impl ReportGenerator {
    /// Create a new report generator with default formats
    pub fn new() -> Self {
        Self {
            output_formats: OutputFormats::default(),
        }
    }

    /// Create with custom output formats
    pub fn with_formats(output_formats: OutputFormats) -> Self {
        Self { output_formats }
    }

    /// Generate a complete analysis report
    pub fn generate_report(
        &self,
        sample_name: &str,
        assembly_results: &AssemblyResults,
        classifications: &[TaxonomicClassification],
        abundance_profile: &AbundanceProfile,
        elapsed_time: std::time::Duration,
    ) -> AnalysisReport {
        let unique_species = classifications
            .iter()
            .map(|c| &c.taxonomy_name)
            .collect::<std::collections::HashSet<_>>()
            .len();

        AnalysisReport {
            sample_name: sample_name.to_string(),
            timestamp: chrono::Utc::now(),
            summary: ReportSummary {
                total_contigs: assembly_results.contigs.len(),
                total_length: assembly_results.assembly_stats.total_length,
                n50: assembly_results.assembly_stats.n50,
                mean_coverage: assembly_results.assembly_stats.coverage_mean,
                unique_species,
                diversity_index: Self::calculate_shannon_diversity(&abundance_profile.abundant_kmers),
            },
            quality_metrics: QualityMetrics {
                assembly_completeness: Self::calculate_assembly_completeness(assembly_results),
                classification_confidence: if !classifications.is_empty() {
                    classifications.iter().map(|c| c.confidence).sum::<f64>()
                        / classifications.len() as f64
                } else {
                    0.0
                },
                coverage_uniformity: Self::calculate_coverage_uniformity(&assembly_results.contigs),
            },
            taxonomic_composition: classifications.to_vec(),
            abundance_data: abundance_profile.clone(),
            performance_metrics: PerformanceMetrics {
                total_processing_time: elapsed_time,
                peak_memory_usage: 0, // TODO: Implement proper memory tracking
                reads_processed: abundance_profile.total_kmers / 4, // Estimate using k=4
                errors_corrected: 0,  // TODO: Track from QC pipeline
                repeats_resolved: 0,  // TODO: Track from assembly
            },
        }
    }

    /// Write all report files
    pub async fn write_report_files(
        &self,
        report: &AnalysisReport,
        output_dir: &Path,
    ) -> Result<()> {
        // JSON report
        if self.output_formats.json {
            let json_path = output_dir.join(format!("{}_report.json", report.sample_name));
            let json_content = serde_json::to_string_pretty(report)?;
            tokio::fs::write(&json_path, json_content).await?;
            info!("📄 JSON report written to: {}", json_path.display());
        }

        // HTML report
        if self.output_formats.html {
            let html_path = output_dir.join(format!("{}_report.html", report.sample_name));
            let html_content = self.generate_html_report(report)?;
            tokio::fs::write(&html_path, html_content).await?;
            info!("🌐 HTML report written to: {}", html_path.display());
        }

        // TSV summary
        if self.output_formats.tsv {
            let tsv_path = output_dir.join(format!("{}_summary.tsv", report.sample_name));
            let tsv_content = self.generate_tsv_summary(report)?;
            tokio::fs::write(&tsv_path, tsv_content).await?;
            info!("📊 TSV summary written to: {}", tsv_path.display());
        }

        // Markdown report
        if self.output_formats.markdown {
            let md_path = output_dir.join(format!("{}_report.md", report.sample_name));
            let md_content = self.generate_markdown_report(report)?;
            tokio::fs::write(&md_path, md_content).await?;
            info!("📝 Markdown report written to: {}", md_path.display());
        }

        // Kraken reports
        if self.output_formats.kraken {
            self.generate_kraken_reports(report, output_dir).await?;
        }

        Ok(())
    }

    /// Generate Kraken-style reports
    async fn generate_kraken_reports(&self, report: &AnalysisReport, output_dir: &Path) -> Result<()> {
        let mut kraken_reporter = KrakenReporter::new(report.sample_name.clone());

        for classification in &report.taxonomic_composition {
            let kraken_classification: KrakenClassification = classification.clone().into();
            kraken_reporter.add_classification(kraken_classification);
        }

        // Standard Kraken output
        let kraken_output_path = output_dir.join(format!("{}.kraken", report.sample_name));
        kraken_reporter.write_kraken_output(&kraken_output_path)?;
        info!("🧬 Kraken output written to: {}", kraken_output_path.display());

        // Kraken report
        let kraken_report_path = output_dir.join(format!("{}_kraken_report.txt", report.sample_name));
        kraken_reporter.write_kraken_report(&kraken_report_path)?;
        info!("📊 Kraken report written to: {}", kraken_report_path.display());

        // Enhanced JSON
        let enhanced_json_path = output_dir.join(format!("{}_enhanced_kraken.json", report.sample_name));
        kraken_reporter.write_enhanced_json_report(&enhanced_json_path)?;
        info!("📄 Enhanced Kraken JSON report written to: {}", enhanced_json_path.display());

        Ok(())
    }

    /// Generate HTML report
    fn generate_html_report(&self, report: &AnalysisReport) -> Result<String> {
        let html = format!(
            r#"<!DOCTYPE html>
<html>
<head>
    <title>Metagenomics Analysis Report - {}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 30px; border-radius: 8px; margin-bottom: 30px; }}
        .header h1 {{ margin: 0 0 10px 0; }}
        .header .timestamp {{ opacity: 0.9; font-size: 14px; }}
        .section {{ margin: 25px 0; padding: 20px; border-left: 4px solid #667eea; background: #f9f9f9; border-radius: 4px; }}
        .section h2 {{ margin-top: 0; color: #333; }}
        .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .metric {{ background: white; padding: 15px; border-radius: 6px; border: 1px solid #e0e0e0; }}
        .metric .label {{ color: #666; font-size: 12px; text-transform: uppercase; margin-bottom: 5px; }}
        .metric .value {{ font-size: 24px; font-weight: bold; color: #667eea; }}
        .taxa-table {{ width: 100%; border-collapse: collapse; margin-top: 15px; }}
        .taxa-table th, .taxa-table td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        .taxa-table th {{ background: #667eea; color: white; }}
        .taxa-table tr:hover {{ background: #f5f5f5; }}
        .footer {{ margin-top: 40px; padding-top: 20px; border-top: 2px solid #e0e0e0; text-align: center; color: #666; font-size: 14px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Metagenomics Analysis Report</h1>
            <div class="timestamp">Sample: {} | Generated: {}</div>
        </div>

        <div class="section">
            <h2>📊 Assembly Summary</h2>
            <div class="metrics">
                <div class="metric">
                    <div class="label">Total Contigs</div>
                    <div class="value">{}</div>
                </div>
                <div class="metric">
                    <div class="label">Total Length</div>
                    <div class="value">{:.2} Mb</div>
                </div>
                <div class="metric">
                    <div class="label">N50</div>
                    <div class="value">{} bp</div>
                </div>
                <div class="metric">
                    <div class="label">Mean Coverage</div>
                    <div class="value">{:.1}x</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>🦠 Taxonomic Composition</h2>
            <div class="metrics">
                <div class="metric">
                    <div class="label">Unique Species</div>
                    <div class="value">{}</div>
                </div>
                <div class="metric">
                    <div class="label">Shannon Diversity</div>
                    <div class="value">{:.3}</div>
                </div>
                <div class="metric">
                    <div class="label">Classification Confidence</div>
                    <div class="value">{:.1}%</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>⚡ Performance Metrics</h2>
            <div class="metrics">
                <div class="metric">
                    <div class="label">Processing Time</div>
                    <div class="value">{:.2}s</div>
                </div>
                <div class="metric">
                    <div class="label">Reads Processed</div>
                    <div class="value">{}</div>
                </div>
            </div>
        </div>

        <div class="footer">
            Generated by MetaForge Metagenomics Pipeline
        </div>
    </div>
</body>
</html>"#,
            report.sample_name,
            report.sample_name,
            report.timestamp.format("%Y-%m-%d %H:%M:%S UTC"),
            report.summary.total_contigs,
            report.summary.total_length as f64 / 1_000_000.0,
            report.summary.n50,
            report.summary.mean_coverage,
            report.summary.unique_species,
            report.summary.diversity_index,
            report.quality_metrics.classification_confidence * 100.0,
            report.performance_metrics.total_processing_time.as_secs_f64(),
            report.performance_metrics.reads_processed,
        );

        Ok(html)
    }

    /// Generate TSV summary
    fn generate_tsv_summary(&self, report: &AnalysisReport) -> Result<String> {
        let mut tsv = String::from("Metric\tValue\n");
        tsv.push_str(&format!("Sample Name\t{}\n", report.sample_name));
        tsv.push_str(&format!("Timestamp\t{}\n", report.timestamp));
        tsv.push_str(&format!("Total Contigs\t{}\n", report.summary.total_contigs));
        tsv.push_str(&format!("Total Length (bp)\t{}\n", report.summary.total_length));
        tsv.push_str(&format!("N50 (bp)\t{}\n", report.summary.n50));
        tsv.push_str(&format!("Mean Coverage\t{:.2}\n", report.summary.mean_coverage));
        tsv.push_str(&format!("Unique Species\t{}\n", report.summary.unique_species));
        tsv.push_str(&format!("Shannon Diversity\t{:.4}\n", report.summary.diversity_index));
        tsv.push_str(&format!(
            "Classification Confidence\t{:.4}\n",
            report.quality_metrics.classification_confidence
        ));
        tsv.push_str(&format!(
            "Processing Time (s)\t{:.2}\n",
            report.performance_metrics.total_processing_time.as_secs_f64()
        ));
        tsv.push_str(&format!("Reads Processed\t{}\n", report.performance_metrics.reads_processed));

        Ok(tsv)
    }

    /// Generate Markdown report
    fn generate_markdown_report(&self, report: &AnalysisReport) -> Result<String> {
        let mut md = format!("# 🧬 Metagenomics Analysis Report\n\n");
        md.push_str(&format!("**Sample:** {}\n", report.sample_name));
        md.push_str(&format!("**Generated:** {}\n\n", report.timestamp));

        md.push_str("## 📊 Assembly Summary\n\n");
        md.push_str(&format!("- **Total Contigs:** {}\n", report.summary.total_contigs));
        md.push_str(&format!(
            "- **Total Length:** {:.2} Mb\n",
            report.summary.total_length as f64 / 1_000_000.0
        ));
        md.push_str(&format!("- **N50:** {} bp\n", report.summary.n50));
        md.push_str(&format!("- **Mean Coverage:** {:.1}x\n\n", report.summary.mean_coverage));

        md.push_str("## 🦠 Taxonomic Composition\n\n");
        md.push_str(&format!("- **Unique Species:** {}\n", report.summary.unique_species));
        md.push_str(&format!("- **Shannon Diversity:** {:.3}\n", report.summary.diversity_index));
        md.push_str(&format!(
            "- **Classification Confidence:** {:.1}%\n\n",
            report.quality_metrics.classification_confidence * 100.0
        ));

        md.push_str("## ⚡ Performance Metrics\n\n");
        md.push_str(&format!(
            "- **Processing Time:** {:.2}s\n",
            report.performance_metrics.total_processing_time.as_secs_f64()
        ));
        md.push_str(&format!("- **Reads Processed:** {}\n\n", report.performance_metrics.reads_processed));

        md.push_str("---\n\n");
        md.push_str("*Generated by MetaForge Metagenomics Pipeline*\n");

        Ok(md)
    }

    /// Calculate Shannon diversity index
    fn calculate_shannon_diversity(abundant_kmers: &HashMap<u64, f64>) -> f64 {
        if abundant_kmers.is_empty() {
            return 0.0;
        }

        let total: f64 = abundant_kmers.values().sum();
        if total == 0.0 {
            return 0.0;
        }

        -abundant_kmers
            .values()
            .map(|&count| {
                let p = count / total;
                if p > 0.0 {
                    p * p.ln()
                } else {
                    0.0
                }
            })
            .sum::<f64>()
    }

    /// Calculate assembly completeness score
    fn calculate_assembly_completeness(assembly_results: &AssemblyResults) -> f64 {
        if assembly_results.contigs.is_empty() {
            return 0.0;
        }

        // Simple heuristic: ratio of N50 to largest contig
        let n50 = assembly_results.assembly_stats.n50 as f64;
        let largest = assembly_results.assembly_stats.largest_contig as f64;

        if largest == 0.0 {
            0.0
        } else {
            (n50 / largest).min(1.0)
        }
    }

    /// Calculate coverage uniformity
    fn calculate_coverage_uniformity(contigs: &[Contig]) -> f64 {
        if contigs.len() <= 1 {
            return 1.0;
        }

        let mean: f64 = contigs.iter().map(|c| c.coverage).sum::<f64>() / contigs.len() as f64;

        if mean == 0.0 {
            return 0.0;
        }

        let std_dev = Self::calculate_std_dev(contigs, mean);

        // Coefficient of variation (CV) - lower is more uniform
        let cv = std_dev / mean;

        // Convert to uniformity score (0-1, higher is better)
        1.0 / (1.0 + cv)
    }

    fn calculate_std_dev(contigs: &[Contig], mean: f64) -> f64 {
        if contigs.len() <= 1 {
            return 0.0;
        }

        let variance: f64 = contigs
            .iter()
            .map(|c| {
                let diff = c.coverage - mean;
                diff * diff
            })
            .sum::<f64>()
            / (contigs.len() - 1) as f64;

        variance.sqrt()
    }
}

impl Default for ReportGenerator {
    fn default() -> Self {
        Self::new()
    }
}
