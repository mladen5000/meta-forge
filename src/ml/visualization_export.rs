//! Visualization data export for classification results
//!
//! Exports classification data in formats suitable for visualization tools
//! like Python matplotlib, R ggplot2, or web-based tools.

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

use crate::ml::classification_reporter::{BinQualityMetrics, ClassificationResults};
use crate::ml::simple_classifier::ContigClassification;
use crate::utils::intermediate_output::{IntermediateOutputManager, PipelineSection};

/// Visualization export manager
pub struct VisualizationExporter {
    output_manager: IntermediateOutputManager,
}

/// Data for bin size distribution plot
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinSizeDistribution {
    pub bin_ids: Vec<usize>,
    pub contig_counts: Vec<usize>,
    pub total_bp: Vec<usize>,
}

/// Data for coverage vs GC content scatter plot
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageGcScatter {
    pub contigs: Vec<ContigPoint>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigPoint {
    pub contig_id: usize,
    pub bin_id: usize,
    pub coverage: f64,
    pub gc_content: f64,
    pub length: usize,
    pub confidence: f64,
}

/// Data for bin quality heatmap
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinQualityHeatmap {
    pub bins: Vec<usize>,
    pub metrics: Vec<QualityMetric>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetric {
    pub metric_name: String,
    pub values: Vec<f64>,
}

/// Contig length distribution per bin
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigLengthDistribution {
    pub bin_id: usize,
    pub lengths: Vec<usize>,
}

impl VisualizationExporter {
    /// Create a new visualization exporter
    pub fn new(output_manager: IntermediateOutputManager) -> Self {
        Self { output_manager }
    }

    /// Export all visualization data
    pub fn export_all(
        &self,
        results: &ClassificationResults,
        classifications: &[ContigClassification],
        contigs_data: &[(usize, f64, f64, usize)], // (contig_id, coverage, gc, length)
    ) -> Result<()> {
        tracing::info!("📊 Exporting visualization data...");

        self.export_bin_size_distribution(results)?;
        self.export_coverage_gc_scatter(classifications, contigs_data)?;
        self.export_bin_quality_heatmap(results)?;
        self.export_contig_length_distributions(classifications, contigs_data)?;
        self.export_python_plotting_script()?;
        self.export_r_plotting_script()?;

        tracing::info!("✅ Visualization data exported successfully");
        Ok(())
    }

    /// Export bin size distribution
    fn export_bin_size_distribution(&self, results: &ClassificationResults) -> Result<()> {
        let bin_ids: Vec<usize> = results.bin_metrics.iter().map(|m| m.bin_id).collect();
        let contig_counts: Vec<usize> = results
            .bin_metrics
            .iter()
            .map(|m| m.contig_count)
            .collect();
        let total_bp: Vec<usize> = results.bin_metrics.iter().map(|m| m.total_bp).collect();

        let distribution = BinSizeDistribution {
            bin_ids,
            contig_counts,
            total_bp,
        };

        self.output_manager.save_intermediate(
            PipelineSection::Classification,
            "viz_bin_size_distribution",
            &distribution,
            serde_json::json!({ "num_bins": results.bin_metrics.len() }),
        )?;

        // Also save as CSV for easy import
        let headers = vec!["bin_id".to_string(), "contig_count".to_string(), "total_bp".to_string()];
        let rows: Vec<Vec<String>> = results
            .bin_metrics
            .iter()
            .map(|m| {
                vec![
                    m.bin_id.to_string(),
                    m.contig_count.to_string(),
                    m.total_bp.to_string(),
                ]
            })
            .collect();

        self.output_manager.save_tsv(
            PipelineSection::Classification,
            "viz_bin_sizes",
            &headers,
            &rows,
            serde_json::json!({}),
        )?;

        Ok(())
    }

    /// Export coverage vs GC content scatter plot data
    fn export_coverage_gc_scatter(
        &self,
        classifications: &[ContigClassification],
        contigs_data: &[(usize, f64, f64, usize)],
    ) -> Result<()> {
        let mut contig_map: std::collections::HashMap<usize, &ContigClassification> =
            classifications
                .iter()
                .map(|c| (c.contig_id, c))
                .collect();

        let contigs: Vec<ContigPoint> = contigs_data
            .iter()
            .filter_map(|(contig_id, coverage, gc, length)| {
                contig_map.get(contig_id).map(|classification| ContigPoint {
                    contig_id: *contig_id,
                    bin_id: classification.bin_id,
                    coverage: *coverage,
                    gc_content: *gc,
                    length: *length,
                    confidence: classification.confidence,
                })
            })
            .collect();

        let scatter = CoverageGcScatter { contigs };

        self.output_manager.save_intermediate(
            PipelineSection::Classification,
            "viz_coverage_gc_scatter",
            &scatter,
            serde_json::json!({ "num_contigs": scatter.contigs.len() }),
        )?;

        // Also save as TSV
        let headers = vec![
            "contig_id".to_string(),
            "bin_id".to_string(),
            "coverage".to_string(),
            "gc_content".to_string(),
            "length".to_string(),
            "confidence".to_string(),
        ];
        let rows: Vec<Vec<String>> = scatter
            .contigs
            .iter()
            .map(|c| {
                vec![
                    c.contig_id.to_string(),
                    c.bin_id.to_string(),
                    format!("{:.2}", c.coverage),
                    format!("{:.4}", c.gc_content),
                    c.length.to_string(),
                    format!("{:.4}", c.confidence),
                ]
            })
            .collect();

        self.output_manager.save_tsv(
            PipelineSection::Classification,
            "viz_coverage_gc",
            &headers,
            &rows,
            serde_json::json!({}),
        )?;

        Ok(())
    }

    /// Export bin quality heatmap data
    fn export_bin_quality_heatmap(&self, results: &ClassificationResults) -> Result<()> {
        let bins: Vec<usize> = results.bin_metrics.iter().map(|m| m.bin_id).collect();

        let metrics = vec![
            QualityMetric {
                metric_name: "Average Coverage".to_string(),
                values: results
                    .bin_metrics
                    .iter()
                    .map(|m| m.avg_coverage)
                    .collect(),
            },
            QualityMetric {
                metric_name: "GC Content".to_string(),
                values: results
                    .bin_metrics
                    .iter()
                    .map(|m| m.avg_gc_content * 100.0)
                    .collect(),
            },
            QualityMetric {
                metric_name: "N50 (kb)".to_string(),
                values: results
                    .bin_metrics
                    .iter()
                    .map(|m| m.n50 as f64 / 1000.0)
                    .collect(),
            },
            QualityMetric {
                metric_name: "Confidence".to_string(),
                values: results
                    .bin_metrics
                    .iter()
                    .map(|m| m.avg_confidence)
                    .collect(),
            },
        ];

        let heatmap = BinQualityHeatmap { bins, metrics };

        self.output_manager.save_intermediate(
            PipelineSection::Classification,
            "viz_quality_heatmap",
            &heatmap,
            serde_json::json!({}),
        )?;

        Ok(())
    }

    /// Export contig length distributions per bin
    fn export_contig_length_distributions(
        &self,
        classifications: &[ContigClassification],
        contigs_data: &[(usize, f64, f64, usize)],
    ) -> Result<()> {
        use std::collections::HashMap;

        // Group lengths by bin
        let mut bin_lengths: HashMap<usize, Vec<usize>> = HashMap::new();

        for classification in classifications {
            if let Some((_id, _cov, _gc, length)) = contigs_data
                .iter()
                .find(|(id, _, _, _)| *id == classification.contig_id)
            {
                bin_lengths
                    .entry(classification.bin_id)
                    .or_insert_with(Vec::new)
                    .push(*length);
            }
        }

        let distributions: Vec<ContigLengthDistribution> = bin_lengths
            .into_iter()
            .map(|(bin_id, mut lengths)| {
                lengths.sort_unstable();
                ContigLengthDistribution { bin_id, lengths }
            })
            .collect();

        self.output_manager.save_intermediate(
            PipelineSection::Classification,
            "viz_length_distributions",
            &distributions,
            serde_json::json!({}),
        )?;

        Ok(())
    }

    /// Export Python plotting script
    fn export_python_plotting_script(&self) -> Result<()> {
        let script = r#"#!/usr/bin/env python3
"""
Visualization script for MetaForge classification results
Requires: matplotlib, seaborn, pandas, numpy
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
sns.set_theme(style="whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def load_json(filename):
    """Load JSON data file"""
    with open(filename) as f:
        return json.load(f)

def plot_bin_size_distribution(data_dir):
    """Plot bin size distribution"""
    df = pd.read_csv(data_dir / "viz_bin_sizes.tsv", sep='\t')

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Contig counts
    axes[0].bar(df['bin_id'], df['contig_count'], color='steelblue')
    axes[0].set_xlabel('Bin ID')
    axes[0].set_ylabel('Number of Contigs')
    axes[0].set_title('Contigs per Bin')

    # Total base pairs
    axes[1].bar(df['bin_id'], df['total_bp'] / 1e6, color='coral')
    axes[1].set_xlabel('Bin ID')
    axes[1].set_ylabel('Total Size (Mbp)')
    axes[1].set_title('Bin Size Distribution')

    plt.tight_layout()
    plt.savefig(data_dir / 'bin_size_distribution.png', dpi=300)
    print(f"✅ Saved: bin_size_distribution.png")

def plot_coverage_gc_scatter(data_dir):
    """Plot coverage vs GC content colored by bin"""
    df = pd.read_csv(data_dir / "viz_coverage_gc.tsv", sep='\t')

    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(
        df['gc_content'] * 100,
        df['coverage'],
        c=df['bin_id'],
        s=df['length'] / 100,
        alpha=0.6,
        cmap='tab20'
    )

    plt.xlabel('GC Content (%)')
    plt.ylabel('Coverage (X)')
    plt.title('Coverage vs GC Content (colored by bin, sized by length)')
    plt.colorbar(scatter, label='Bin ID')

    plt.tight_layout()
    plt.savefig(data_dir / 'coverage_gc_scatter.png', dpi=300)
    print(f"✅ Saved: coverage_gc_scatter.png")

def plot_bin_quality_metrics(data_dir):
    """Plot bin quality metrics TSV"""
    df = pd.read_csv(data_dir / "bin_quality_metrics.tsv", sep='\t')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # N50
    axes[0, 0].bar(df['bin_id'], df['n50'] / 1000, color='green')
    axes[0, 0].set_xlabel('Bin ID')
    axes[0, 0].set_ylabel('N50 (kb)')
    axes[0, 0].set_title('Bin N50 Statistics')

    # Coverage
    axes[0, 1].bar(df['bin_id'], df['avg_coverage'], color='blue')
    axes[0, 1].set_xlabel('Bin ID')
    axes[0, 1].set_ylabel('Average Coverage')
    axes[0, 1].set_title('Bin Coverage Distribution')

    # GC Content
    axes[1, 0].bar(df['bin_id'], df['avg_gc_content'] * 100, color='orange')
    axes[1, 0].set_xlabel('Bin ID')
    axes[1, 0].set_ylabel('GC Content (%)')
    axes[1, 0].set_title('Bin GC Content')

    # Confidence
    axes[1, 1].bar(df['bin_id'], df['avg_confidence'], color='purple')
    axes[1, 1].set_xlabel('Bin ID')
    axes[1, 1].set_ylabel('Average Confidence')
    axes[1, 1].set_title('Bin Classification Confidence')
    axes[1, 1].set_ylim([0, 1])

    plt.tight_layout()
    plt.savefig(data_dir / 'bin_quality_metrics.png', dpi=300)
    print(f"✅ Saved: bin_quality_metrics.png")

if __name__ == "__main__":
    # Find most recent run directory
    import sys
    if len(sys.argv) > 1:
        data_dir = Path(sys.argv[1])
    else:
        # Auto-detect latest run
        output_dir = Path("output")
        runs = sorted(output_dir.glob("run_*"))
        if not runs:
            print("❌ No run directories found")
            sys.exit(1)
        data_dir = runs[-1] / "classification"

    print(f"📊 Generating plots from: {data_dir}")

    plot_bin_size_distribution(data_dir)
    plot_coverage_gc_scatter(data_dir)
    plot_bin_quality_metrics(data_dir)

    print("✅ All plots generated successfully!")
"#;

        let script_path = self
            .output_manager
            .get_section_dir(&PipelineSection::Classification)
            .join("plot_results.py");

        std::fs::write(&script_path, script)
            .with_context(|| format!("Failed to write Python script: {}", script_path.display()))?;

        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perms = std::fs::metadata(&script_path)?.permissions();
            perms.set_mode(0o755);
            std::fs::set_permissions(&script_path, perms)?;
        }

        tracing::info!("📜 Saved Python plotting script: plot_results.py");

        Ok(())
    }

    /// Export R plotting script
    fn export_r_plotting_script(&self) -> Result<()> {
        let script = r#"#!/usr/bin/env Rscript
# Visualization script for MetaForge classification results
# Requires: ggplot2, dplyr, readr

library(ggplot2)
library(dplyr)
library(readr)

# Find data directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  data_dir <- args[1]
} else {
  # Auto-detect latest run
  runs <- list.dirs("output", recursive = FALSE)
  runs <- runs[grepl("run_", basename(runs))]
  if (length(runs) == 0) {
    stop("❌ No run directories found")
  }
  data_dir <- file.path(tail(sort(runs), 1), "classification")
}

cat(sprintf("📊 Generating plots from: %s\n", data_dir))

# Load data
bin_sizes <- read_tsv(file.path(data_dir, "viz_bin_sizes.tsv"))
coverage_gc <- read_tsv(file.path(data_dir, "viz_coverage_gc.tsv"))
quality_metrics <- read_tsv(file.path(data_dir, "bin_quality_metrics.tsv"))

# Plot 1: Bin size distribution
p1 <- ggplot(bin_sizes, aes(x = factor(bin_id), y = contig_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Bin ID", y = "Number of Contigs", title = "Contigs per Bin") +
  theme_minimal()

ggsave(file.path(data_dir, "bin_sizes_r.png"), p1, width = 10, height = 6, dpi = 300)
cat("✅ Saved: bin_sizes_r.png\n")

# Plot 2: Coverage vs GC scatter
p2 <- ggplot(coverage_gc, aes(x = gc_content * 100, y = coverage, color = factor(bin_id))) +
  geom_point(aes(size = length), alpha = 0.6) +
  labs(x = "GC Content (%)", y = "Coverage (X)",
       title = "Coverage vs GC Content by Bin",
       color = "Bin ID", size = "Length") +
  theme_minimal()

ggsave(file.path(data_dir, "coverage_gc_r.png"), p2, width = 12, height = 8, dpi = 300)
cat("✅ Saved: coverage_gc_r.png\n")

cat("✅ All R plots generated successfully!\n")
"#;

        let script_path = self
            .output_manager
            .get_section_dir(&PipelineSection::Classification)
            .join("plot_results.R");

        std::fs::write(&script_path, script)
            .with_context(|| format!("Failed to write R script: {}", script_path.display()))?;

        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mut perms = std::fs::metadata(&script_path)?.permissions();
            perms.set_mode(0o755);
            std::fs::set_permissions(&script_path, perms)?;
        }

        tracing::info!("📜 Saved R plotting script: plot_results.R");

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_point_creation() {
        let point = ContigPoint {
            contig_id: 1,
            bin_id: 0,
            coverage: 15.5,
            gc_content: 0.52,
            length: 1500,
            confidence: 0.9,
        };

        assert_eq!(point.contig_id, 1);
        assert_eq!(point.bin_id, 0);
        assert!((point.coverage - 15.5).abs() < 0.01);
    }
}
