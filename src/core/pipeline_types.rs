//! Pipeline Data Types
//!
//! Shared data structures used across the metagenomics pipeline.
//! These types are used by all pipeline coordinators and represent
//! intermediate and final results of the analysis workflow.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::Duration;

use super::data_structures::{AssemblyStats, Contig, GraphFragment};
use crate::features::extraction::FeatureVector;
use crate::utils::configuration::PipelineConfiguration;

/// Abundance profile for k-mer frequency analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AbundanceProfile {
    pub unique_kmers: u64,
    pub abundant_kmers: HashMap<u64, f64>,
    pub total_kmers: u64,
}

/// Complete analysis results from the pipeline
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisResults {
    pub sample_name: String,
    pub assembly_results: AssemblyResults,
    pub classifications: Vec<TaxonomicClassification>,
    pub abundance_profile: AbundanceProfile,
    pub features: FeatureCollection,
    pub report: AnalysisReport,
    pub processing_time: Duration,
    pub config_used: PipelineConfiguration,
}

/// Results from the assembly phase
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyResults {
    pub contigs: Vec<Contig>,
    pub assembly_stats: AssemblyStats,
    pub graph_fragment: GraphFragment,
}

/// Taxonomic classification of a contig
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonomicClassification {
    pub contig_id: usize,
    pub taxonomy_id: u32,
    pub taxonomy_name: String,
    pub confidence: f64,
    pub lineage: String,
    pub method: String,
}

/// Collection of extracted features
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureCollection {
    pub sequence_features: HashMap<usize, FeatureVector>,
    pub graph_features: Option<ndarray::Array1<f64>>,
}

impl Default for FeatureCollection {
    fn default() -> Self {
        Self::new()
    }
}

impl FeatureCollection {
    pub fn new() -> Self {
        Self {
            sequence_features: HashMap::new(),
            graph_features: None,
        }
    }

    pub fn add_sequence_features(&mut self, contig_id: usize, features: FeatureVector) {
        self.sequence_features.insert(contig_id, features);
    }

    pub fn set_graph_features(&mut self, features: ndarray::Array1<f64>) {
        self.graph_features = Some(features);
    }
}

/// Final analysis report
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

/// Performance metrics from the pipeline run
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
    pub total_processing_time: Duration,
    pub peak_memory_usage: usize,
    pub reads_processed: u64,
    pub errors_corrected: u64,
    pub repeats_resolved: u64,
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

/// Quality metrics for the analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
    pub assembly_completeness: f64,
    pub classification_confidence: f64,
    pub coverage_uniformity: f64,
}

/// File format detection
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fastq,
    FastqGz,
    Fasta,
    FastaGz,
    Unknown,
}

impl FileFormat {
    /// Detect file format from path extension
    pub fn from_path(path: &std::path::Path) -> Self {
        let file_name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_lowercase();

        if file_name.ends_with(".fastq") || file_name.ends_with(".fq") {
            FileFormat::Fastq
        } else if file_name.ends_with(".fastq.gz") || file_name.ends_with(".fq.gz") {
            FileFormat::FastqGz
        } else if file_name.ends_with(".fasta") || file_name.ends_with(".fa") {
            FileFormat::Fasta
        } else if file_name.ends_with(".fasta.gz")
            || file_name.ends_with(".fa.gz")
            || file_name.ends_with(".fna.gz")
        {
            FileFormat::FastaGz
        } else {
            FileFormat::Unknown
        }
    }
}
