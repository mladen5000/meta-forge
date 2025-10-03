//! Quality Control Module
//!
//! Provides comprehensive QC for metagenomic reads:
//! - Quality-based trimming (sliding window)
//! - Length filtering
//! - Adapter detection and removal
//! - Complexity filtering
//! - QC statistics and reporting

pub mod quality_filter;
pub mod adapter_trimmer;
pub mod qc_stats;
pub mod qc_pipeline;

pub use quality_filter::{QualityFilter, QualityFilterConfig};
pub use adapter_trimmer::{AdapterTrimmer, AdapterConfig};
pub use qc_stats::{QCStats, QCReport};
pub use qc_pipeline::{QCPipeline, QCPipelineConfig};
