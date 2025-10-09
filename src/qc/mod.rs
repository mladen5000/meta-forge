//! Quality Control Module
//!
//! Provides comprehensive QC for metagenomic reads:
//! - Quality-based trimming (sliding window)
//! - Length filtering
//! - Adapter detection and removal
//! - QC statistics and reporting
//!
//! # Examples
//!
//! ## Basic Usage
//!
//! ```rust,ignore
//! use meta_forge::qc::{QCPipeline, QCPipelineConfig};
//!
//! let config = QCPipelineConfig::default();
//! let mut pipeline = QCPipeline::new(config);
//! let filtered = pipeline.process_reads(&raw_reads);
//! ```
//!
//! ## Iterator-based Processing
//!
//! ```rust,ignore
//! use meta_forge::qc::iter::QCIterExt;
//!
//! let filtered: Vec<_> = raw_reads
//!     .into_iter()
//!     .qc_default()
//!     .collect();
//! ```
//!
//! ## Using Presets
//!
//! ```rust,ignore
//! use meta_forge::qc::{QCPipeline, QCPreset};
//!
//! // Use preset configurations
//! let config = QCPreset::Strict.into();
//! let mut pipeline = QCPipeline::new(config);
//! let filtered = pipeline.process_reads(&raw_reads);
//! ```

pub mod quality_filter;
pub mod adapter_trimmer;
pub mod qc_stats;
pub mod qc_pipeline;
pub mod iter;
pub mod presets;
pub mod preprocessing;

pub use quality_filter::{QualityFilter, QualityFilterConfig};
pub use adapter_trimmer::{AdapterTrimmer, AdapterConfig, AdapterMatch};
pub use qc_stats::{QCStats, QCReport, FailureReason};
pub use qc_pipeline::{QCPipeline, QCPipelineConfig};
pub use iter::{QCFilterIter, QCIterExt};
pub use presets::QCPreset;
pub use preprocessing::{Preprocessor, FileFormat};
