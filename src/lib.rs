//! # MetaForge - Metagenomic Analysis Pipeline
//!
//! A high-performance metagenomics analysis pipeline built in Rust.
//! Implements advanced algorithms for DNA sequence assembly, taxonomic classification,
//! and abundance estimation using machine learning and graph neural networks.

pub mod assembly;
pub mod core;
pub mod database;
pub mod features;
pub mod ml;
pub mod pipeline;
pub mod utils;

// Re-export commonly used types at crate level
pub use crate::core::data_structures::*;
pub use crate::pipeline::complete_integration::MetagenomicsPipeline;

/// Result type used throughout the crate
pub type Result<T> = anyhow::Result<T>;

/// Error type used throughout the crate
pub type Error = anyhow::Error;