//! Optimized Assembly Components
//! ============================
//!
//! Experimental data structures for benchmarking and research.
//!
//! **Production assembler**: Use `laptop_assembly::LaptopAssembler` instead.
//! This module contains experimental alternatives not used in the main pipeline.

pub mod csr_graph;
pub mod streaming_pipeline;
pub mod resource_manager;

// Re-export experimental components (for benchmarking/research only)
pub use csr_graph::{CSRAssemblyGraph, NeighborIterator};
pub use streaming_pipeline::{StreamingAssemblyPipeline, PipelineStage};
pub use resource_manager::{AdaptiveResourceManager, SystemMonitor};