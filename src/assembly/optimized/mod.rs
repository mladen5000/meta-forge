//! Optimized Assembly Architecture
//! =============================
//!
//! Experimental high-performance assembly components.
//!
//! NOTE: Production code uses `laptop_assembly.rs` which has all optimizations
//! integrated. This module contains experimental alternatives for benchmarking.

pub mod csr_graph;
pub mod streaming_pipeline;
pub mod resource_manager;

// Re-export experimental components (for benchmarking only)
pub use csr_graph::{CSRAssemblyGraph, NeighborIterator};
pub use streaming_pipeline::{StreamingAssemblyPipeline, PipelineStage};
pub use resource_manager::{AdaptiveResourceManager, SystemMonitor};

// Experimental optimized assembler (DISABLED - needs refactoring to use CompactKmer)
// pub mod optimized_assembler;
// pub use optimized_assembler::OptimizedAssembler;