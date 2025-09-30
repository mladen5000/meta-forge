//! Optimized Assembly Architecture
//! =============================
//!
//! High-performance assembly implementation using advanced data structures
//! and algorithms optimized for speed, memory efficiency, and correctness.
//!
//! This module implements the architectural optimizations outlined in the
//! assembly_architecture_optimization.md document.

pub mod bit_packed_kmer;
pub mod csr_graph;
pub mod hyperlog_counter;
pub mod streaming_pipeline;
pub mod memory_pool;
pub mod resource_manager;

// Temporarily disabled - requires SIMD/architecture-specific support
// pub mod zero_copy_kmer;
// pub mod fast_memory_pool;
// pub mod fast_contig_builder;

// Re-export main optimized components
pub use bit_packed_kmer::{BitPackedKmer, SIMDKmerExtractor};
pub use csr_graph::{CSRAssemblyGraph, NeighborIterator};
pub use hyperlog_counter::{HyperLogKmerCounter, CountMinSketch};
pub use streaming_pipeline::{StreamingAssemblyPipeline, PipelineStage};
pub use memory_pool::{AssemblyMemoryPool, PooledKmer};
pub use resource_manager::{AdaptiveResourceManager, SystemMonitor};

// Disabled temporarily
// pub use zero_copy_kmer::{ZeroCopyKmerIterator, RollingHashKmerCounter};
// pub use fast_memory_pool::{FastAssemblyMemoryPool, PoolConfig};
// pub use fast_contig_builder::{FastContigBuilder, ContigBuilderConfig};

/// Optimized assembler using new architecture
pub mod optimized_assembler;
pub use optimized_assembler::OptimizedAssembler;