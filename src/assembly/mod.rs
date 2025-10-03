//! Assembly modules for laptop-optimized processing
//!
//! **Production Code**: `laptop_assembly.rs` - Complete assembler with all optimizations
//! **Experimental**: `optimized/` - Research data structures for benchmarking

// Production assembler
pub mod laptop_assembly;
pub mod adaptive_k;

// Experimental data structures (not used in production pipeline)
pub mod optimized;

// Re-export production components
pub use laptop_assembly::{
    LaptopAssembler,
    LaptopConfig,
    LaptopAssemblyGraph,
    CompactKmer,
    RollingKmerHash,
};
