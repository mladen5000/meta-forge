//! Assembly modules for laptop-optimized processing
//!
//! **Production Code**: `laptop_assembly.rs` - Optimized assembler with all features integrated
//! **Experimental**: `optimized/` - Alternative implementations for benchmarking

// Production assembler
pub mod laptop_assembly;
pub mod adaptive_k; // Adaptive k-mer selection

// Experimental optimized architecture (for benchmarking only)
pub mod optimized;

// Re-export production components
pub use laptop_assembly::{
    LaptopAssembler,
    LaptopConfig,
    LaptopAssemblyGraph,
    CompactKmer,
    RollingKmerHash,
};

// Re-export experimental components (OptimizedAssembler disabled - needs refactoring)
// pub use optimized::OptimizedAssembler;
