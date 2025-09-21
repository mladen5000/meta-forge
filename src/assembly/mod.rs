//! Assembly modules for laptop-optimized processing
//! Includes both the original laptop-optimized implementation and new architectural optimizations

// Original laptop-optimized implementation
pub mod laptop_assembly;
pub mod graph_construction; // Basic graph construction utilities
pub mod adaptive_k; // Adaptive k-mer selection

// New optimized architecture (experimental/high-performance)
pub mod optimized; // Re-enabling optimized modules

// Re-export main components for convenient access
pub use laptop_assembly::{LaptopAssembler, LaptopConfig, LaptopAssemblyGraph};

// Re-export optimized components for advanced users
pub use optimized::OptimizedAssembler;
