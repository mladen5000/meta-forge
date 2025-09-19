// Essential assembly modules for laptop-optimized processing
pub mod laptop_assembly;
pub mod graph_construction; // Basic graph construction utilities
pub mod adaptive_k; // Adaptive k-mer selection

// Re-export main laptop assembler for convenient access
pub use laptop_assembly::{LaptopAssembler, LaptopConfig, LaptopAssemblyGraph};
