pub mod data_structures;
pub mod paired_reads;

// Re-export key types for assembly integration
pub use data_structures::{
    // Core k-mer and graph structures
    CanonicalKmer, GraphNode, GraphEdge, GraphFragment,
    
    // Assembly-specific types
    AssemblyGraph, AssemblyStats, Contig, ContigType,
    EdgeWeight, EulerianPathType,
    
    // Read and correction types
    CorrectedRead, CorrectionMetadata, BaseCorrection, CorrectionType,
    
    // Graph analysis types
    BubbleStructure, BubbleType, NodeType, EdgeType,
    
    // Processing types
    AssemblyChunk, ProcessingStats, GraphUpdate,
    
    // Utility types
    Minimizer, MinimizerExtractor, ReadPosition, Strand,
    CoverageStats, NodeUpdate, EdgeUpdate,
};