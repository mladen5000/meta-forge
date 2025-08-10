pub mod data_structures;
pub mod paired_reads;

// Re-export key types for assembly integration
pub use data_structures::{
    // Processing types
    AssemblyChunk,
    // Assembly-specific types
    AssemblyGraph,
    AssemblyStats,
    BaseCorrection,
    // Graph analysis types
    BubbleStructure,
    BubbleType,
    // Core k-mer and graph structures
    CanonicalKmer,
    Contig,
    ContigType,
    // Read and correction types
    CorrectedRead,
    CorrectionMetadata,
    CorrectionType,

    CoverageStats,
    EdgeType,

    EdgeUpdate,
    EdgeWeight,
    EulerianPathType,

    GraphEdge,
    GraphFragment,

    GraphNode,
    GraphUpdate,

    // Utility types
    Minimizer,
    MinimizerExtractor,
    NodeType,
    NodeUpdate,
    ProcessingStats,
    ReadPosition,
    Strand,
};
