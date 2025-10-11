pub mod data_structures;
pub mod paired_reads;
pub mod pipeline_types;

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

// Re-export pipeline types for easy access
pub use pipeline_types::{
    AbundanceProfile, AnalysisReport, AnalysisResults, AssemblyResults, FeatureCollection,
    FileFormat, PerformanceMetrics, QualityMetrics, ReportSummary, TaxonomicClassification,
};
