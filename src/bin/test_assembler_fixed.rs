//! Simple test of the optimized assembler to verify it generates contigs
//! ===================================================================

use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};

/// Test the assembly graph builder with corrected reads
fn main() -> anyhow::Result<()> {
    println!("🧬 Testing Optimized Assembly Graph Construction");
    
    // Create test reads with overlapping sequences - longer sequences for proper k-mer extraction
    let test_reads = vec![
        CorrectedRead {
            id: 0,
            original: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
            corrected: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
            corrections: vec![],
            quality_scores: vec![40; 40],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 1,
            original: "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT".to_string(),
            corrected: "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT".to_string(),
            corrections: vec![],
            quality_scores: vec![40; 40],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 2,
            original: "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATC".to_string(),
            corrected: "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATC".to_string(),
            corrections: vec![],
            quality_scores: vec![40; 40],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        },
    ];

    println!("📊 Created {} test reads", test_reads.len());

    // Build the assembly graph
    let builder = AssemblyGraphBuilder::new(
        11,  // base k-mer size
        15,  // max k-mer size  
        1    // minimum coverage
    );

    println!("🔧 Building assembly graph...");
    let graph = builder.build(&test_reads)?;

    println!("✅ Assembly graph built successfully!");
    println!("   - Nodes: {}", graph.graph_fragment.nodes.len());
    println!("   - Edges: {}", graph.graph_fragment.edges.len());
    println!("   - Contigs: {}", graph.contigs.len());

    if !graph.contigs.is_empty() {
        println!("🧬 First contig sequence: {}", graph.contigs[0].sequence);
    }

    Ok(())
}