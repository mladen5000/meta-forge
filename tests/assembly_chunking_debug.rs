//! Debug test for the chunking and parallel fragment construction process
//! This will help identify why chunks have graph fragments but the final merge results in 0 nodes.

use meta_forge::core::data_structures::*;
use meta_forge::assembly::graph_construction::AdvancedAssemblyGraphBuilder;

#[test]
fn test_chunking_and_merging() {
    println!("🧪 Testing chunking and merging process");
    
    // Create test reads - same as our CLI test
    let test_reads = vec![
        CorrectedRead {
            id: 0,
            original: "ATCGATCGATCGATCGATCG".to_string(),
            corrected: "ATCGATCGATCGATCGATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 20],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 1,
            original: "TCGATCGATCGATCGATCGA".to_string(),
            corrected: "TCGATCGATCGATCGATCGA".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 20],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 2,
            original: "CGATCGATCGATCGATCGAT".to_string(),
            corrected: "CGATCGATCGATCGATCGAT".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 20],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
    ];
    
    println!("📊 Test data: {} reads", test_reads.len());
    for (i, read) in test_reads.iter().enumerate() {
        println!("  Read {}: '{}'", i+1, read.corrected);
    }
    
    // Create builder with same parameters as CLI
    let builder = match AdvancedAssemblyGraphBuilder::new(15, 31, 2, 8) {
        Ok(b) => b,
        Err(e) => {
            println!("❌ Failed to create builder: {}", e);
            panic!("Builder creation failed");
        }
    };
    
    println!("\n🔧 Testing graph construction...");
    
    // Try to build the graph using the same method as the CLI
    match builder.build_graph(&test_reads) {
        Ok(assembly_graph) => {
            println!("✅ Graph construction succeeded!");
            println!("📊 Final results:");
            println!("  Graph fragment nodes: {}", assembly_graph.graph_fragment.nodes.len());
            println!("  Graph fragment edges: {}", assembly_graph.graph_fragment.edges.len());
            println!("  Contigs: {}", assembly_graph.contigs.len());
            println!("  Assembly stats - total length: {}", assembly_graph.assembly_stats.total_length);
            
            if assembly_graph.graph_fragment.nodes.is_empty() {
                println!("🚨 BUG CONFIRMED: Graph construction produces empty GraphFragment!");
                println!("This is the root cause of zero contigs.");
                
                // The issue might be in coverage filtering - let's check if nodes exist but get filtered out
                println!("\nDEBUG: This could be caused by:");
                println!("1. Coverage filtering (min_coverage=2) removing all nodes");
                println!("2. Hierarchical merging losing graph fragments");
                println!("3. Parallel processing issues");
                println!("4. Graph optimization removing all nodes");
            } else {
                println!("✅ Graph construction working correctly");
                
                // Print detailed node information
                for (hash, node) in &assembly_graph.graph_fragment.nodes {
                    println!("  Node {}: k-mer='{}', coverage={}", 
                            hash, node.kmer.sequence, node.coverage);
                }
                
                // Print detailed edge information
                for (i, edge) in assembly_graph.graph_fragment.edges.iter().enumerate() {
                    println!("  Edge {}: {} -> {} (weight={})", 
                            i+1, edge.from_hash, edge.to_hash, edge.weight);
                }
                
                // Print contigs
                for (i, contig) in assembly_graph.contigs.iter().enumerate() {
                    println!("  Contig {}: {}bp, sequence='{}'", 
                            i+1, contig.length, contig.sequence);
                }
            }
        }
        Err(e) => {
            println!("❌ Graph construction failed: {}", e);
            panic!("Graph construction error: {}", e);
        }
    }
}