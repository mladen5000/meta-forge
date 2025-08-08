//! Simple test of the optimized assembler to verify it generates contigs
//! ===================================================================

use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::core::data_structures::CorrectedRead;

fn main() -> anyhow::Result<()> {
    println!("🧬 Testing Optimized Assembly Graph Construction");
    
    // Create test reads with overlapping sequences - longer sequences for proper k-mer extraction
    let test_reads = vec![
        CorrectedRead {
            id: 0,
            original: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
            corrected: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 40],
        },
        CorrectedRead {
            id: 1,
            original: "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA".to_string(),
            corrected: "TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 40],
        },
        CorrectedRead {
            id: 2,
            original: "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT".to_string(),
            corrected: "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 40],
        },
        CorrectedRead {
            id: 3,
            original: "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATC".to_string(),
            corrected: "GATCGATCGATCGATCGATCGATCGATCGATCGATCGATC".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 40],
        },
        CorrectedRead {
            id: 4,
            original: "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA".to_string(),
            corrected: "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 40],
        },
        CorrectedRead {
            id: 5,
            original: "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG".to_string(),
            corrected: "CTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 40],
        },
    ];
    
    println!("📊 Test data: {} reads", test_reads.len());
    
    // Create assembly graph builder with smaller k-mer sizes
    let builder = AssemblyGraphBuilder::new(
        11,  // base k-mer size
        15,  // max k-mer size  
        1    // minimum coverage
    );
    
    println!("🏗️  Building assembly graph...");
    let mut graph = builder.build_graph(&test_reads)?;
    
    println!("📈 Graph statistics:");
    println!("   Nodes: {}", graph.graph_fragment.nodes.len());
    println!("   Edges: {}", graph.graph_fragment.edges.len());
    
    if graph.graph_fragment.nodes.is_empty() {
        println!("❌ No nodes in graph - check k-mer extraction!");
        return Ok(());
    }
    
    if graph.graph_fragment.edges.is_empty() {
        println!("⚠️  No edges in graph - will create singleton contigs");
    }
    
    println!("🧬 Generating contigs...");
    graph.generate_contigs()?;
    
    println!("📊 Assembly results:");
    println!("   Contigs generated: {}", graph.contigs.len());
    
    if graph.contigs.is_empty() {
        println!("❌ NO CONTIGS GENERATED - This indicates the bug is still present!");
        
        // Debug information
        println!("\n🔍 Debug information:");
        println!("   Graph fragment nodes: {}", graph.graph_fragment.nodes.len());
        println!("   Graph fragment edges: {}", graph.graph_fragment.edges.len());
        println!("   Petgraph nodes: {}", graph.petgraph.node_count());
        println!("   Petgraph edges: {}", graph.petgraph.edge_count());
        
        // Print sample node information
        if !graph.graph_fragment.nodes.is_empty() {
            println!("   Sample nodes:");
            for (i, (hash, node)) in graph.graph_fragment.nodes.iter().take(5).enumerate() {
                println!("     {}. Hash: {}, Sequence: {}, Coverage: {}", 
                         i + 1, hash, node.kmer.sequence, node.coverage);
            }
        }
        
        return Ok(());
    }
    
    println!("✅ SUCCESS! Contigs generated:");
    for (i, contig) in graph.contigs.iter().enumerate() {
        println!("   Contig {}: {} bp, coverage: {:.2}, type: {:?}", 
                 i + 1, contig.length, contig.coverage, contig.contig_type);
        if contig.sequence.len() <= 100 {
            println!("     Sequence: {}", contig.sequence);
        } else {
            println!("     Sequence: {}...{} (truncated)", 
                     &contig.sequence[..50], 
                     &contig.sequence[contig.sequence.len()-50..]);
        }
    }
    
    println!("📊 Assembly statistics:");
    println!("   Total contigs: {}", graph.assembly_stats.total_contigs);
    println!("   Total length: {} bp", graph.assembly_stats.total_length);
    println!("   Longest contig: {} bp", graph.assembly_stats.longest_contig);
    println!("   N50: {} bp", graph.assembly_stats.n50);
    println!("   Mean coverage: {:.2}", graph.assembly_stats.mean_coverage);
    println!("   GC content: {:.3}", graph.assembly_stats.gc_content);
    
    println!("\n🎉 Assembler test completed successfully!");
    
    Ok(())
}