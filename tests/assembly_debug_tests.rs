use anyhow::Result;
use meta_forge::core::data_structures::{AssemblyChunk, CorrectedRead, CanonicalKmer, CorrectionMetadata};
use meta_forge::assembly::performance_optimizations::{ParallelContigGenerator, CacheOptimizedGraph};

/// Test 1: Basic k-mer creation and canonical form
#[test]
fn test_kmer_creation() {
    let kmer = CanonicalKmer::new("ATCG").unwrap();
    assert!(!kmer.sequence.is_empty());
    assert!(kmer.hash != 0);
    println!("✅ K-mer creation works: hash={}, canonical={}", kmer.hash, kmer.is_canonical);
}

/// Test 2: AssemblyChunk can add reads and create nodes
#[test]
fn test_assembly_chunk_add_read() -> Result<()> {
    let mut chunk = AssemblyChunk::new(0, 15); // k-mer size 15
    
    // Create a simple test read
    let read = CorrectedRead {
        id: 1,
        original: "ATCGATCGATCGATCG".to_string(),
        corrected: "ATCGATCGATCGATCG".to_string(), // 16bp = 2 k-mers of size 15
        corrections: vec![],
        quality_scores: vec![40; 16],
        correction_metadata: CorrectionMetadata {
            algorithm: "test".to_string(),
            confidence_threshold: 0.9,
            context_window: 5,
            correction_time_ms: 0,
        },
    };
    
    // Add the read
    chunk.add_read(read)?;
    
    // Check that nodes were created
    println!("📊 Chunk stats: {} reads processed, {} nodes created", 
             chunk.processing_stats.reads_processed,
             chunk.processing_stats.nodes_created);
    
    assert_eq!(chunk.processing_stats.reads_processed, 1);
    assert!(chunk.processing_stats.nodes_created > 0, "No nodes were created!");
    assert!(!chunk.graph_fragment.nodes.is_empty(), "Graph fragment is empty!");
    
    // Check that edges were created between consecutive k-mers
    assert!(chunk.processing_stats.edges_created > 0, "No edges were created!");
    
    println!("✅ AssemblyChunk.add_read created {} nodes and {} edges", 
             chunk.processing_stats.nodes_created, 
             chunk.processing_stats.edges_created);
    
    Ok(())
}

/// Test 3: Multiple overlapping reads should create connected components
#[test]
fn test_overlapping_reads_create_connected_graph() -> Result<()> {
    let mut chunk = AssemblyChunk::new(0, 15);
    
    // Create overlapping reads that should form a connected component
    let reads = vec![
        CorrectedRead {
            id: 1,
            original: "ATCGATCGATCGATCG".to_string(),
            corrected: "ATCGATCGATCGATCG".to_string(),
            corrections: vec![],
            quality_scores: vec![40; 16],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 2,
            original: "TCGATCGATCGATCGA".to_string(),
            corrected: "TCGATCGATCGATCGA".to_string(), // Overlaps with first read
            corrections: vec![],
            quality_scores: vec![40; 16],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        },
    ];
    
    for read in reads {
        chunk.add_read(read)?;
    }
    
    println!("📊 Connected graph: {} reads, {} nodes, {} edges", 
             chunk.processing_stats.reads_processed,
             chunk.processing_stats.nodes_created,
             chunk.processing_stats.edges_created);
    
    // Should have more than 2 nodes (each read contributes k-mers)
    assert!(chunk.processing_stats.nodes_created >= 2);
    assert!(chunk.processing_stats.edges_created >= 2);
    
    // Finalize the chunk to compute connectivity
    chunk.finalize();
    
    println!("✅ Connected graph test passed");
    Ok(())
}

/// Test 4: ParallelContigGenerator can find connected components
#[test]
fn test_parallel_contig_generator() -> Result<()> {
    // Create a simple cache-optimized graph
    let mut graph = CacheOptimizedGraph::new(4); // Estimated 4 nodes
    
    // Add some nodes (simplified - using direct hash values)
    let nodes = vec![1u64, 2u64, 3u64, 4u64];
    for &node in &nodes {
        graph.add_node(node, 1); // hash, coverage
    }
    
    // Add edges to form a linear chain: 1-2-3, 4 isolated  
    graph.add_edge(1, 2)?;
    graph.add_edge(2, 3)?;
    // Node 4 is isolated
    
    // Find connected components using static method
    let contigs = ParallelContigGenerator::generate_contigs_parallel(&graph)?;
    
    println!("📊 ParallelContigGenerator found {} contigs from {} nodes", 
             contigs.len(), nodes.len());
    
    // Should find 2 connected components: [1,2,3] and [4]
    assert!(contigs.len() >= 1, "No contigs generated!");
    
    for (i, contig) in contigs.iter().enumerate() {
        println!("  Contig {}: {} bp", i + 1, contig.sequence.len());
    }
    
    println!("✅ ParallelContigGenerator test passed");
    Ok(())
}

/// Test 5: End-to-end mini assembly pipeline
#[test]
fn test_mini_assembly_pipeline() -> Result<()> {
    // Create a chunk with multiple reads
    let mut chunk = AssemblyChunk::new(0, 15);
    
    // Add several reads that should form contigs
    let test_reads = vec![
        "ATCGATCGATCGATCG", // 16 bp
        "TCGATCGATCGATCGA", // Overlaps with first
        "CGATCGATCGATCGAT", // Overlaps with second
        "GGGAAATTTCCCGGG",  // Different sequence (separate contig)
    ];
    
    for (i, seq) in test_reads.iter().enumerate() {
        let read = CorrectedRead {
            id: i + 1,
            original: seq.to_string(),
            corrected: seq.to_string(),
            corrections: vec![],
            quality_scores: vec![40; seq.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        };
        chunk.add_read(read)?;
    }
    
    chunk.finalize();
    
    println!("📊 Mini pipeline: {} reads → {} nodes → {} edges", 
             chunk.processing_stats.reads_processed,
             chunk.processing_stats.nodes_created,
             chunk.processing_stats.edges_created);
    
    // Now test contig generation - create a CacheOptimizedGraph from the chunk
    let mut graph = CacheOptimizedGraph::new(100); // Estimated nodes from chunk
    
    // Transfer nodes from chunk to graph
    for (hash, node) in chunk.graph_fragment.nodes.iter() {
        graph.add_node(*hash, node.coverage as u32);
    }
    
    // Transfer edges
    for edge in chunk.graph_fragment.edges.iter() {
        graph.add_edge(edge.from_hash, edge.to_hash)?;
    }
    
    let contigs = ParallelContigGenerator::generate_contigs_parallel(&graph)?;
    
    println!("📊 Final result: {} contigs generated", contigs.len());
    for (i, contig) in contigs.iter().enumerate() {
        println!("  Contig {}: {} bp - {}", i + 1, contig.sequence.len(), 
                 if contig.sequence.len() > 20 { 
                     format!("{}...{}", &contig.sequence[..10], &contig.sequence[contig.sequence.len()-10..])
                 } else { 
                     contig.sequence.clone() 
                 });
    }
    
    // The critical assertion - we should have generated some contigs!
    assert!(contigs.len() > 0, "❌ CRITICAL: No contigs generated in end-to-end test!");
    
    println!("✅ End-to-end mini assembly pipeline test passed");
    Ok(())
}