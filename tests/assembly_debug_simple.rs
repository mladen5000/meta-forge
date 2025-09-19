//! Simple Assembly Debug Test
//!
//! This test focuses on the core issue: why are zero contigs being generated?
//! We'll test each step individually to identify the failure point.

use meta_forge::assembly::performance_optimizations::*;
use meta_forge::core::data_structures::*;

/// Test the ParallelContigGenerator directly with minimal data
#[test]
fn test_parallel_contig_generator_basic() {
    println!("🧪 Testing ParallelContigGenerator with minimal setup");

    // Create a very simple graph with just a few connected nodes
    let mut graph = CacheOptimizedGraph::new(3);

    // Add 3 nodes with different hashes
    let node1 = graph.add_node(100, 5); // hash=100, coverage=5
    let node2 = graph.add_node(200, 4); // hash=200, coverage=4
    let node3 = graph.add_node(300, 3); // hash=300, coverage=3

    println!("Added 3 nodes: {} {} {}", node1, node2, node3);

    // Add edges to connect them: node1 -> node2 -> node3
    match graph.add_edge(100, 200) {
        Ok(()) => println!("✅ Added edge 100 -> 200"),
        Err(e) => println!("❌ Failed to add edge 100 -> 200: {}", e),
    }

    match graph.add_edge(200, 300) {
        Ok(()) => println!("✅ Added edge 200 -> 300"),
        Err(e) => println!("❌ Failed to add edge 200 -> 300: {}", e),
    }

    // Check graph statistics
    let (nodes, edges, memory, cache_rate) = graph.get_statistics();
    println!(
        "📊 Graph stats: {} nodes, {} edges, {}KB memory, {:.1}% cache hit",
        nodes,
        edges,
        memory / 1024,
        cache_rate * 100.0
    );

    // Now try to generate contigs
    println!("🧬 Generating contigs...");
    match ParallelContigGenerator::generate_contigs_parallel(&graph) {
        Ok(contigs) => {
            println!("✅ SUCCESS! Generated {} contigs:", contigs.len());
            for (i, contig) in contigs.iter().enumerate() {
                println!(
                    "  Contig {}: {}bp, coverage={:.1}x, sequence='{}'",
                    i + 1,
                    contig.length,
                    contig.coverage,
                    if contig.sequence.len() > 50 {
                        format!("{}...", &contig.sequence[..47])
                    } else {
                        contig.sequence.clone()
                    }
                );
            }

            // CRITICAL: Check for the bug
            if contigs.is_empty() {
                println!(
                    "🚨 BUG DETECTED: Zero contigs generated from {} connected nodes",
                    nodes
                );
                panic!("Zero contigs bug confirmed!");
            }
        }
        Err(e) => {
            println!("❌ FAILED to generate contigs: {}", e);
            panic!("Contig generation completely failed: {}", e);
        }
    }
}

/// Test with unconnected nodes (should generate multiple single-node contigs)
#[test]
fn test_parallel_contig_generator_unconnected() {
    println!("🧪 Testing ParallelContigGenerator with unconnected nodes");

    let mut graph = CacheOptimizedGraph::new(3);

    // Add 3 unconnected nodes
    graph.add_node(111, 10);
    graph.add_node(222, 8);
    graph.add_node(333, 6);

    // Do NOT add any edges - nodes are isolated

    let (nodes, edges, _, _) = graph.get_statistics();
    println!("📊 Unconnected graph: {} nodes, {} edges", nodes, edges);

    match ParallelContigGenerator::generate_contigs_parallel(&graph) {
        Ok(contigs) => {
            println!(
                "✅ Generated {} contigs from {} unconnected nodes",
                contigs.len(),
                nodes
            );

            // Each isolated node should become its own contig
            assert_eq!(
                contigs.len(),
                nodes,
                "Unconnected nodes should produce {} contigs, got {}",
                nodes,
                contigs.len()
            );

            for (i, contig) in contigs.iter().enumerate() {
                println!(
                    "  Isolated contig {}: {}bp, sequence='{}'",
                    i + 1,
                    contig.length,
                    contig.sequence
                );
            }
        }
        Err(e) => {
            println!("❌ Failed on unconnected nodes: {}", e);
        }
    }
}

/// Test k-mer extraction and canonical representation
#[test]
fn test_kmer_extraction_basics() {
    println!("🧪 Testing basic k-mer extraction");

    let sequence = "ATCGATCG"; // 8bp sequence
    let k = 4;

    println!("Input sequence: '{}' ({}bp)", sequence, sequence.len());
    println!("k-mer size: {}", k);

    // Manual k-mer extraction
    let expected_kmers = vec!["ATCG", "TCGA", "CGAT", "GATC", "ATCG"];
    println!("Expected k-mers: {:?}", expected_kmers);

    // Test canonical k-mer creation
    let mut canonical_kmers = Vec::new();
    for (i, expected) in expected_kmers.iter().enumerate() {
        match CanonicalKmer::new(expected) {
            Ok(canonical) => {
                println!(
                    "  K-mer {}: '{}' -> canonical '{}' (hash={})",
                    i + 1,
                    expected,
                    canonical.sequence,
                    canonical.hash
                );
                canonical_kmers.push(canonical);
            }
            Err(e) => {
                println!("  ❌ Failed to create k-mer '{}': {}", expected, e);
                panic!("K-mer creation failed");
            }
        }
    }

    assert_eq!(canonical_kmers.len(), expected_kmers.len());
    println!("✅ K-mer extraction successful");
}

/// Test AssemblyChunk processing (the critical step)
#[test]
fn test_assembly_chunk_processing() {
    println!("🧪 Testing AssemblyChunk processing");

    // Create a simple read
    let read = CorrectedRead {
        id: 0,
        original: "ATCGATCGATCG".to_string(),
        corrected: "ATCGATCGATCG".to_string(),
        corrections: Vec::new(),
        quality_scores: vec![30; 12],
        correction_metadata: CorrectionMetadata {
            algorithm: "test".to_string(),
            confidence_threshold: 0.95,
            context_window: 3,
            correction_time_ms: 0,
        },
    };

    let k = 4;
    println!("Processing read: '{}' with k={}", read.corrected, k);

    // Create AssemblyChunk and process the read
    let mut chunk = meta_forge::assembly::graph_construction::AssemblyChunk::new(0, k);

    match chunk.add_read(read) {
        Ok(()) => println!("✅ Read added to chunk"),
        Err(e) => {
            println!("❌ Failed to add read: {}", e);
            panic!("Read addition failed");
        }
    }

    // Finalize the chunk (this should build the graph)
    println!("Finalizing chunk...");
    chunk.finalize();

    // Check the resulting graph fragment
    let nodes = chunk.graph_fragment.nodes.len();
    let edges = chunk.graph_fragment.edges.len();

    println!(
        "📊 GraphFragment after finalization: {} nodes, {} edges",
        nodes, edges
    );

    if nodes == 0 {
        println!("🚨 BUG FOUND: No nodes created from read processing!");
        println!("This indicates the AssemblyChunk.finalize() is not working correctly");
        panic!("Zero nodes created from read - this is the bug!");
    }

    if edges == 0 {
        println!("🚨 BUG FOUND: No edges created from read processing!");
        println!("Nodes exist but are not connected - this causes zero contigs");
        panic!("Zero edges created from read - nodes are isolated!");
    }

    // Print detailed node information
    for (hash, node) in &chunk.graph_fragment.nodes {
        println!(
            "  Node {}: coverage={}, k-mer='{}'",
            hash, node.coverage, node.kmer.sequence
        );
    }

    // Print edge information
    for (i, edge) in chunk.graph_fragment.edges.iter().enumerate() {
        println!(
            "  Edge {}: {} -> {} (weight={})",
            i + 1,
            edge.from_hash,
            edge.to_hash,
            edge.weight
        );
    }

    println!("✅ AssemblyChunk processing completed successfully");
}

/// Integration test: Read -> Chunk -> Graph -> Contigs
#[test]
fn test_full_pipeline_simple() {
    println!("🧪 Testing full pipeline with simple data");

    // Create overlapping reads
    let reads = vec![
        CorrectedRead {
            id: 0,
            original: "ATCGATCG".to_string(),
            corrected: "ATCGATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 8],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 1,
            original: "TCGATCGA".to_string(), // Overlaps by 7bp
            corrected: "TCGATCGA".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 8],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
    ];

    let k = 4;
    println!("Processing {} reads with k={}", reads.len(), k);

    // Step 1: Process reads into AssemblyChunk
    let mut chunk = meta_forge::assembly::graph_construction::AssemblyChunk::new(0, k);

    for (i, read) in reads.iter().enumerate() {
        println!("Adding read {}: '{}'", i + 1, read.corrected);
        chunk
            .add_read(read.clone())
            .expect("Read addition should succeed");
    }

    chunk.finalize();

    let nodes = chunk.graph_fragment.nodes.len();
    let edges = chunk.graph_fragment.edges.len();
    println!(
        "📊 After chunk processing: {} nodes, {} edges",
        nodes, edges
    );

    if nodes == 0 {
        panic!(
            "🚨 PIPELINE FAILURE: No nodes created from {} reads",
            reads.len()
        );
    }

    // Step 2: Convert to CacheOptimizedGraph
    let mut cache_graph = CacheOptimizedGraph::new(nodes);

    for node in chunk.graph_fragment.nodes.values() {
        cache_graph.add_node(node.kmer.hash, node.coverage);
    }

    for edge in &chunk.graph_fragment.edges {
        match cache_graph.add_edge(edge.from_hash, edge.to_hash) {
            Ok(()) => {}
            Err(e) => println!("Warning: Failed to add edge: {}", e),
        }
    }

    let (final_nodes, final_edges, _, _) = cache_graph.get_statistics();
    println!(
        "📊 CacheOptimizedGraph: {} nodes, {} edges",
        final_nodes, final_edges
    );

    // Step 3: Generate contigs
    println!("🧬 Generating contigs...");
    match ParallelContigGenerator::generate_contigs_parallel(&cache_graph) {
        Ok(contigs) => {
            println!(
                "✅ PIPELINE SUCCESS! Generated {} contigs from {} reads",
                contigs.len(),
                reads.len()
            );

            if contigs.is_empty() {
                println!(
                    "🚨 CRITICAL BUG: Zero contigs generated despite {} nodes and {} edges",
                    final_nodes, final_edges
                );
                panic!("Zero contig bug confirmed in full pipeline!");
            }

            // Verify contigs
            for (i, contig) in contigs.iter().enumerate() {
                println!(
                    "  Contig {}: {}bp, coverage={:.1}x, sequence='{}'",
                    i + 1,
                    contig.length,
                    contig.coverage,
                    contig.sequence
                );

                assert!(contig.length > 0, "Contig should have positive length");
                assert!(!contig.sequence.is_empty(), "Contig should have sequence");
            }

            // Success metric: fewer contigs than reads
            assert!(
                contigs.len() < reads.len(),
                "Should produce fewer contigs ({}) than reads ({})",
                contigs.len(),
                reads.len()
            );
        }
        Err(e) => {
            println!("🚨 PIPELINE FAILURE at contig generation: {}", e);
            panic!("Contig generation failed: {}", e);
        }
    }
}
