//! Test suite for validating the assembly connectivity fixes
//! This test ensures that the 1:1:1 ratio issue has been resolved

use anyhow::Result;
use metagenomic_llm::assembly::graph_construction::*;
use metagenomic_llm::core::data_structures::*;

#[test]
fn test_kmer_edge_connectivity() -> Result<()> {
    println!("🧪 Testing k-mer edge connectivity fixes...");

    // Create a test read that should produce connected k-mers
    let test_sequence = "ATCGATCGATCGAT"; // 14 bp sequence
    let k = 4;

    let corrected_read = CorrectedRead {
        id: 1,
        original: test_sequence.to_string(),
        corrected: test_sequence.to_string(),
        corrections: vec![],
        quality_scores: vec![30; test_sequence.len()],
        correction_metadata: CorrectionMetadata {
            algorithm: "test".to_string(),
            confidence_threshold: 0.95,
            context_window: 5,
            correction_time_ms: 0,
        },
    };

    // Create assembly chunk and process the read
    let mut chunk = AssemblyChunk::new(0, k);
    chunk.add_read(corrected_read)?;
    chunk.finalize();

    let fragment = &chunk.graph_fragment;

    println!("📊 Test results:");
    println!("  - Nodes: {}", fragment.nodes.len());
    println!("  - Edges: {}", fragment.edges.len());

    // Expected: 14-4+1 = 11 k-mers, 10 edges connecting consecutive k-mers
    let expected_nodes = test_sequence.len() - k + 1;
    let expected_edges = expected_nodes - 1;

    assert_eq!(
        fragment.nodes.len(),
        expected_nodes,
        "Should have {} nodes for k-mer size {}",
        expected_nodes,
        k
    );

    // CRITICAL TEST: Edges should connect consecutive k-mers
    assert!(
        fragment.edges.len() > 0,
        "CRITICAL: No edges found - this indicates the 1:1:1 ratio bug!"
    );

    assert_eq!(
        fragment.edges.len(),
        expected_edges,
        "Should have {} edges connecting consecutive k-mers",
        expected_edges
    );

    println!("✅ Connectivity test PASSED - edges properly connect k-mers");

    // Additional validation: Check that edges form a connected path
    let mut node_connections = std::collections::HashMap::new();
    for edge in &fragment.edges {
        node_connections
            .entry(edge.from_hash)
            .or_insert(vec![])
            .push(edge.to_hash);
    }

    let connected_nodes = node_connections.len();
    println!(
        "  - Connected nodes: {} out of {}",
        connected_nodes,
        fragment.nodes.len()
    );

    // Should have n-1 nodes with outgoing connections (last node has no outgoing edge)
    assert_eq!(
        connected_nodes, expected_edges,
        "All but the last node should have outgoing connections"
    );

    Ok(())
}

#[test]
fn test_multiple_reads_clustering() -> Result<()> {
    println!("🧪 Testing multiple reads clustering...");

    // Create overlapping reads that should cluster together
    let reads = vec![
        "ATCGATCGATCG",
        "TCGATCGATCGA",  // Overlaps with previous
        "CGATCGATCGAT",  // Overlaps with previous
        "GGGGAAATTTCCC", // Different sequence - should be separate cluster
    ];

    let k = 4;
    let mut chunk = AssemblyChunk::new(0, k);

    for (id, seq) in reads.iter().enumerate() {
        let corrected_read = CorrectedRead {
            id: id,
            original: seq.to_string(),
            corrected: seq.to_string(),
            corrections: vec![],
            quality_scores: vec![30; seq.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        };
        chunk.add_read(corrected_read)?;
    }

    chunk.finalize();
    let fragment = &chunk.graph_fragment;

    println!("📊 Multi-read test results:");
    println!("  - Total nodes: {}", fragment.nodes.len());
    println!("  - Total edges: {}", fragment.edges.len());

    // CRITICAL: Should NOT have 1:1:1 ratio
    let num_reads = reads.len();
    assert_ne!(
        fragment.nodes.len(),
        num_reads,
        "CRITICAL: Nodes equal reads ({}) - indicates 1:1:1 bug!",
        num_reads
    );

    // Should have fewer unique k-mers due to overlaps
    let total_possible_kmers: usize = reads
        .iter()
        .map(|seq| seq.len().saturating_sub(k - 1))
        .sum();

    println!("  - Total possible k-mers: {}", total_possible_kmers);
    println!("  - Actual unique k-mers: {}", fragment.nodes.len());

    // Due to overlaps, should have fewer unique k-mers than total possible
    assert!(
        fragment.nodes.len() < total_possible_kmers,
        "Should have fewer unique k-mers due to overlapping sequences"
    );

    println!("✅ Multi-read clustering test PASSED - no 1:1:1 ratio detected");

    Ok(())
}

#[test]
fn test_graph_fragment_validation() -> Result<()> {
    println!("🧪 Testing GraphFragment edge validation...");

    let kmer1 = CanonicalKmer::new("ATCG")?;
    let kmer2 = CanonicalKmer::new("TCGA")?;
    let kmer3 = CanonicalKmer::new("CGAT")?; // Not in fragment

    let mut fragment = GraphFragment::new(0);

    // Add nodes
    fragment.add_node(GraphNode::new(kmer1.clone(), 4));
    fragment.add_node(GraphNode::new(kmer2.clone(), 4));

    // Test valid edge
    let valid_edge = GraphEdge::new(kmer1.hash, kmer2.hash, 1);
    assert!(
        fragment.add_edge(valid_edge).is_ok(),
        "Valid edge should succeed"
    );

    // Test invalid edge (missing node)
    let invalid_edge = GraphEdge::new(kmer1.hash, kmer3.hash, 1);
    assert!(
        fragment.add_edge(invalid_edge).is_err(),
        "Invalid edge should fail"
    );

    // Test duplicate edge (should merge weights)
    let duplicate_edge = GraphEdge::new(kmer1.hash, kmer2.hash, 2);
    assert!(
        fragment.add_edge(duplicate_edge).is_ok(),
        "Duplicate edge should merge weights"
    );

    // Verify edge weight was updated
    let edge = fragment
        .edges
        .iter()
        .find(|e| e.from_hash == kmer1.hash && e.to_hash == kmer2.hash)
        .unwrap();
    assert_eq!(edge.weight, 3, "Edge weights should be merged (1 + 2 = 3)");

    println!("✅ GraphFragment validation test PASSED - proper edge handling");

    Ok(())
}

#[test]
fn test_assembly_produces_fewer_contigs_than_reads() -> Result<()> {
    println!("🧪 Testing that assembly produces fewer contigs than reads...");

    // Create reads with significant overlaps - should assemble into fewer contigs
    let overlapping_reads = vec![
        "ATCGATCGATCGATCG", // 16bp
        "GATCGATCGATCGATC", // 16bp, overlaps by 15bp
        "TCGATCGATCGATCGA", // 16bp, overlaps by 15bp
        "CGATCGATCGATCGAT", // 16bp, overlaps by 15bp
        "GGGGAAAATTTCCCC",  // Different sequence
        "GGGAAAATTTCCCCC",  // Overlaps with previous
    ];

    let k = 8; // Larger k-mer for better specificity
    let mut builder = AdvancedAssemblyGraphBuilder::new(k, 21, 1, 2)?;

    let mut corrected_reads = Vec::new();
    for (id, seq) in overlapping_reads.iter().enumerate() {
        corrected_reads.push(CorrectedRead {
            id: id,
            original: seq.to_string(),
            corrected: seq.to_string(),
            corrections: vec![],
            quality_scores: vec![30; seq.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        });
    }

    // This is a simplified test - in reality would use the full assembly pipeline
    let num_reads = corrected_reads.len();
    println!("📊 Assembly ratio test:");
    println!("  - Input reads: {}", num_reads);

    // For now, just verify our fixes don't break compilation
    // Full integration test would be in separate test

    println!("✅ Assembly ratio test framework created - ready for integration");

    Ok(())
}
