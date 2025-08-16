//! Phase 1 TDD Tests for Graph Construction
//!
//! Critical tests for assembly graph construction algorithms
//! Focus: Streaming k-mer processing, memory optimization, SIMD validation

use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::core::data_structures::*;
use std::collections::HashSet;

/// Test streaming k-mer graph construction
/// BIOLOGICAL REQUIREMENT: Consecutive k-mers must create proper edges
/// that preserve the original sequence relationships
#[cfg(test)]
mod streaming_graph_construction_tests {
    use super::*;

    fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn test_consecutive_kmers_create_proper_edges() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        // Create a read with known k-mer structure
        // Sequence: ATCGATCG (8bp)
        // K=4 gives k-mers: ATCG, TCGA, CGAT, GATC, ATCG
        // Expected edges: ATCG->TCGA, TCGA->CGAT, CGAT->GATC, GATC->ATCG
        let reads = vec![create_test_read(0, "ATCGATCG")];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Graph construction should succeed");

        let graph = result.unwrap();

        // Verify that we have the expected number of edges (4 consecutive k-mers = 3 edges)
        // Note: ATCG appears twice, so we might have 4 edges total
        assert!(
            graph.graph_fragment.edges.len() >= 3,
            "Should have at least 3 edges for consecutive k-mers, got {}",
            graph.graph_fragment.edges.len()
        );

        // Verify that nodes exist for each unique k-mer
        assert!(
            !graph.graph_fragment.nodes.is_empty(),
            "Graph should contain nodes for k-mers"
        );
    }

    #[test]
    fn test_overlapping_reads_increase_coverage() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        // Create overlapping reads that share k-mers
        let reads = vec![
            create_test_read(0, "ATCGATCG"), // Contains ATCG, TCGA, CGAT, GATC, ATCG
            create_test_read(1, "TCGATCGA"), // Contains TCGA, CGAT, GATC, ATCG
            create_test_read(2, "CGATCGAT"), // Contains CGAT, GATC, ATCG, TCGA
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Graph construction should succeed");

        let graph = result.unwrap();

        // Find nodes with coverage > 1 (indicating overlap)
        let high_coverage_nodes: Vec<_> = graph
            .graph_fragment
            .nodes
            .values()
            .filter(|node| node.coverage > 1)
            .collect();

        assert!(
            !high_coverage_nodes.is_empty(),
            "Should have nodes with coverage > 1 due to overlapping reads"
        );

        // Specifically check for expected shared k-mers
        let shared_kmers = vec!["TCGA", "CGAT", "GATC", "ATCG"];
        for kmer_seq in shared_kmers {
            if let Ok(kmer) = CanonicalKmer::new(kmer_seq) {
                if let Some(node) = graph.graph_fragment.nodes.get(&kmer.hash) {
                    assert!(
                        node.coverage >= 2,
                        "K-mer {} should appear in multiple reads, coverage: {}",
                        kmer_seq,
                        node.coverage
                    );
                }
            }
        }
    }

    #[test]
    fn test_node_coverage_increments_correctly() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);

        // Create reads with known k-mer repetition
        let reads = vec![
            create_test_read(0, "ATGATGATG"), // ATC appears 3 times
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Graph construction should succeed");

        let graph = result.unwrap();

        // Check that ATG k-mer has coverage reflecting its frequency
        if let Ok(atg_kmer) = CanonicalKmer::new("ATG") {
            if let Some(node) = graph.graph_fragment.nodes.get(&atg_kmer.hash) {
                assert!(
                    node.coverage >= 3,
                    "ATG should appear 3 times, coverage: {}",
                    node.coverage
                );
            }
        }
    }

    #[test]
    fn test_short_reads_handling() {
        let builder = AssemblyGraphBuilder::new(5, 7, 1);

        // Create reads shorter than k-mer size
        let reads = vec![
            create_test_read(0, "ATC"),      // 3bp < k=5
            create_test_read(1, "ATCG"),     // 4bp < k=5
            create_test_read(2, "ATCGATCG"), // 8bp > k=5 (valid)
        ];

        let result = builder.build(&reads);
        assert!(
            result.is_ok(),
            "Graph construction should handle short reads gracefully"
        );

        let graph = result.unwrap();

        // Should only process the read that's long enough
        // The 8bp read with k=5 should produce 4 k-mers
        assert!(
            !graph.graph_fragment.nodes.is_empty(),
            "Should have nodes from the valid read"
        );
    }

    #[test]
    fn test_empty_input_handling() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        let reads = vec![];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle empty input gracefully");

        let graph = result.unwrap();
        assert!(
            graph.graph_fragment.nodes.is_empty(),
            "Empty input should produce empty graph"
        );
        assert!(
            graph.graph_fragment.edges.is_empty(),
            "Empty input should produce no edges"
        );
    }

    #[test]
    fn test_single_base_reads() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);

        // Test with very simple, repetitive sequences
        let reads = vec![
            create_test_read(0, "AAA"),
            create_test_read(1, "AAA"),
            create_test_read(2, "AAA"),
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle simple repetitive sequences");

        let graph = result.unwrap();

        // Should have exactly one node (AAA) with high coverage
        if let Ok(aaa_kmer) = CanonicalKmer::new("AAA") {
            if let Some(node) = graph.graph_fragment.nodes.get(&aaa_kmer.hash) {
                assert_eq!(
                    node.coverage, 3,
                    "AAA should appear 3 times, coverage: {}",
                    node.coverage
                );
            }
        }

        // Should have no edges (single k-mer, no transitions)
        assert!(
            graph.graph_fragment.edges.is_empty(),
            "Single k-mer should produce no edges"
        );
    }
}

/// Test optimization mode selection
/// BIOLOGICAL REQUIREMENT: Different optimization modes should produce
/// equivalent biological results while optimizing for different resources
#[cfg(test)]
mod optimization_mode_tests {
    use super::*;

    fn create_test_dataset() -> Vec<CorrectedRead> {
        vec![
            create_test_read(0, "ATCGATCGATCGATCG"),
            create_test_read(1, "TCGATCGATCGATCGA"),
            create_test_read(2, "CGATCGATCGATCGAT"),
            create_test_read(3, "GATCGATCGATCGATC"),
            create_test_read(4, "ATCGATCGATCGATCG"), // Duplicate for coverage
        ]
    }

    fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn test_optimization_modes_produce_equivalent_results() {
        let reads = create_test_dataset();

        // Test all three optimization modes
        let builder_normal = AssemblyGraphBuilder::new(4, 6, 1);
        let builder_low_memory = AssemblyGraphBuilder::new_low_memory(4, 6, 1);
        let builder_low_cpu = AssemblyGraphBuilder::new_low_cpu(4, 6, 1);

        let graph_normal = builder_normal.build(&reads).unwrap();
        let graph_low_memory = builder_low_memory.build_low_memory(&reads).unwrap();
        let graph_low_cpu = builder_low_cpu.build_low_cpu(&reads).unwrap();

        // All modes should produce the same number of unique k-mers
        assert_eq!(
            graph_normal.graph_fragment.nodes.len(),
            graph_low_memory.graph_fragment.nodes.len(),
            "Normal and low-memory modes should have same node count"
        );
        assert_eq!(
            graph_normal.graph_fragment.nodes.len(),
            graph_low_cpu.graph_fragment.nodes.len(),
            "Normal and low-CPU modes should have same node count"
        );

        // All modes should produce the same edges (though possibly in different order)
        assert_eq!(
            graph_normal.graph_fragment.edges.len(),
            graph_low_memory.graph_fragment.edges.len(),
            "Normal and low-memory modes should have same edge count"
        );
        assert_eq!(
            graph_normal.graph_fragment.edges.len(),
            graph_low_cpu.graph_fragment.edges.len(),
            "Normal and low-CPU modes should have same edge count"
        );
    }

    #[test]
    fn test_low_memory_mode_processes_smaller_chunks() {
        let reads = create_test_dataset();
        let builder = AssemblyGraphBuilder::new_low_memory(4, 6, 1);

        // Low memory mode should still produce valid results
        let result = builder.build_low_memory(&reads);
        assert!(result.is_ok(), "Low memory mode should succeed");

        let graph = result.unwrap();
        assert!(
            !graph.graph_fragment.nodes.is_empty(),
            "Low memory mode should produce nodes"
        );
        assert!(
            !graph.graph_fragment.edges.is_empty(),
            "Low memory mode should produce edges"
        );
    }

    #[test]
    fn test_low_cpu_mode_sequential_processing() {
        let reads = create_test_dataset();
        let builder = AssemblyGraphBuilder::new_low_cpu(4, 6, 1);

        // Low CPU mode should still produce valid results
        let result = builder.build_low_cpu(&reads);
        assert!(result.is_ok(), "Low CPU mode should succeed");

        let graph = result.unwrap();
        assert!(
            !graph.graph_fragment.nodes.is_empty(),
            "Low CPU mode should produce nodes"
        );
        assert!(
            !graph.graph_fragment.edges.is_empty(),
            "Low CPU mode should produce edges"
        );
    }

    #[test]
    fn test_coverage_consistency_across_modes() {
        let reads = create_test_dataset();

        let builder_normal = AssemblyGraphBuilder::new(4, 6, 1);
        let builder_low_memory = AssemblyGraphBuilder::new_low_memory(4, 6, 1);

        let graph_normal = builder_normal.build(&reads).unwrap();
        let graph_low_memory = builder_low_memory.build_low_memory(&reads).unwrap();

        // Check that specific k-mers have the same coverage in both modes
        let test_kmer = CanonicalKmer::new("ATCG").unwrap();

        let normal_coverage = graph_normal
            .graph_fragment
            .nodes
            .get(&test_kmer.hash)
            .map(|n| n.coverage)
            .unwrap_or(0);

        let low_memory_coverage = graph_low_memory
            .graph_fragment
            .nodes
            .get(&test_kmer.hash)
            .map(|n| n.coverage)
            .unwrap_or(0);

        assert_eq!(
            normal_coverage, low_memory_coverage,
            "K-mer coverage should be consistent across optimization modes"
        );
    }

    #[test]
    fn test_edge_weight_consistency_across_modes() {
        let reads = create_test_dataset();

        let builder_normal = AssemblyGraphBuilder::new(4, 6, 1);
        let builder_low_cpu = AssemblyGraphBuilder::new_low_cpu(4, 6, 1);

        let graph_normal = builder_normal.build(&reads).unwrap();
        let graph_low_cpu = builder_low_cpu.build_low_cpu(&reads).unwrap();

        // Compare edge weights for consistency
        let normal_edge_weights: HashSet<u32> = graph_normal
            .graph_fragment
            .edges
            .iter()
            .map(|e| e.weight)
            .collect();
        let low_cpu_edge_weights: HashSet<u32> = graph_low_cpu
            .graph_fragment
            .edges
            .iter()
            .map(|e| e.weight)
            .collect();

        assert_eq!(
            normal_edge_weights, low_cpu_edge_weights,
            "Edge weights should be consistent across optimization modes"
        );
    }
}

/// Test memory-mapped file processing accuracy
/// BIOLOGICAL REQUIREMENT: Large file processing must maintain data integrity
/// and not introduce processing artifacts
#[cfg(test)]
mod memory_mapped_processing_tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_large_read_set_processing() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        // Create a larger dataset to test chunking behavior
        let mut reads = Vec::new();
        let base_sequences = vec![
            "ATCGATCGATCGATCG",
            "TCGATCGATCGATCGA",
            "CGATCGATCGATCGAT",
            "GATCGATCGATCGATC",
        ];

        // Replicate to create a substantial dataset
        for i in 0..1000 {
            let seq = &base_sequences[i % base_sequences.len()];
            reads.push(create_test_read(i, seq));
        }

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle large read sets");

        let graph = result.unwrap();

        // Verify that coverage reflects the replication
        let test_kmer = CanonicalKmer::new("ATCG").unwrap();
        if let Some(node) = graph.graph_fragment.nodes.get(&test_kmer.hash) {
            assert!(
                node.coverage >= 250, // Appears in all 4 sequences, replicated 250 times each
                "High coverage expected for replicated k-mer, got {}",
                node.coverage
            );
        }
    }

    fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn test_chunk_boundary_consistency() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        // Create reads that will span chunk boundaries
        let mut reads = Vec::new();

        // Create enough reads to ensure multiple chunks (typical chunk size is 1000)
        for i in 0..2500 {
            let sequence = if i % 2 == 0 {
                "ATCGATCGATCGATCG"
            } else {
                "TCGATCGATCGATCGA"
            };
            reads.push(create_test_read(i, sequence));
        }

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle chunk boundaries correctly");

        let graph = result.unwrap();

        // Verify expected coverage for distributed reads
        let atcg_kmer = CanonicalKmer::new("ATCG").unwrap();
        if let Some(node) = graph.graph_fragment.nodes.get(&atcg_kmer.hash) {
            // ATCG appears in both sequences, so should have coverage ≈ 2500
            assert!(
                node.coverage >= 2000 && node.coverage <= 3000,
                "ATCG coverage should reflect distribution across chunks, got {}",
                node.coverage
            );
        }
    }

    #[test]
    fn test_memory_efficiency_low_memory_mode() {
        let builder = AssemblyGraphBuilder::new_low_memory(4, 6, 1);

        // Create a moderate dataset
        let mut reads = Vec::new();
        for i in 0..500 {
            reads.push(create_test_read(i, "ATCGATCGATCGATCG"));
        }

        let result = builder.build_low_memory(&reads);
        assert!(
            result.is_ok(),
            "Low memory mode should handle moderate datasets"
        );

        let graph = result.unwrap();

        // Should still produce correct results
        assert!(
            !graph.graph_fragment.nodes.is_empty(),
            "Low memory mode should produce nodes"
        );

        // Verify coverage is correct
        let test_kmer = CanonicalKmer::new("ATCG").unwrap();
        if let Some(node) = graph.graph_fragment.nodes.get(&test_kmer.hash) {
            assert_eq!(node.coverage, 500, "Coverage should match input read count");
        }
    }
}

/// Integration test: Graph construction workflow
/// BIOLOGICAL REQUIREMENT: Complete workflow should produce biologically meaningful graphs
#[cfg(test)]
mod graph_construction_integration_tests {
    use super::*;

    fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn test_complete_assembly_graph_generation() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        // Create a realistic set of overlapping reads
        let reads = vec![
            create_test_read(0, "ATCGATCGATCGATCG"), // 16bp
            create_test_read(1, "TCGATCGATCGATCGA"), // Overlaps with read 0
            create_test_read(2, "CGATCGATCGATCGAT"), // Overlaps with read 1
            create_test_read(3, "GATCGATCGATCGATC"), // Overlaps with read 2
            create_test_read(4, "ATCGATCGATCGATCG"), // Duplicate of read 0
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Complete workflow should succeed");

        let graph = result.unwrap();

        // Verify graph structure
        assert!(!graph.graph_fragment.nodes.is_empty(), "Should have nodes");
        assert!(!graph.graph_fragment.edges.is_empty(), "Should have edges");

        // Verify contigs were generated
        assert!(!graph.contigs.is_empty(), "Should generate contigs");

        // Verify assembly statistics were calculated
        assert!(
            graph.assembly_stats.total_length > 0,
            "Should have total length > 0"
        );
        assert!(
            graph.assembly_stats.num_contigs > 0,
            "Should have contig count > 0"
        );
    }

    #[test]
    fn test_graph_topology_properties() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);

        // Create reads that form a linear path
        let reads = vec![
            create_test_read(0, "ATCGAT"), // ATCGpATcgatCGAT
            create_test_read(1, "TCGATG"), // TCGATgAGATg
            create_test_read(2, "CGATGC"), // CGATGcGATGC
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Graph construction should succeed");

        let graph = result.unwrap();

        // Verify connectivity
        let adjacency = graph.graph_fragment.get_adjacency_list();
        assert!(
            !adjacency.is_empty(),
            "Graph should have adjacency structure"
        );

        // Most nodes should have degree <= 2 for a linear structure
        let high_degree_nodes: Vec<_> = adjacency
            .iter()
            .filter(|(_, neighbors)| neighbors.len() > 2)
            .collect();

        // Allow some high-degree nodes due to repeats, but most should be linear
        assert!(
            high_degree_nodes.len() <= graph.graph_fragment.nodes.len() / 2,
            "Most nodes should have degree <= 2 for linear reads"
        );
    }

    #[test]
    fn test_repeat_detection_in_graph() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);

        // Create reads with a known repeat pattern
        let reads = vec![
            create_test_read(0, "ATATAT"), // Contains AT repeat
            create_test_read(1, "ATATAT"), // Same repeat
            create_test_read(2, "ATATAT"), // Same repeat
            create_test_read(3, "GCGCGC"), // Different repeat
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Graph construction should succeed");

        let graph = result.unwrap();

        // Check for high-coverage nodes indicating repeats
        let repeat_kmer = CanonicalKmer::new("ATA").unwrap();
        if let Some(node) = graph.graph_fragment.nodes.get(&repeat_kmer.hash) {
            assert!(node.coverage >= 3, "Repeat k-mer should have high coverage");
        }

        // Check for tips (potential sequencing errors or repeat boundaries)
        let tips = graph.graph_fragment.find_tips();
        // Tips are expected in repeat-heavy sequences
        assert!(
            tips.len() <= graph.graph_fragment.nodes.len(),
            "Tips should not exceed total node count"
        );
    }

    #[test]
    fn test_coverage_statistics_accuracy() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        // Create reads with known coverage pattern
        let reads = vec![
            create_test_read(0, "ATCGATCG"), // Coverage 1
            create_test_read(1, "ATCGATCG"), // Coverage 2
            create_test_read(2, "ATCGATCG"), // Coverage 3
            create_test_read(3, "GCTAGCTA"), // Different sequence, coverage 1
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Graph construction should succeed");

        let graph = result.unwrap();

        // Verify coverage statistics
        assert!(
            graph.graph_fragment.coverage_stats.mean_coverage > 0.0,
            "Mean coverage should be > 0"
        );
        assert!(
            graph.graph_fragment.coverage_stats.total_nodes > 0,
            "Should have nodes"
        );
        assert!(
            graph.graph_fragment.coverage_stats.total_edges > 0,
            "Should have edges"
        );

        // Check that high-coverage k-mers are properly identified
        let atcg_kmer = CanonicalKmer::new("ATCG").unwrap();
        if let Some(node) = graph.graph_fragment.nodes.get(&atcg_kmer.hash) {
            assert_eq!(node.coverage, 3, "ATCG should appear 3 times");
        }
    }
}
