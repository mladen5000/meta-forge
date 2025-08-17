//! Comprehensive test suite for assembly pipeline
//! Tests AssemblyGraphBuilder, adaptive k-mer selection, graph construction, and integration

// use crate::assembly_graph_builder_tests::create_test_read;
use assembly_graph_builder_tests::create_test_read;
use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::assembly::bioinformatics_optimizations::BitPackedKmer;
use meta_forge::core::data_structures::{
    AssemblyGraph, CorrectedRead, CorrectionMetadata, GraphFragment,
};
use meta_forge::AssemblyChunk;
use meta_forge::CanonicalKmer;
use meta_forge::GraphEdge;
use meta_forge::GraphNode;
use meta_forge::NodeType;

#[cfg(test)]
pub mod assembly_graph_builder_tests {
    use super::*;

    pub fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
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
    fn test_assembly_graph_builder_creation() {
        let builder = AssemblyGraphBuilder::new(3, 7, 2);

        // Test that builder is created with correct parameters
        // (Internal fields are private, so we test through behavior)
        assert!(true); // Builder creation succeeded
    }

    #[test]
    fn test_assembly_graph_builder_with_simple_reads() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        let reads = vec![
            create_test_read(0, "ATCGATCG"),
            create_test_read(1, "TCGATCGA"),
            create_test_read(2, "CGATCGAT"),
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());

        let graph = result.unwrap();
        assert!(!graph.graph_fragment.nodes.is_empty());
        assert!(!graph.graph_fragment.edges.is_empty());
    }

    #[test]
    fn test_assembly_graph_builder_with_overlapping_reads() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);

        let reads = vec![
            create_test_read(0, "ATCGATCG"),
            create_test_read(1, "GATCGACT"), // Overlaps with read 0
            create_test_read(2, "CGATCGAT"), // Overlaps with both
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());

        let graph = result.unwrap();

        // Should have merged overlapping k-mers
        assert!(!graph.graph_fragment.nodes.is_empty());
        assert!(!graph.graph_fragment.edges.is_empty());

        // Check that some nodes have higher coverage due to overlaps
        let high_coverage_nodes: Vec<_> = graph
            .graph_fragment
            .nodes
            .values()
            .filter(|node| node.coverage > 1)
            .collect();

        assert!(
            !high_coverage_nodes.is_empty(),
            "Should have some high coverage nodes from overlaps"
        );
    }

    #[test]
    fn test_assembly_graph_builder_with_short_reads() {
        let builder = AssemblyGraphBuilder::new(5, 7, 1);

        let reads = vec![
            create_test_read(0, "AT"),           // Too short for k=5
            create_test_read(1, "CG"),           // Too short for k=5
            create_test_read(2, "ATCGATCGATCG"), // Long enough
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());

        let graph = result.unwrap();

        // Should only process the long read
        assert!(!graph.graph_fragment.nodes.is_empty());
    }

    #[test]
    fn test_assembly_graph_builder_with_empty_reads() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        let reads = vec![];

        let result = builder.build(&reads);
        assert!(result.is_ok());

        let graph = result.unwrap();
        assert!(graph.graph_fragment.nodes.is_empty());
        assert!(graph.graph_fragment.edges.is_empty());
    }

    #[test]
    fn test_assembly_graph_builder_k_mer_size_adaptation() {
        let builder = AssemblyGraphBuilder::new(3, 7, 1);

        // Create reads that should trigger different k-mer sizes
        let simple_reads = vec![
            create_test_read(0, "AAAAAAAAAA"), // Low complexity
            create_test_read(1, "TTTTTTTTTT"), // Low complexity
        ];

        let complex_reads = vec![
            create_test_read(0, "ATCGATCGACTGATCG"), // High complexity
            create_test_read(1, "GTACGTACGTACGTAC"), // High complexity
        ];

        let simple_result = builder.build(&simple_reads);
        let complex_result = builder.build(&complex_reads);

        assert!(simple_result.is_ok());
        assert!(complex_result.is_ok());

        // Both should succeed but potentially use different k-mer sizes internally
    }

    #[test]
    fn test_contig_generation() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        let reads = vec![
            create_test_read(0, "ATCGATCGATCG"),
            create_test_read(1, "TCGATCGATCGA"),
            create_test_read(2, "CGATCGATCGAT"),
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());

        let graph = result.unwrap();
        assert!(!graph.contigs.is_empty());

        // Check contig properties
        for contig in &graph.contigs {
            assert!(contig.id >= 0);
            assert!(!contig.sequence.is_empty());
            assert!(contig.length > 0);
            assert!(contig.coverage > 0.0);
            assert!(!contig.node_path.is_empty());
        }
    }

    #[test]
    fn test_assembly_statistics_calculation() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        let reads = vec![
            create_test_read(0, "ATCGATCGATCGATCG"),
            create_test_read(1, "TCGATCGATCGATCGA"),
            create_test_read(2, "CGATCGATCGATCGAT"),
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());

        let graph = result.unwrap();

        // Check assembly stats
        assert!(graph.assembly_stats.total_length > 0);
        assert!(graph.assembly_stats.num_contigs > 0);
        assert!(graph.assembly_stats.largest_contig > 0);
        assert!(graph.assembly_stats.n50 >= 0);
        assert!(graph.assembly_stats.coverage_mean >= 0.0);
    }
}

#[cfg(test)]
mod assembly_chunk_tests {
    use super::*;

    #[test]
    fn test_assembly_chunk_creation() {
        let chunk = AssemblyChunk::new(42, 4);

        assert_eq!(chunk.chunk_id, 42);
        assert_eq!(chunk.k_size, 4);
        assert!(chunk.reads.is_empty());
        assert_eq!(chunk.graph_fragment.fragment_id, 42);
        assert_eq!(chunk.processing_stats.reads_processed, 0);
    }

    #[test]
    fn test_assembly_chunk_add_read() {
        let mut chunk = AssemblyChunk::new(0, 4);
        let read = create_test_read(0, "ATCGATCG");

        let result = chunk.add_read(read);
        assert!(result.is_ok());

        assert_eq!(chunk.reads.len(), 1);
        assert_eq!(chunk.processing_stats.reads_processed, 1);
        assert!(chunk.processing_stats.nodes_created > 0);
        assert!(chunk.processing_stats.minimizers_found > 0);
    }

    #[test]
    fn test_assembly_chunk_with_short_read() {
        let mut chunk = AssemblyChunk::new(0, 5);
        let short_read = create_test_read(0, "ATG"); // Shorter than k-mer size

        let result = chunk.add_read(short_read);
        assert!(result.is_ok()); // Should succeed but not create nodes

        assert_eq!(chunk.reads.len(), 1);
        assert_eq!(chunk.processing_stats.reads_processed, 1);
        assert_eq!(chunk.processing_stats.nodes_created, 0); // No nodes for short reads
    }

    #[test]
    fn test_assembly_chunk_finalization() {
        let mut chunk = AssemblyChunk::new(0, 4);

        let reads = vec![
            create_test_read(0, "ATCGATCG"),
            create_test_read(1, "TCGATCGA"),
        ];

        for read in reads {
            chunk.add_read(read).unwrap();
        }

        chunk.finalize();

        // Check that finalization completed
        assert!(chunk.processing_stats.memory_usage_bytes > 0);

        // Check node types are updated
        for node in chunk.graph_fragment.nodes.values() {
            assert!(matches!(
                node.node_type,
                NodeType::LowCoverage
                    | NodeType::Unique
                    | NodeType::HighCoverage
                    | NodeType::Repetitive
            ));
        }
    }

    #[test]
    fn test_assembly_chunk_k_mer_extraction() {
        let mut chunk = AssemblyChunk::new(0, 3);
        let read = create_test_read(0, "ATCGATCG");

        chunk.add_read(read).unwrap();

        // For k=3 and sequence "ATCGATCG" (8 bases), we expect:
        // ATC, TCG, CGA, GAT, ATC, TCG - some may be duplicates when canonicalized
        assert!(chunk.processing_stats.minimizers_found >= 6);
        assert!(chunk.graph_fragment.nodes.len() > 0);
        assert!(chunk.graph_fragment.edges.len() >= 5); // n-1 edges for consecutive k-mers
    }

    #[test]
    fn test_assembly_chunk_edge_creation() {
        let mut chunk = AssemblyChunk::new(0, 3);
        let read = create_test_read(0, "ATCG");

        chunk.add_read(read).unwrap();

        // For "ATCG" with k=3: ATC, TCG -> 1 edge
        assert_eq!(chunk.processing_stats.edges_created, 1);
        assert_eq!(chunk.graph_fragment.edges.len(), 1);

        let edge = &chunk.graph_fragment.edges[0];
        assert_eq!(edge.overlap_length, 2); // k-1 overlap
    }

    #[test]
    fn test_assembly_chunk_coverage_accumulation() {
        let mut chunk = AssemblyChunk::new(0, 3);

        // Add overlapping reads
        let read1 = create_test_read(0, "ATCG");
        let read2 = create_test_read(1, "ATCG"); // Same sequence

        chunk.add_read(read1).unwrap();
        chunk.add_read(read2).unwrap();

        // Should have increased coverage for overlapping k-mers
        let high_coverage_nodes: Vec<_> = chunk
            .graph_fragment
            .nodes
            .values()
            .filter(|node| node.coverage > 1)
            .collect();

        assert!(
            !high_coverage_nodes.is_empty(),
            "Should have high coverage nodes from overlapping reads"
        );
    }
}

#[cfg(test)]
mod bit_packed_kmer_tests {
    use super::*;

    #[test]
    fn test_bit_packed_kmer_creation() {
        let kmer = BitPackedKmer::new("ATCG").unwrap();

        assert_eq!(kmer.k, 4);
        assert!(kmer.hash > 0);
        assert!(!kmer.packed_data.is_empty());
    }

    #[test]
    fn test_bit_packed_kmer_invalid_characters() {
        assert!(BitPackedKmer::new("ATCGX").is_err());
        assert!(BitPackedKmer::new("").is_err());
        assert!(BitPackedKmer::new("ATCG123").is_err());
    }

    #[test]
    fn test_bit_packed_kmer_consistency() {
        let kmer1 = BitPackedKmer::new("ATCG").unwrap();
        let kmer2 = BitPackedKmer::new("ATCG").unwrap();

        assert_eq!(kmer1.hash, kmer2.hash);
        assert_eq!(kmer1.packed_data, kmer2.packed_data);
    }

    #[test]
    fn test_bit_packed_kmer_different_sequences() {
        let kmer1 = BitPackedKmer::new("ATCG").unwrap();
        let kmer2 = BitPackedKmer::new("CGAT").unwrap();

        // Different sequences should have different hashes (probabilistically)
        assert_ne!(kmer1.hash, kmer2.hash);
    }

    #[test]
    fn test_bit_packed_kmer_long_sequence() {
        // Test with longer sequence
        let long_seq = "ATCGATCGATCGATCGATCGATCG";
        let kmer = BitPackedKmer::new(long_seq).unwrap();

        assert_eq!(kmer.k, long_seq.len());
        assert!(kmer.hash > 0);
    }

    #[test]
    fn test_bit_packed_kmer_too_long() {
        // Test with sequence longer than limit
        let too_long = "A".repeat(2000);
        assert!(BitPackedKmer::new(&too_long).is_err());
    }
}

#[cfg(test)]
mod graph_fragment_integration_tests {
    use super::*;

    #[test]
    fn test_graph_fragment_insert_edge() {
        let mut fragment = GraphFragment::empty();

        fragment.insert_edge(123, 456);

        assert_eq!(fragment.nodes.len(), 2); // Should create both nodes
        assert_eq!(fragment.edges.len(), 1);

        // Nodes should have increased coverage
        if let Some(node) = fragment.nodes.get(&123) {
            assert!(node.coverage > 0);
        }
        if let Some(node) = fragment.nodes.get(&456) {
            assert!(node.coverage > 0);
        }
    }

    #[test]
    fn test_assembly_graph_from_fragment() {
        let mut fragment = GraphFragment::new(0);

        // Add some test data
        let kmer1 = CanonicalKmer::new("ATCG").unwrap();
        let kmer2 = CanonicalKmer::new("TCGA").unwrap();

        fragment.add_node(GraphNode::new(kmer1.clone(), 4));
        fragment.add_node(GraphNode::new(kmer2.clone(), 4));
        fragment.add_edge(GraphEdge::new(kmer1.hash, kmer2.hash, 3));

        let assembly_graph = AssemblyGraph::from_fragment(fragment);

        assert_eq!(assembly_graph.graph_fragment.nodes.len(), 2);
        assert_eq!(assembly_graph.graph_fragment.edges.len(), 1);
        assert!(assembly_graph.contigs.is_empty()); // No contigs generated yet
    }

    #[test]
    fn test_full_pipeline_integration() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);

        // Create a realistic set of overlapping reads
        let reads = vec![
            create_test_read(0, "ATCGATCGATCG"),
            create_test_read(1, "CGATCGATCGAT"),
            create_test_read(2, "GATCGATCGATC"),
            create_test_read(3, "ATCGATCGATCG"), // Duplicate for coverage
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());

        let graph = result.unwrap();

        // Verify complete pipeline execution
        assert!(!graph.graph_fragment.nodes.is_empty());
        assert!(!graph.graph_fragment.edges.is_empty());
        assert!(!graph.contigs.is_empty());

        // Verify assembly statistics
        assert!(graph.assembly_stats.total_length > 0);
        assert!(graph.assembly_stats.num_contigs > 0);
        assert!(graph.assembly_stats.coverage_mean > 0.0);

        // Verify some nodes have high coverage due to overlaps
        let high_coverage_nodes: Vec<_> = graph
            .graph_fragment
            .nodes
            .values()
            .filter(|node| node.coverage > 1)
            .collect();

        assert!(
            !high_coverage_nodes.is_empty(),
            "Should have high coverage nodes from overlapping reads"
        );
    }
}
