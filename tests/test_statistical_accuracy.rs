//! Statistical accuracy tests for abundance estimation algorithms
//! Tests coverage calculation, complexity scoring, statistical metrics, and quality assessment

use abundance_estimation_tests::create_test_read_with_quality;

use meta_forge::assembly::adaptive_k::AssemblyGraph;
use meta_forge::calculate_gc_content;
use meta_forge::calculate_sequence_complexity;
use meta_forge::core::data_structures::{
    CorrectedRead, CorrectionMetadata, GraphEdge, GraphFragment, GraphNode,
};
use meta_forge::AssemblyChunk;
use meta_forge::CanonicalKmer;
use meta_forge::Contig;
use meta_forge::ContigType;
use meta_forge::Strand;
use std::collections::HashMap;

#[cfg(test)]
pub mod abundance_estimation_tests {
    use super::*;

    pub fn create_test_read_with_quality(
        id: usize,
        sequence: &str,
        qualities: Vec<u8>,
    ) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: qualities,
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
            kmer_hash_cache: Vec::new(),
        }
    }

    #[test]
    fn test_coverage_calculation_accuracy() {
        let mut fragment = GraphFragment::new(0);

        // Create nodes with known coverage values
        let test_cases = vec![
            ("ATCG", 1),
            ("TCGA", 5),
            ("CGAT", 10),
            ("GATC", 15),
            ("ATGC", 25),
        ];

        for (seq, coverage) in &test_cases {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let mut node = GraphNode::new(kmer, 4);
            node.coverage = *coverage;
            fragment.add_node(node);
        }

        // Test statistics calculation
        let expected_mean =
            test_cases.iter().map(|(_, c)| *c as f64).sum::<f64>() / test_cases.len() as f64;
        let mut sorted_coverages: Vec<u32> = test_cases.iter().map(|(_, c)| *c).collect();
        sorted_coverages.sort_unstable();
        let expected_median = sorted_coverages[sorted_coverages.len() / 2];

        assert_eq!(fragment.coverage_stats.mean_coverage, expected_mean);
        assert_eq!(fragment.coverage_stats.median_coverage, expected_median);
        assert_eq!(fragment.coverage_stats.total_nodes, test_cases.len());
    }

    #[test]
    fn test_coverage_distribution_accuracy() {
        let mut fragment = GraphFragment::new(0);

        // Create nodes with specific coverage distribution
        let coverage_counts = HashMap::from([
            (1, 3),  // 3 nodes with coverage 1
            (5, 2),  // 2 nodes with coverage 5
            (10, 1), // 1 node with coverage 10
        ]);

        let mut node_id = 0;
        for (&coverage, &count) in &coverage_counts {
            for _ in 0..count {
                let seq = format!("A{node_id:03}");
                let kmer = CanonicalKmer::new(&seq).unwrap();
                let mut node = GraphNode::new(kmer, 4);
                node.coverage = coverage;
                fragment.add_node(node);
                node_id += 1;
            }
        }

        // Verify distribution matches
        for (&expected_coverage, &expected_count) in &coverage_counts {
            let actual_count = fragment
                .coverage_stats
                .coverage_distribution
                .get(&expected_coverage)
                .copied()
                .unwrap_or(0);
            assert_eq!(
                actual_count, expected_count,
                "Coverage {expected_coverage} should have {expected_count} nodes"
            );
        }

        // Verify total counts match
        let total_nodes: usize = coverage_counts.values().sum();
        assert_eq!(fragment.coverage_stats.total_nodes, total_nodes);
    }

    #[test]
    fn test_path_coverage_estimation() {
        let mut fragment = GraphFragment::new(0);

        // Create path with varying coverage
        let path_data = vec![("ATCG", 10), ("TCGA", 20), ("CGAT", 30), ("GATC", 40)];

        let mut path_hashes = Vec::new();
        for (seq, coverage) in &path_data {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let mut node = GraphNode::new(kmer.clone(), 4);
            node.coverage = *coverage;
            path_hashes.push(kmer.hash);
            fragment.add_node(node);
        }

        let avg_coverage = fragment.calculate_path_coverage_from_hashes(&path_hashes);
        let expected_avg =
            path_data.iter().map(|(_, c)| *c as f64).sum::<f64>() / path_data.len() as f64;

        assert_eq!(avg_coverage, expected_avg);
    }

    #[test]
    fn test_coverage_with_overlapping_reads() {
        let mut chunk = AssemblyChunk::new(0, 3);

        // Add overlapping reads to increase coverage
        let base_sequence = "ATCGATCG";
        let overlapping_reads = vec![
            create_test_read_with_quality(0, base_sequence, vec![30; base_sequence.len()]),
            create_test_read_with_quality(
                1,
                &base_sequence[1..],
                vec![35; base_sequence.len() - 1],
            ), // Offset by 1
            create_test_read_with_quality(
                2,
                &base_sequence[2..],
                vec![25; base_sequence.len() - 2],
            ), // Offset by 2
        ];

        for read in overlapping_reads {
            chunk.add_read(read).unwrap();
        }

        chunk.finalize();

        // Some k-mers should have higher coverage due to overlaps
        let coverages: Vec<u32> = chunk
            .graph_fragment
            .nodes
            .values()
            .map(|n| n.coverage)
            .collect();
        let max_coverage = coverages.iter().max().unwrap();
        let min_coverage = coverages.iter().min().unwrap();

        // With overlapping reads, we should see varying coverage
        assert!(
            max_coverage > min_coverage,
            "Overlapping reads should create coverage variation"
        );
        assert!(
            *max_coverage > 1,
            "Some k-mers should have coverage > 1 from overlaps"
        );
    }

    #[test]
    fn test_abundance_estimation_with_quality_scores() {
        let mut chunk = AssemblyChunk::new(0, 4);

        // Reads with different quality profiles
        let reads = vec![
            create_test_read_with_quality(0, "ATCGATCG", vec![40, 40, 40, 40, 40, 40, 40, 40]), // High quality
            create_test_read_with_quality(1, "ATCGATCG", vec![20, 20, 20, 20, 20, 20, 20, 20]), // Low quality
            create_test_read_with_quality(2, "ATCGATCG", vec![30, 35, 25, 40, 30, 35, 25, 40]), // Mixed quality
        ];

        for read in reads {
            chunk.add_read(read).unwrap();
        }

        chunk.finalize();

        // All reads contribute to coverage regardless of quality
        // (In a more sophisticated implementation, quality might weight the contribution)
        let total_reads = chunk.processing_stats.reads_processed;
        assert_eq!(total_reads, 3);

        // Verify nodes were created and have appropriate coverage
        assert!(!chunk.graph_fragment.nodes.is_empty());
        let high_coverage_nodes: Vec<_> = chunk
            .graph_fragment
            .nodes
            .values()
            .filter(|n| n.coverage > 2) // Should have coverage from multiple identical reads
            .collect();
        assert!(!high_coverage_nodes.is_empty());
    }
}

#[cfg(test)]
mod complexity_scoring_tests {
    use super::*;

    #[test]
    fn test_sequence_complexity_extremes() {
        // Test minimum complexity (homopolymer)
        let min_complexity = calculate_sequence_complexity("AAAAAAAAAA");
        assert!(
            min_complexity < 0.1,
            "Homopolymer should have very low complexity"
        );

        // Test maximum complexity (perfectly balanced)
        let max_complexity = calculate_sequence_complexity("ATCGATCGATCG");
        assert!(
            max_complexity > 0.9,
            "Balanced sequence should have high complexity"
        );

        // Test intermediate complexity
        let mid_complexity = calculate_sequence_complexity("AAAATCGATCG");
        assert!(mid_complexity > min_complexity && mid_complexity < max_complexity);
    }

    #[test]
    fn test_complexity_calculation_mathematical_properties() {
        // Empty sequence
        assert_eq!(calculate_sequence_complexity(""), 0.0);

        // Single character
        assert_eq!(calculate_sequence_complexity("A"), 0.0);

        // Two different characters (should have complexity > 0)
        let two_char_complexity = calculate_sequence_complexity("AT");
        assert!(two_char_complexity > 0.0 && two_char_complexity <= 1.0);

        // Three different characters
        let three_char_complexity = calculate_sequence_complexity("ATC");
        assert!(three_char_complexity > two_char_complexity);

        // All four characters equally represented (maximum entropy)
        let four_char_complexity = calculate_sequence_complexity("ATCGATCGATCG");
        assert!(four_char_complexity > three_char_complexity);
    }

    #[test]
    fn test_complexity_invariant_to_case() {
        let upper = calculate_sequence_complexity("ATCGATCG");
        let lower = calculate_sequence_complexity("atcgatcg");
        let mixed = calculate_sequence_complexity("AtCgAtCg");

        assert_eq!(upper, lower);
        assert_eq!(upper, mixed);
    }

    #[test]
    fn test_complexity_with_different_lengths() {
        let short = "ATCG";
        let medium = "ATCGATCGATCG";
        let long = "ATCGATCGATCGATCGATCGATCG";

        let short_complexity = calculate_sequence_complexity(short);
        let medium_complexity = calculate_sequence_complexity(medium);
        let long_complexity = calculate_sequence_complexity(long);

        // For perfectly balanced sequences, complexity should be similar regardless of length
        let tolerance = 0.1;
        assert!((short_complexity - medium_complexity).abs() < tolerance);
        assert!((medium_complexity - long_complexity).abs() < tolerance);
    }

    #[test]
    fn test_node_complexity_calculation() {
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let mut node = GraphNode::new(kmer, 4);

        // Test with balanced sequence
        node.calculate_complexity("ATCGATCGATCG");
        let balanced_complexity = node.complexity_score;
        assert!(balanced_complexity > 0.8);

        // Test with unbalanced sequence
        node.calculate_complexity("AAAAAAAAAAA");
        let unbalanced_complexity = node.complexity_score;
        assert!(unbalanced_complexity < 0.2);

        assert!(balanced_complexity > unbalanced_complexity);
    }
}

#[cfg(test)]
mod gc_content_analysis_tests {
    use super::*;

    #[test]
    fn test_gc_content_calculation_accuracy() {
        // Test cases with known GC content
        let test_cases = vec![
            ("AAAA", 0.0),      // 0% GC
            ("GGGG", 1.0),      // 100% GC
            ("ATCG", 0.5),      // 50% GC
            ("AATGCC", 0.5),    // 3/6 = 50% GC
            ("ATATATATA", 0.0), // 0% GC (all AT)
            ("GCGCGCGC", 1.0),  // 100% GC
        ];

        for (sequence, expected) in test_cases {
            let actual = calculate_gc_content(sequence);
            assert!(
                (actual - expected).abs() < 1e-10,
                "GC content of '{sequence}' should be {expected}, got {actual}"
            );
        }
    }

    #[test]
    fn test_gc_content_case_insensitive() {
        let sequences = vec!["atcg", "ATCG", "AtCg", "atCG"];

        for seq in sequences {
            let gc_content = calculate_gc_content(seq);
            assert_eq!(gc_content, 0.5, "GC content should be case insensitive");
        }
    }

    #[test]
    fn test_gc_content_with_assembly_stats() {
        let fragment = GraphFragment::default();
        let mut assembly_graph = AssemblyGraph::from_fragment(fragment);

        // Create contigs with known GC content
        let contigs = vec![
            Contig {
                id: 0,
                sequence: "AAAAAAAAAA".to_string(), // 0% GC
                coverage: 10.0,
                length: 10,
                node_path: vec![],
                contig_type: ContigType::Linear,
            },
            Contig {
                id: 1,
                sequence: "GGGGGGGGGG".to_string(), // 100% GC
                coverage: 15.0,
                length: 10,
                node_path: vec![],
                contig_type: ContigType::Linear,
            },
        ];

        assembly_graph.contigs = contigs;
        assembly_graph.calculate_assembly_stats();

        // Overall GC content should be 50% (10 A/T bases + 10 G/C bases)
        assert_eq!(assembly_graph.assembly_stats.gc_content, 0.5);
    }

    #[test]
    fn test_gc_content_empty_sequences() {
        assert_eq!(calculate_gc_content(""), 0.0);

        // Assembly with empty contigs
        let fragment = GraphFragment::default();
        let mut assembly_graph = AssemblyGraph::from_fragment(fragment);
        assembly_graph.contigs = vec![Contig {
            id: 0,
            sequence: "".to_string(),
            coverage: 0.0,
            length: 0,
            node_path: vec![],
            contig_type: ContigType::Linear,
        }];

        assembly_graph.calculate_assembly_stats();
        assert_eq!(assembly_graph.assembly_stats.gc_content, 0.0);
    }
}

#[cfg(test)]
mod assembly_statistics_tests {
    use super::*;

    #[test]
    fn test_n50_calculation() {
        let fragment = GraphFragment::default();
        let mut assembly_graph = AssemblyGraph::from_fragment(fragment);

        // Create contigs with specific lengths for N50 calculation
        let lengths = [100, 200, 300, 400, 500]; // Total: 1500, N50 should be 400
        let mut contigs = Vec::new();

        for (i, length) in lengths.iter().enumerate() {
            contigs.push(Contig {
                id: i,
                sequence: "A".repeat(*length),
                coverage: 10.0,
                length: *length,
                node_path: vec![],
                contig_type: ContigType::Linear,
            });
        }

        assembly_graph.contigs = contigs;
        assembly_graph.calculate_assembly_stats();

        assert_eq!(assembly_graph.assembly_stats.total_length, 1500);
        assert_eq!(assembly_graph.assembly_stats.num_contigs, 5);
        assert_eq!(assembly_graph.assembly_stats.largest_contig, 500);

        // N50: cumulative length reaches 750 (half of 1500) at contig of length 400
        // Sorted lengths: [500, 400, 300, 200, 100]
        // Cumulative: 500, 900 (>= 750), so N50 = 400
        assert_eq!(assembly_graph.assembly_stats.n50, 400);
    }

    #[test]
    fn test_n50_edge_cases() {
        let fragment = GraphFragment::default();
        let mut assembly_graph = AssemblyGraph::from_fragment(fragment);

        // Single contig
        assembly_graph.contigs = vec![Contig {
            id: 0,
            sequence: "A".repeat(1000),
            coverage: 10.0,
            length: 1000,
            node_path: vec![],
            contig_type: ContigType::Linear,
        }];

        assembly_graph.calculate_assembly_stats();
        assert_eq!(assembly_graph.assembly_stats.n50, 1000);

        // Two equal contigs
        assembly_graph.contigs = vec![
            Contig {
                id: 0,
                sequence: "A".repeat(500),
                coverage: 10.0,
                length: 500,
                node_path: vec![],
                contig_type: ContigType::Linear,
            },
            Contig {
                id: 1,
                sequence: "T".repeat(500),
                coverage: 10.0,
                length: 500,
                node_path: vec![],
                contig_type: ContigType::Linear,
            },
        ];

        assembly_graph.calculate_assembly_stats();
        assert_eq!(assembly_graph.assembly_stats.n50, 500);
    }

    #[test]
    fn test_coverage_statistics_calculation() {
        let fragment = GraphFragment::default();
        let mut assembly_graph = AssemblyGraph::from_fragment(fragment);

        let contigs = vec![
            Contig {
                id: 0,
                sequence: "A".repeat(100),
                coverage: 10.0,
                length: 100,
                node_path: vec![],
                contig_type: ContigType::Linear,
            },
            Contig {
                id: 1,
                sequence: "T".repeat(200),
                coverage: 20.0,
                length: 200,
                node_path: vec![],
                contig_type: ContigType::Linear,
            },
            Contig {
                id: 2,
                sequence: "G".repeat(300),
                coverage: 30.0,
                length: 300,
                node_path: vec![],
                contig_type: ContigType::Linear,
            },
        ];

        assembly_graph.contigs = contigs;
        assembly_graph.calculate_assembly_stats();

        // Weighted average coverage: (10*100 + 20*200 + 30*300) / (100+200+300)
        // = (1000 + 4000 + 9000) / 600 = 14000 / 600 = 23.33...
        let expected_coverage = 14000.0 / 600.0;
        assert!((assembly_graph.assembly_stats.coverage_mean - expected_coverage).abs() < 1e-10);
    }

    #[test]
    fn test_empty_assembly_statistics() {
        let fragment = GraphFragment::default();
        let mut assembly_graph = AssemblyGraph::from_fragment(fragment);

        // No contigs
        assembly_graph.calculate_assembly_stats();

        assert_eq!(assembly_graph.assembly_stats.total_length, 0);
        assert_eq!(assembly_graph.assembly_stats.num_contigs, 0);
        assert_eq!(assembly_graph.assembly_stats.n50, 0);
        assert_eq!(assembly_graph.assembly_stats.largest_contig, 0);
        assert_eq!(assembly_graph.assembly_stats.gc_content, 0.0);
        assert_eq!(assembly_graph.assembly_stats.coverage_mean, 0.0);
    }
}

#[cfg(test)]
mod statistical_validation_tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_read_position_tracking_accuracy() {
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let mut node = GraphNode::new(kmer, 4);

        // Add read positions from different reads
        let positions = vec![
            (1, 10, Strand::Forward),
            (2, 15, Strand::Reverse),
            (3, 20, Strand::Forward),
            (1, 25, Strand::Forward), // Same read, different position
        ];

        for (read_id, pos, strand) in positions {
            node.add_read_position(read_id, pos, strand);
        }

        assert_eq!(node.read_positions.len(), 4);
        assert_eq!(node.coverage, 5); // Initial 1 + 4 additions

        // Verify all positions were recorded
        let read_ids: HashSet<usize> = node.read_positions.iter().map(|rp| rp.read_id).collect();
        assert!(read_ids.contains(&1));
        assert!(read_ids.contains(&2));
        assert!(read_ids.contains(&3));
        assert_eq!(read_ids.len(), 3); // Three unique reads
    }

    #[test]
    fn test_edge_weight_accumulation() {
        let mut edge = GraphEdge::new(123, 456, 3);

        // Simulate multiple reads supporting this edge
        let supporting_reads = vec![1, 2, 3, 4, 5];
        for read_id in supporting_reads {
            edge.add_support(read_id);
        }

        assert_eq!(edge.weight, 5);
        assert_eq!(edge.supporting_reads.len(), 5);

        // Adding duplicate support shouldn't change weight
        edge.add_support(1);
        assert_eq!(edge.weight, 5);

        // Confidence should increase with more support
        assert!(edge.confidence > 0.1);
        assert!(edge.confidence <= 1.0);
    }

    #[test]
    fn test_processing_statistics_accuracy() {
        let mut chunk = AssemblyChunk::new(0, 3);

        let reads = vec![
            create_test_read_with_quality(0, "ATCGATCG", vec![30; 8]),
            create_test_read_with_quality(1, "TCGATCGA", vec![35; 8]),
        ];

        for read in reads {
            chunk.add_read(read).unwrap();
        }

        chunk.finalize();

        // Verify statistics accuracy
        assert_eq!(chunk.processing_stats.reads_processed, 2);
        assert!(chunk.processing_stats.nodes_created > 0);
        assert!(chunk.processing_stats.edges_created > 0);
        assert!(chunk.processing_stats.minimizers_found > 0);
        assert!(chunk.processing_stats.processing_time_ms >= 0);
        assert!(chunk.processing_stats.memory_usage_bytes > 0);

        // Verify relationships between statistics
        assert!(chunk.processing_stats.nodes_created <= chunk.processing_stats.minimizers_found);
        assert_eq!(
            chunk.graph_fragment.nodes.len(),
            chunk.processing_stats.nodes_created
        );
        assert_eq!(
            chunk.graph_fragment.edges.len(),
            chunk.processing_stats.edges_created
        );
    }

    #[test]
    fn test_memory_usage_estimation() {
        let mut chunk = AssemblyChunk::new(0, 4);

        let initial_memory = chunk.estimate_memory_usage();

        // Add some reads
        let reads = vec![
            create_test_read_with_quality(0, "ATCGATCGATCG", vec![30; 12]),
            create_test_read_with_quality(1, "TCGATCGATCGA", vec![30; 12]),
        ];

        for read in reads {
            chunk.add_read(read).unwrap();
        }

        chunk.finalize();

        let final_memory = chunk.processing_stats.memory_usage_bytes;

        // Memory usage should increase after adding reads
        assert!(final_memory > initial_memory);
        assert!(final_memory > 0);

        // Should account for nodes, edges, and read data
        let expected_minimum = chunk.graph_fragment.nodes.len() * std::mem::size_of::<GraphNode>()
            + chunk.graph_fragment.edges.len() * std::mem::size_of::<GraphEdge>();

        assert!(final_memory >= expected_minimum);
    }
}
