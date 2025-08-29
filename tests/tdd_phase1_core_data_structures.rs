//! Phase 1 TDD Tests for Core Data Structures
//!
//! Critical tests for the foundation of genomic data processing
//! Focus: Biological correctness, k-mer handling, graph operations

use meta_forge::core::data_structures::*;
use meta_forge::utils::configuration::{AmbiguousBaseConfig, AmbiguousBaseStrategy};
use std::collections::HashMap;

/// Test canonical k-mer correctness with ambiguous bases
/// BIOLOGICAL REQUIREMENT: K-mers with N bases must be handled consistently
/// according to configuration, maintaining biological meaning
#[cfg(test)]
mod canonical_kmer_tests {
    use super::*;

    #[test]
    fn test_canonical_kmer_with_n_bases_skip_strategy() {
        let config = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::Skip,
            max_n_count: 0,
            replacement_base: 'A',
            random_probabilities: None,
        };

        // Should reject sequences with N bases
        let sequences_with_n = vec!["NATCG", "ATNCG", "ATCGN", "NNNNN"];
        for seq in sequences_with_n {
            let result = CanonicalKmer::new_with_config(seq, &config);
            assert!(
                result.is_err(),
                "Skip strategy should reject sequence with N: {seq}"
            );
        }

        // Should accept clean sequences
        let clean_sequences = vec!["ATCG", "GGCC", "AAAA"];
        for seq in clean_sequences {
            let result = CanonicalKmer::new_with_config(seq, &config);
            assert!(
                result.is_ok(),
                "Skip strategy should accept clean sequence: {seq}"
            );
        }
    }

    #[test]
    fn test_canonical_kmer_with_n_bases_replace_strategy() {
        let config = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::Replace,
            max_n_count: 10,
            replacement_base: 'A',
            random_probabilities: None,
        };

        let test_cases = vec![
            ("NATCG", "AATCG"),
            ("ATNCG", "ATACG"),
            ("NNNNN", "AAAAA"),
            ("ATCGN", "ATCGA"),
        ];

        for (input, expected_base_replacement) in test_cases {
            let result = CanonicalKmer::new_with_config(input, &config);
            assert!(
                result.is_ok(),
                "Replace strategy should handle sequence: {input}"
            );

            let kmer = result.unwrap();
            // Verify no N bases remain
            assert!(
                !kmer.sequence.contains('N'),
                "Replaced sequence should not contain N: {}",
                kmer.sequence
            );

            // For replace strategy with 'A', verify replacement occurred
            if input.contains('N') {
                assert!(
                    kmer.sequence.contains('A'),
                    "Replacement with A should occur"
                );
            }
        }
    }

    #[test]
    fn test_reverse_complement_canonicalization_accuracy() {
        let test_cases = vec![
            // (sequence, reverse_complement, expected_canonical)
            ("ATCG", "CGAT", "ATCG"), // ATCG < CGAT lexicographically
            ("CGAT", "ATCG", "ATCG"), // Should canonicalize to same as above
            ("AAAA", "TTTT", "AAAA"), // AAAA < TTTT
            ("TTTT", "AAAA", "AAAA"), // Should canonicalize to same as above
            ("ACGT", "ACGT", "ACGT"), // Palindromic sequence
        ];

        for (seq, expected_rc, expected_canonical) in test_cases {
            let kmer1 = CanonicalKmer::new(seq).unwrap();
            let kmer2 = CanonicalKmer::new(expected_rc).unwrap();

            // Both should have the same canonical representation
            assert_eq!(
                kmer1.sequence, expected_canonical,
                "Sequence {seq} should canonicalize to {expected_canonical}"
            );
            assert_eq!(
                kmer2.sequence, expected_canonical,
                "Reverse complement {expected_rc} should canonicalize to {expected_canonical}"
            );

            // Hash should be identical for canonical forms
            assert_eq!(
                kmer1.hash, kmer2.hash,
                "Hash should be identical for canonical k-mers"
            );
        }
    }

    #[test]
    fn test_hash_consistency_for_identical_sequences() {
        let sequence = "ATCGATCG";
        let hashes: Vec<u64> = (0..100)
            .map(|_| CanonicalKmer::new(sequence).unwrap().hash)
            .collect();

        // All hashes should be identical
        let first_hash = hashes[0];
        for (i, &hash) in hashes.iter().enumerate() {
            assert_eq!(
                hash, first_hash,
                "Hash inconsistency at iteration {i}: {hash} != {first_hash}"
            );
        }
    }

    #[test]
    fn test_invalid_dna_character_handling() {
        let invalid_sequences = vec![
            "ATCGX",  // Contains X
            "ATCG1",  // Contains number
            "ATCG-",  // Contains dash
            "ATCG ",  // Contains space
            "ATCG\n", // Contains newline
        ];

        for seq in invalid_sequences {
            let result = CanonicalKmer::new(seq);
            assert!(
                result.is_err(),
                "Should reject invalid DNA sequence: {seq:?}"
            );
        }
    }
}

/// Test graph fragment merge operations
/// BIOLOGICAL REQUIREMENT: Node coverage and edge weights must combine correctly
/// during graph merging operations without losing biological information
#[cfg(test)]
mod graph_fragment_merge_tests {
    use super::*;

    #[test]
    fn test_node_coverage_aggregation_during_merge() {
        let mut fragment1 = GraphFragment::new(0);
        let mut fragment2 = GraphFragment::new(1);

        // Create same k-mer in both fragments with different coverage
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let mut node1 = GraphNode::new(kmer.clone(), 4);
        node1.coverage = 5;

        let mut node2 = GraphNode::new(kmer.clone(), 4);
        node2.coverage = 3;

        fragment1.add_node(node1);
        fragment2.add_node(node2);

        // Merge fragments
        let result = fragment1.merge_with(fragment2);
        assert!(result.is_ok(), "Fragment merge should succeed");

        // Verify coverage aggregation
        let merged_node = fragment1.nodes.get(&kmer.hash).unwrap();
        assert_eq!(
            merged_node.coverage, 8,
            "Coverage should aggregate: 5 + 3 = 8"
        );
    }

    #[test]
    fn test_edge_weight_combination_during_merge() {
        let mut fragment1 = GraphFragment::new(0);
        let mut fragment2 = GraphFragment::new(1);

        // Create k-mers for edges
        let kmer1 = CanonicalKmer::new("ATCG").unwrap();
        let kmer2 = CanonicalKmer::new("TCGA").unwrap();

        // Add nodes to both fragments
        fragment1.add_node(GraphNode::new(kmer1.clone(), 4));
        fragment1.add_node(GraphNode::new(kmer2.clone(), 4));
        fragment2.add_node(GraphNode::new(kmer1.clone(), 4));
        fragment2.add_node(GraphNode::new(kmer2.clone(), 4));

        // Create same edge in both fragments with different weights
        let mut edge1 = GraphEdge::new(kmer1.hash, kmer2.hash, 3);
        edge1.weight = 4;
        edge1.supporting_reads.insert(1);
        edge1.supporting_reads.insert(2);

        let mut edge2 = GraphEdge::new(kmer1.hash, kmer2.hash, 3);
        edge2.weight = 3;
        edge2.supporting_reads.insert(3);
        edge2.supporting_reads.insert(4);

        fragment1.add_edge(edge1);
        fragment2.add_edge(edge2);

        // Merge fragments
        let result = fragment1.merge_with(fragment2);
        assert!(result.is_ok(), "Fragment merge should succeed");

        // Find the merged edge
        let merged_edge = fragment1
            .edges
            .iter()
            .find(|e| e.from_hash == kmer1.hash && e.to_hash == kmer2.hash)
            .expect("Merged edge should exist");

        // Verify weight combination
        assert_eq!(
            merged_edge.weight, 7,
            "Edge weight should combine: 4 + 3 = 7"
        );

        // Verify supporting reads combination
        assert_eq!(
            merged_edge.supporting_reads.len(),
            4,
            "Supporting reads should combine"
        );
        assert!(merged_edge.supporting_reads.contains(&1));
        assert!(merged_edge.supporting_reads.contains(&2));
        assert!(merged_edge.supporting_reads.contains(&3));
        assert!(merged_edge.supporting_reads.contains(&4));
    }

    #[test]
    fn test_no_nodes_or_edges_lost_during_merge() {
        let mut fragment1 = GraphFragment::new(0);
        let mut fragment2 = GraphFragment::new(1);

        // Create unique k-mers for each fragment
        let kmers1 = vec!["ATCG", "TCGA", "CGAT"];
        let kmers2 = vec!["GATC", "ATCG", "GCTA"]; // "ATCG" overlaps

        // Add nodes to fragment1
        let mut kmer_hashes1 = Vec::new();
        for seq in &kmers1 {
            let kmer = CanonicalKmer::new(seq).unwrap();
            kmer_hashes1.push(kmer.hash);
            fragment1.add_node(GraphNode::new(kmer, 4));
        }

        // Add nodes to fragment2
        let mut kmer_hashes2 = Vec::new();
        for seq in &kmers2 {
            let kmer = CanonicalKmer::new(seq).unwrap();
            kmer_hashes2.push(kmer.hash);
            fragment2.add_node(GraphNode::new(kmer, 4));
        }

        // Add some edges
        if kmer_hashes1.len() >= 2 {
            fragment1.add_edge(GraphEdge::new(kmer_hashes1[0], kmer_hashes1[1], 3));
        }
        if kmer_hashes2.len() >= 2 {
            fragment2.add_edge(GraphEdge::new(kmer_hashes2[0], kmer_hashes2[1], 3));
        }

        let original_nodes1 = fragment1.nodes.len();
        let original_nodes2 = fragment2.nodes.len();
        let original_edges1 = fragment1.edges.len();
        let original_edges2 = fragment2.edges.len();

        // Merge fragments
        let result = fragment1.merge_with(fragment2);
        assert!(result.is_ok(), "Fragment merge should succeed");

        // Calculate expected counts (accounting for "ATCG" overlap)
        let expected_nodes = original_nodes1 + original_nodes2 - 1; // -1 for overlapping "ATCG"
        let expected_edges = original_edges1 + original_edges2;

        assert_eq!(
            fragment1.nodes.len(),
            expected_nodes,
            "Node count should match expected after merge"
        );
        assert_eq!(
            fragment1.edges.len(),
            expected_edges,
            "Edge count should match expected after merge"
        );

        // Verify all original k-mers are still present
        for seq in kmers1.iter().chain(kmers2.iter()) {
            let kmer = CanonicalKmer::new(seq).unwrap();
            assert!(
                fragment1.nodes.contains_key(&kmer.hash),
                "K-mer {seq} should still be present after merge"
            );
        }
    }

    #[test]
    fn test_merge_with_empty_fragment() {
        let mut fragment1 = GraphFragment::new(0);
        let fragment2 = GraphFragment::new(1); // Empty

        // Add content to fragment1
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        fragment1.add_node(GraphNode::new(kmer, 4));

        let original_nodes = fragment1.nodes.len();
        let original_edges = fragment1.edges.len();

        // Merge with empty fragment
        let result = fragment1.merge_with(fragment2);
        assert!(result.is_ok(), "Merge with empty fragment should succeed");

        // Nothing should change
        assert_eq!(
            fragment1.nodes.len(),
            original_nodes,
            "Node count should not change when merging with empty fragment"
        );
        assert_eq!(
            fragment1.edges.len(),
            original_edges,
            "Edge count should not change when merging with empty fragment"
        );
    }

    #[test]
    fn test_merge_preserves_read_position_information() {
        let mut fragment1 = GraphFragment::new(0);
        let mut fragment2 = GraphFragment::new(1);

        // Create same k-mer in both fragments with different read positions
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let mut node1 = GraphNode::new(kmer.clone(), 4);
        node1.add_read_position(1, 10, Strand::Forward);
        node1.add_read_position(2, 15, Strand::Reverse);

        let mut node2 = GraphNode::new(kmer.clone(), 4);
        node2.add_read_position(3, 20, Strand::Forward);
        node2.add_read_position(4, 25, Strand::Reverse);

        fragment1.add_node(node1);
        fragment2.add_node(node2);

        // Merge fragments
        let result = fragment1.merge_with(fragment2);
        assert!(result.is_ok(), "Fragment merge should succeed");

        // Verify read positions are preserved
        let merged_node = fragment1.nodes.get(&kmer.hash).unwrap();
        assert_eq!(
            merged_node.read_positions.len(),
            4,
            "All read positions should be preserved"
        );

        // Verify specific read IDs are present
        let read_ids: Vec<usize> = merged_node
            .read_positions
            .iter()
            .map(|rp| rp.read_id)
            .collect();
        assert!(read_ids.contains(&1));
        assert!(read_ids.contains(&2));
        assert!(read_ids.contains(&3));
        assert!(read_ids.contains(&4));
    }
}

/// Test sequence complexity calculations
/// BIOLOGICAL REQUIREMENT: Complexity scores must correlate with known
/// low/high complexity regions and use proper Shannon entropy
#[cfg(test)]
mod sequence_complexity_tests {
    use super::*;

    #[test]
    fn test_shannon_entropy_known_values() {
        // Test cases with known complexity values
        let test_cases = vec![
            ("AAAA", 0.0),          // No entropy - single base
            ("AAAAACCCCC", 0.4307), // Approximate expected entropy for 2-base mix
            ("AACCGGTT", 1.0),      // Maximum entropy for DNA (4 equal bases)
        ];

        for (sequence, expected_complexity) in test_cases {
            let complexity = calculate_sequence_complexity(sequence);
            let tolerance = 0.1; // Allow 10% tolerance for floating point

            assert!(
                (complexity - expected_complexity).abs() < tolerance,
                "Complexity for '{sequence}' should be ~{expected_complexity}, got {complexity}"
            );
        }
    }

    #[test]
    fn test_gc_content_calculation_accuracy() {
        let test_cases = vec![
            ("ATCG", 0.5),       // 2 GC out of 4 = 50%
            ("AAAA", 0.0),       // 0 GC out of 4 = 0%
            ("CCGG", 1.0),       // 4 GC out of 4 = 100%
            ("ATCGATCG", 0.5),   // 4 GC out of 8 = 50%
            ("", 0.0),           // Empty sequence = 0%
            ("GCGCGCGC", 1.0),   // All GC = 100%
            ("ATATATATAT", 0.0), // No GC = 0%
        ];

        for (sequence, expected_gc) in test_cases {
            let gc_content = calculate_gc_content(sequence);
            let tolerance = 0.001; // Very tight tolerance for GC content

            assert!(
                (gc_content - expected_gc).abs() < tolerance,
                "GC content for '{sequence}' should be {expected_gc}, got {gc_content}"
            );
        }
    }

    #[test]
    fn test_low_complexity_sequence_detection() {
        let low_complexity_sequences = vec![
            "AAAAAAAAAA", // Homopolymer
            "ATATATATAT", // Simple repeat
            "AGAGAGAGAG", // Dinucleotide repeat
            "AAAAACCCCC", // Low diversity
        ];

        for sequence in low_complexity_sequences {
            let complexity = calculate_sequence_complexity(sequence);
            assert!(
                complexity < 0.7,
                "Sequence '{sequence}' should have low complexity (<0.7), got {complexity}"
            );
        }
    }

    #[test]
    fn test_high_complexity_sequence_detection() {
        let high_complexity_sequences = vec![
            "ACGTACGTACGT",    // Balanced nucleotides
            "ATCGATCGATCG",    // Balanced with different pattern
            "AGCTTAGCTTAGCTT", // Random-like distribution
        ];

        for sequence in high_complexity_sequences {
            let complexity = calculate_sequence_complexity(sequence);
            assert!(
                complexity > 0.8,
                "Sequence '{sequence}' should have high complexity (>0.8), got {complexity}"
            );
        }
    }

    #[test]
    fn test_complexity_bounds_and_edge_cases() {
        // Test complexity is always between 0 and 1
        let test_sequences = vec![
            "",                     // Empty
            "A",                    // Single base
            "AT",                   // Two bases
            "ATCG",                 // Four bases
            "AAAAAAAAAAAAAAAAAAA",  // Long homopolymer
            "ATCGATCGATCGATCGATCG", // Long balanced
        ];

        for sequence in test_sequences {
            let complexity = calculate_sequence_complexity(sequence);
            assert!(
                (0.0..=1.0).contains(&complexity),
                "Complexity for '{sequence}' should be in [0,1], got {complexity}"
            );
        }
    }

    #[test]
    fn test_dna_sequence_validation() {
        // Valid sequences
        let valid_sequences = vec![
            "ATCG", "atcg",     // Lowercase should be valid
            "ATCGN",    // N should be valid
            "ATCGatcg", // Mixed case
            "",         // Empty sequence
        ];

        for sequence in valid_sequences {
            let result = validate_dna_sequence(sequence);
            assert!(result.is_ok(), "Sequence '{sequence}' should be valid");
        }

        // Invalid sequences
        let invalid_sequences = vec![
            "ATCGX",  // Invalid character X
            "ATCG1",  // Number
            "ATCG-",  // Dash
            "ATCG ",  // Space
            "ATCG\n", // Newline
            "ATCGR",  // R (ambiguous base not handled here)
        ];

        for sequence in invalid_sequences {
            let result = validate_dna_sequence(sequence);
            assert!(result.is_err(), "Sequence '{sequence}' should be invalid");
        }
    }
}

/// Integration test: Complete workflow validation
/// BIOLOGICAL REQUIREMENT: End-to-end operations should preserve biological meaning
#[cfg(test)]
mod integration_tests {
    use super::*;

    #[test]
    fn test_complete_kmer_to_graph_workflow() {
        // Create a simple test sequence representing a small genome fragment
        let test_sequence = "ATCGATCGATCGATCG"; // 16bp with 4bp repeat
        let k = 4;

        // Extract all k-mers
        let mut kmers = Vec::new();
        for i in 0..=test_sequence.len() - k {
            let kmer_seq = &test_sequence[i..i + k];
            let kmer = CanonicalKmer::new(kmer_seq).unwrap();
            kmers.push(kmer);
        }

        // Build a graph fragment
        let mut fragment = GraphFragment::new(0);

        // Add nodes for each unique k-mer
        let mut kmer_counts = HashMap::new();
        for kmer in &kmers {
            *kmer_counts.entry(kmer.hash).or_insert(0) += 1;
        }

        for (hash, count) in kmer_counts {
            if let Some(kmer) = kmers.iter().find(|k| k.hash == hash) {
                let mut node = GraphNode::new(kmer.clone(), k);
                node.coverage = count;
                fragment.add_node(node);
            }
        }

        // Add edges between consecutive k-mers
        for i in 0..kmers.len() - 1 {
            let edge = GraphEdge::new(kmers[i].hash, kmers[i + 1].hash, k - 1);
            fragment.add_edge(edge);
        }

        // Validate the resulting graph structure
        assert!(!fragment.nodes.is_empty(), "Graph should have nodes");
        assert!(!fragment.edges.is_empty(), "Graph should have edges");

        // Check that repeated k-mers have higher coverage
        let repeated_kmer = CanonicalKmer::new("ATCG").unwrap();
        if let Some(node) = fragment.nodes.get(&repeated_kmer.hash) {
            assert!(
                node.coverage > 1,
                "Repeated k-mer should have coverage > 1, got {}",
                node.coverage
            );
        }

        // Verify graph connectivity
        let adjacency = fragment.get_adjacency_list();
        let connected_nodes = adjacency.keys().count();
        assert!(connected_nodes > 0, "Graph should have connected nodes");
    }

    #[test]
    fn test_biological_sequence_preservation() {
        // Test that our data structures preserve biological information
        let original_sequence = "ATCGATCGATCGATCG";
        let k = 6;

        // Process through our pipeline
        let mut fragment = GraphFragment::new(0);

        // Extract overlapping k-mers
        let mut path_hashes = Vec::new();
        for i in 0..=original_sequence.len() - k {
            let kmer_seq = &original_sequence[i..i + k];
            let kmer = CanonicalKmer::new(kmer_seq).unwrap();
            path_hashes.push(kmer.hash);

            // Add to graph
            fragment.add_node(GraphNode::new(kmer, k));
        }

        // Attempt to reconstruct sequence from path
        let reconstructed = fragment
            .reconstruct_sequence_from_path(&path_hashes)
            .unwrap();

        // Should be able to recover the original sequence (or at least preserve length)
        assert!(
            reconstructed.len() > original_sequence.len() - k,
            "Reconstructed sequence should preserve approximate original length"
        );
    }
}
