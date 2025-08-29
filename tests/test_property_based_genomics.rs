//! Property-based tests for genomic data structures
//! Tests invariants and properties that should hold across all valid inputs

use meta_forge::assembly::adaptive_k::*;
use meta_forge::core::data_structures::*;
use std::collections::{HashMap, HashSet};

#[cfg(test)]
pub mod canonical_kmer_properties {
    use super::*;

    fn generate_valid_dna_char() -> char {
        match fastrand::usize(0..4) {
            0 => 'A',
            1 => 'T',
            2 => 'C',
            _ => 'G',
        }
    }

    pub fn generate_random_dna_sequence(length: usize) -> String {
        (0..length).map(|_| generate_valid_dna_char()).collect()
    }

    fn generate_reverse_complement(seq: &str) -> String {
        seq.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => c,
            })
            .collect()
    }

    #[test]
    fn property_canonical_kmer_idempotent() {
        // Property: Canonicalizing a k-mer multiple times should give the same result
        for _ in 0..100 {
            let seq = generate_random_dna_sequence(fastrand::usize(1..20));

            if let (Ok(kmer1), Ok(kmer2)) = (CanonicalKmer::new(&seq), CanonicalKmer::new(&seq)) {
                assert_eq!(kmer1.sequence, kmer2.sequence);
                assert_eq!(kmer1.hash, kmer2.hash);
                assert_eq!(kmer1.is_canonical, kmer2.is_canonical);
            }
        }
    }

    #[test]
    fn property_reverse_complement_canonical_equivalence() {
        // Property: A k-mer and its reverse complement should have the same canonical form
        for _ in 0..100 {
            let seq = generate_random_dna_sequence(fastrand::usize(3..15));
            let rev_comp = generate_reverse_complement(&seq);

            if let (Ok(kmer1), Ok(kmer2)) =
                (CanonicalKmer::new(&seq), CanonicalKmer::new(&rev_comp))
            {
                assert_eq!(
                    kmer1.sequence, kmer2.sequence,
                    "Canonical form should be same for {seq} and {rev_comp}"
                );
                assert_eq!(kmer1.hash, kmer2.hash);
            }
        }
    }

    #[test]
    fn property_canonical_kmer_ordering() {
        // Property: Canonical k-mer should be lexicographically <= its reverse complement
        for _ in 0..100 {
            let seq = generate_random_dna_sequence(fastrand::usize(3..15));

            if let Ok(kmer) = CanonicalKmer::new(&seq) {
                let rev_comp = generate_reverse_complement(&kmer.sequence);
                assert!(
                    kmer.sequence <= rev_comp,
                    "Canonical k-mer '{}' should be <= reverse complement '{}'",
                    kmer.sequence,
                    rev_comp
                );
            }
        }
    }

    #[test]
    fn property_hash_consistency() {
        // Property: Same sequence should always produce same hash
        let test_sequences = vec!["ATCG", "GGCC", "ATAT", "CGCG"];

        for seq in test_sequences {
            let mut hashes = Vec::new();

            for _ in 0..10 {
                if let Ok(kmer) = CanonicalKmer::new(seq) {
                    hashes.push(kmer.hash);
                }
            }

            // All hashes should be identical
            for hash in &hashes {
                assert_eq!(
                    *hash, hashes[0],
                    "Hash should be consistent for sequence {seq}"
                );
            }
        }
    }

    #[test]
    fn property_case_insensitive_canonicalization() {
        // Property: Case should not affect canonicalization
        for _ in 0..50 {
            let seq = generate_random_dna_sequence(fastrand::usize(3..12));
            let seq_lower = seq.to_lowercase();
            let seq_upper = seq.to_uppercase();

            let mut seq_mixed = String::new();
            for (i, c) in seq.chars().enumerate() {
                if i % 2 == 0 {
                    seq_mixed.push(c.to_ascii_uppercase());
                } else {
                    seq_mixed.push(c.to_ascii_lowercase());
                }
            }

            if let (Ok(kmer_lower), Ok(kmer_upper), Ok(kmer_mixed)) = (
                CanonicalKmer::new(&seq_lower),
                CanonicalKmer::new(&seq_upper),
                CanonicalKmer::new(&seq_mixed),
            ) {
                assert_eq!(kmer_lower.sequence, kmer_upper.sequence);
                assert_eq!(kmer_upper.sequence, kmer_mixed.sequence);
                assert_eq!(kmer_lower.hash, kmer_upper.hash);
                assert_eq!(kmer_upper.hash, kmer_mixed.hash);
            }
        }
    }

    #[test]
    fn property_length_preservation() {
        // Property: Canonical k-mer should preserve original length
        for _ in 0..100 {
            let length = fastrand::usize(1..25);
            let seq = generate_random_dna_sequence(length);

            if let Ok(kmer) = CanonicalKmer::new(&seq) {
                assert_eq!(kmer.len(), length, "Length should be preserved");
                assert_eq!(kmer.sequence.len(), length, "Sequence length should match");
            }
        }
    }
}

#[cfg(test)]
mod minimizer_properties {
    use super::*;

    #[test]
    fn property_minimizer_positions_valid() {
        // Property: All minimizer positions should be valid for the sequence
        for _ in 0..50 {
            let k = fastrand::usize(3..8);
            let w = fastrand::usize(2..6);
            let seq_len = fastrand::usize(k..50);
            let sequence = canonical_kmer_properties::generate_random_dna_sequence(seq_len);

            let extractor = MinimizerExtractor::new(k, w);
            if let Ok(minimizers) = extractor.extract_minimizers(&sequence) {
                for minimizer in &minimizers {
                    assert!(
                        minimizer.position < sequence.len(),
                        "Position {} should be < sequence length {}",
                        minimizer.position,
                        sequence.len()
                    );

                    assert!(
                        minimizer.position + k <= sequence.len(),
                        "Position {} + k {} should be <= sequence length {}",
                        minimizer.position,
                        k,
                        sequence.len()
                    );

                    assert!(
                        minimizer.window_start <= minimizer.window_end,
                        "Window start {} should be <= window end {}",
                        minimizer.window_start,
                        minimizer.window_end
                    );
                }
            }
        }
    }

    #[test]
    fn property_minimizer_k_mer_length() {
        // Property: All minimizer k-mers should have the specified length
        for _ in 0..50 {
            let k = fastrand::usize(3..10);
            let w = fastrand::usize(2..8);
            let sequence =
                canonical_kmer_properties::generate_random_dna_sequence(fastrand::usize(k..100));

            let extractor = MinimizerExtractor::new(k, w);
            if let Ok(minimizers) = extractor.extract_minimizers(&sequence) {
                for minimizer in &minimizers {
                    assert_eq!(
                        minimizer.kmer.len(),
                        k,
                        "Minimizer k-mer should have length {k}"
                    );
                }
            }
        }
    }

    #[test]
    fn property_minimizer_deduplication() {
        // Property: Consecutive identical minimizers should be deduplicated
        for _ in 0..30 {
            let k = 4;
            let w = 3;
            let repeat_count = fastrand::usize(5..15);
            let sequence = "ATCG".repeat(repeat_count); // Repeated pattern

            let extractor = MinimizerExtractor::new(k, w);
            if let Ok(minimizers) = extractor.extract_minimizers(&sequence) {
                // Should have fewer minimizers than possible due to deduplication
                let max_possible = sequence.len() - k + 1;
                assert!(
                    minimizers.len() <= max_possible,
                    "Should have <= {} minimizers, got {}",
                    max_possible,
                    minimizers.len()
                );

                // Check for consecutive duplicates
                for window in minimizers.windows(2) {
                    // Allow consecutive identical if they're from different windows
                    // but positions should be different
                    if window[0].kmer.hash == window[1].kmer.hash {
                        assert!(
                            window[0].position != window[1].position,
                            "Consecutive identical minimizers should have different positions"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn property_minimizer_window_coverage() {
        // Property: Every window should contribute at least one minimizer (if possible)
        for _ in 0..30 {
            let k = 3;
            let w = 4;
            let sequence = canonical_kmer_properties::generate_random_dna_sequence(20);

            let extractor = MinimizerExtractor::new(k, w);
            if let Ok(minimizers) = extractor.extract_minimizers(&sequence) {
                if sequence.len() >= k {
                    assert!(!minimizers.is_empty(), "Should have at least one minimizer");

                    // Check that minimizers span a reasonable portion of the sequence
                    if minimizers.len() > 1 {
                        let first_pos = minimizers.first().unwrap().position;
                        let last_pos = minimizers.last().unwrap().position;
                        let span = last_pos.saturating_sub(first_pos);

                        // Should cover a reasonable span
                        assert!(
                            span <= sequence.len() - k,
                            "Minimizer span should be reasonable"
                        );
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod graph_structure_properties {
    use super::*;

    #[test]
    fn property_graph_fragment_node_edge_consistency() {
        // Property: All edges should reference nodes that exist in the fragment
        for _ in 0..20 {
            let mut fragment = GraphFragment::new(0);
            let node_count = fastrand::usize(3..10);
            let mut node_hashes = Vec::new();

            // Add nodes
            for i in 0..node_count {
                let seq = format!("{}TCG", (b'A' + (i % 4) as u8) as char);
                if let Ok(kmer) = CanonicalKmer::new(&seq) {
                    let node = GraphNode::new(kmer.clone(), 4);
                    node_hashes.push(kmer.hash);
                    fragment.add_node(node);
                }
            }

            // Add edges between random nodes
            for _ in 0..fastrand::usize(1..node_count) {
                if node_hashes.len() >= 2 {
                    let from_idx = fastrand::usize(0..node_hashes.len());
                    let to_idx = fastrand::usize(0..node_hashes.len());
                    let edge = GraphEdge::new(node_hashes[from_idx], node_hashes[to_idx], 3);
                    fragment.add_edge(edge);
                }
            }

            // Check consistency
            for edge in &fragment.edges {
                assert!(
                    fragment.nodes.contains_key(&edge.from_hash),
                    "Edge source node should exist in fragment"
                );
                assert!(
                    fragment.nodes.contains_key(&edge.to_hash),
                    "Edge target node should exist in fragment"
                );
            }
        }
    }

    #[test]
    fn property_coverage_statistics_consistency() {
        // Property: Coverage statistics should be consistent with actual nodes
        for _ in 0..20 {
            let mut fragment = GraphFragment::new(0);
            let mut expected_total = 0u32;
            let mut coverages = Vec::new();

            let node_count = fastrand::usize(1..15);
            for i in 0..node_count {
                let seq = format!(
                    "{}{}CG",
                    (b'A' + (i % 4) as u8) as char,
                    (b'A' + ((i + 1) % 4) as u8) as char
                );

                if let Ok(kmer) = CanonicalKmer::new(&seq) {
                    let mut node = GraphNode::new(kmer, 4);
                    let coverage = fastrand::u32(1..50);
                    node.coverage = coverage;

                    expected_total += coverage;
                    coverages.push(coverage);
                    fragment.add_node(node);
                }
            }

            // Check statistics
            let stats = &fragment.coverage_stats;
            assert_eq!(stats.total_nodes, fragment.nodes.len());

            if !coverages.is_empty() {
                let expected_mean = expected_total as f64 / coverages.len() as f64;
                assert!(
                    (stats.mean_coverage - expected_mean).abs() < 1e-10,
                    "Mean coverage should match calculated value"
                );

                coverages.sort_unstable();
                let expected_median = coverages[coverages.len() / 2];
                assert_eq!(
                    stats.median_coverage, expected_median,
                    "Median coverage should match"
                );
            }
        }
    }

    #[test]
    fn property_adjacency_list_symmetry() {
        // Property: Adjacency list should correctly represent all edges
        for _ in 0..20 {
            let mut fragment = GraphFragment::new(0);
            let mut added_edges = HashSet::new();

            // Add nodes
            let node_seqs = vec!["ATCG", "TCGA", "CGAT", "GATC", "ATGG"];
            let mut node_hashes = Vec::new();

            for seq in &node_seqs {
                if let Ok(kmer) = CanonicalKmer::new(seq) {
                    let node = GraphNode::new(kmer.clone(), 4);
                    node_hashes.push(kmer.hash);
                    fragment.add_node(node);
                }
            }

            // Add random edges
            for _ in 0..fastrand::usize(2..8) {
                if node_hashes.len() >= 2 {
                    let from = node_hashes[fastrand::usize(0..node_hashes.len())];
                    let to = node_hashes[fastrand::usize(0..node_hashes.len())];

                    let edge = GraphEdge::new(from, to, 3);
                    added_edges.insert((from, to));
                    fragment.add_edge(edge);
                }
            }

            // Check adjacency list
            let adj = fragment.get_adjacency_list();
            for (from_hash, neighbors) in &adj {
                for &to_hash in neighbors {
                    assert!(
                        added_edges.contains(&(*from_hash, to_hash)),
                        "Adjacency list should contain all added edges"
                    );
                }
            }

            // Check that all edges are represented
            let adj_edge_count: usize = adj.values().map(|v| v.len()).sum();
            assert_eq!(
                adj_edge_count,
                fragment.edges.len(),
                "Adjacency list should represent all edges"
            );
        }
    }

    #[test]
    fn property_tip_detection_correctness() {
        // Property: Tips should be nodes with in-degree=0 OR out-degree=0
        for _ in 0..15 {
            let mut fragment = GraphFragment::new(0);

            // Create a simple structure with known tips
            let sequences = vec!["AAAA", "TTTT", "GGGG", "CCCC"];
            let mut nodes = Vec::new();

            for seq in &sequences {
                if let Ok(kmer) = CanonicalKmer::new(seq) {
                    let node = GraphNode::new(kmer.clone(), 4);
                    nodes.push(kmer.hash);
                    fragment.add_node(node);
                }
            }

            // Create linear chain: 0->1->2->3 (0 and 3 are tips)
            for i in 0..nodes.len() - 1 {
                let edge = GraphEdge::new(nodes[i], nodes[i + 1], 3);
                fragment.add_edge(edge);
            }

            let tips = fragment.find_tips();
            let adjacency = fragment.get_adjacency_list();

            // Calculate in-degrees
            let mut in_degree: HashMap<u64, usize> = HashMap::new();
            for neighbors in adjacency.values() {
                for &neighbor in neighbors {
                    *in_degree.entry(neighbor).or_insert(0) += 1;
                }
            }

            // Verify tip detection
            for &node_hash in &nodes {
                let out_deg = adjacency.get(&node_hash).map_or(0, |v| v.len());
                let in_deg = in_degree.get(&node_hash).copied().unwrap_or(0);

                let is_tip = in_deg == 0 || out_deg == 0;
                let detected_as_tip = tips.contains(&node_hash);

                assert_eq!(is_tip, detected_as_tip,
                          "Node {node_hash:?} tip detection mismatch: in_deg={in_deg}, out_deg={out_deg}, is_tip={is_tip}, detected={detected_as_tip}");
            }
        }
    }
}

#[cfg(test)]
mod sequence_complexity_properties {
    use super::*;

    #[test]
    fn property_complexity_bounds() {
        // Property: Complexity should always be in [0, 1] range
        for _ in 0..100 {
            let length = fastrand::usize(1..50);
            let sequence = canonical_kmer_properties::generate_random_dna_sequence(length);

            let complexity = calculate_sequence_complexity(&sequence);
            assert!(
                (0.0..=1.0).contains(&complexity),
                "Complexity {complexity} should be in [0,1] for sequence {sequence}"
            );
        }
    }

    #[test]
    fn property_complexity_monotonicity() {
        // Property: More diverse sequences should have higher complexity

        // Homopolymer (lowest complexity)
        let homopolymer = "A".repeat(20);
        let homo_complexity = calculate_sequence_complexity(&homopolymer);

        // Two-base sequence
        let two_base = "AT".repeat(10);
        let two_complexity = calculate_sequence_complexity(&two_base);

        // Three-base sequence
        let three_base = "ATC".repeat(7);
        let three_complexity = calculate_sequence_complexity(&three_base);

        // Four-base sequence (maximum diversity)
        let four_base = "ATCG".repeat(5);
        let four_complexity = calculate_sequence_complexity(&four_base);

        assert!(
            homo_complexity <= two_complexity,
            "Two-base should have higher complexity than homopolymer"
        );
        assert!(
            two_complexity <= three_complexity,
            "Three-base should have higher complexity than two-base"
        );
        assert!(
            three_complexity <= four_complexity,
            "Four-base should have highest complexity"
        );
    }

    #[test]
    fn property_complexity_invariant_to_order() {
        // Property: Complexity should depend only on base frequencies, not order
        let bases = ['A', 'T', 'C', 'G'];

        for _ in 0..20 {
            // Create two sequences with same base composition but different order
            let mut counts = [0; 4];
            for _ in 0..fastrand::usize(10..30) {
                counts[fastrand::usize(0..4)] += 1;
            }

            let mut seq1 = String::new();
            let mut seq2 = String::new();

            for (i, &count) in counts.iter().enumerate() {
                seq1.extend(std::iter::repeat_n(bases[i], count));
                seq2.extend(std::iter::repeat_n(bases[i], count));
            }

            // Shuffle seq2 differently
            let mut seq2_chars: Vec<char> = seq2.chars().collect();
            for i in 0..seq2_chars.len() {
                let j = fastrand::usize(0..seq2_chars.len());
                seq2_chars.swap(i, j);
            }
            let seq2_shuffled: String = seq2_chars.into_iter().collect();

            let complexity1 = calculate_sequence_complexity(&seq1);
            let complexity2 = calculate_sequence_complexity(&seq2_shuffled);

            assert!(
                (complexity1 - complexity2).abs() < 1e-10,
                "Sequences with same base composition should have same complexity"
            );
        }
    }

    #[test]
    fn property_gc_content_bounds() {
        // Property: GC content should always be in [0, 1] range
        for _ in 0..100 {
            let length = fastrand::usize(1..100);
            let sequence = canonical_kmer_properties::generate_random_dna_sequence(length);

            let gc_content = calculate_gc_content(&sequence);
            assert!(
                (0.0..=1.0).contains(&gc_content),
                "GC content {gc_content} should be in [0,1] for sequence {sequence}"
            );
        }
    }

    #[test]
    fn property_gc_content_accuracy() {
        // Property: GC content should match manual calculation
        for _ in 0..50 {
            let sequence =
                canonical_kmer_properties::generate_random_dna_sequence(fastrand::usize(5..50));

            let calculated_gc = calculate_gc_content(&sequence);

            // Manual calculation
            let gc_count = sequence.chars().filter(|&c| c == 'G' || c == 'C').count();
            let expected_gc = gc_count as f64 / sequence.len() as f64;

            assert!(
                (calculated_gc - expected_gc).abs() < 1e-10,
                "GC content calculation should be accurate"
            );
        }
    }
}

#[cfg(test)]
mod assembly_invariants {
    use super::*;

    pub fn create_property_test_read(id: usize, sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "property_test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn property_assembly_chunk_invariants() {
        // Property: Assembly chunk should maintain consistent internal state
        for _ in 0..20 {
            let k_size = fastrand::usize(3..8);
            let mut chunk = AssemblyChunk::new(0, k_size);

            let read_count = fastrand::usize(1..10);
            let mut expected_reads = 0;

            for i in 0..read_count {
                let seq_len = fastrand::usize(k_size..30);
                let sequence = canonical_kmer_properties::generate_random_dna_sequence(seq_len);
                let read = create_property_test_read(i, &sequence);

                if chunk.add_read(read).is_ok() {
                    expected_reads += 1;
                }
            }

            chunk.finalize();

            // Invariants
            assert_eq!(chunk.processing_stats.reads_processed, expected_reads);
            assert_eq!(chunk.reads.len(), expected_reads);
            assert_eq!(
                chunk.graph_fragment.nodes.len(),
                chunk.processing_stats.nodes_created
            );
            assert_eq!(
                chunk.graph_fragment.edges.len(),
                chunk.processing_stats.edges_created
            );

            // Memory usage should be reasonable
            assert!(chunk.processing_stats.memory_usage_bytes > 0);

            // All nodes should have valid properties
            for node in chunk.graph_fragment.nodes.values() {
                assert!(node.coverage > 0);
                assert_eq!(node.kmer_size, k_size);
                assert!(!node.kmer.sequence.is_empty());
            }
        }
    }

    #[test]
    fn property_contig_generation_consistency() {
        // Property: Generated contigs should be consistent with graph structure
        for _ in 0..10 {
            let builder = AssemblyGraphBuilder::new(4, 6, 1);

            // Generate overlapping reads
            let base_sequences = ["ATCGATCGATCGATCG", "TCGATCGATCGATCGA", "CGATCGATCGATCGAT"];

            let mut reads = Vec::new();
            for (i, seq) in base_sequences.iter().enumerate() {
                reads.push(create_property_test_read(i, seq));
                // Add duplicate for coverage
                reads.push(create_property_test_read(i + 100, seq));
            }

            if let Ok(assembly) = builder.build(&reads) {
                // Invariants
                assert!(!assembly.contigs.is_empty(), "Should generate contigs");

                for contig in &assembly.contigs {
                    assert!(contig.id >= 0);
                    assert!(contig.length > 0);
                    assert_eq!(contig.length, contig.sequence.len());
                    assert!(contig.coverage > 0.0);
                    assert!(!contig.node_path.is_empty());

                    // All nodes in path should exist in graph
                    for &node_hash in &contig.node_path {
                        assert!(
                            assembly.graph_fragment.nodes.contains_key(&node_hash),
                            "Contig path node should exist in graph"
                        );
                    }
                }

                // Assembly statistics should be consistent
                let total_length: usize = assembly.contigs.iter().map(|c| c.length).sum();
                assert_eq!(assembly.assembly_stats.total_length, total_length);
                assert_eq!(assembly.assembly_stats.num_contigs, assembly.contigs.len());

                if !assembly.contigs.is_empty() {
                    let max_length = assembly.contigs.iter().map(|c| c.length).max().unwrap();
                    assert_eq!(assembly.assembly_stats.largest_contig, max_length);
                }
            }
        }
    }

    #[test]
    fn property_graph_merge_associativity() {
        // Property: Graph merging should be associative: (A ∪ B) ∪ C = A ∪ (B ∪ C)
        for _ in 0..10 {
            let mut frag_a = GraphFragment::new(0);
            let mut frag_b = GraphFragment::new(1);
            let mut frag_c = GraphFragment::new(2);

            // Add distinct nodes to each fragment
            let seqs_a = vec!["AAAA", "TTTT"];
            let seqs_b = vec!["GGGG", "CCCC"];
            let seqs_c = vec!["ATCG", "CGAT"];
            let common_seq = "ATAT"; // Common node for interesting merge behavior

            for seq in &seqs_a {
                if let Ok(kmer) = CanonicalKmer::new(seq) {
                    frag_a.add_node(GraphNode::new(kmer, 4));
                }
            }

            for seq in &seqs_b {
                if let Ok(kmer) = CanonicalKmer::new(seq) {
                    frag_b.add_node(GraphNode::new(kmer, 4));
                }
            }

            for seq in &seqs_c {
                if let Ok(kmer) = CanonicalKmer::new(seq) {
                    frag_c.add_node(GraphNode::new(kmer, 4));
                }
            }

            // Add common node to all fragments
            if let Ok(common_kmer) = CanonicalKmer::new(common_seq) {
                frag_a.add_node(GraphNode::new(common_kmer.clone(), 4));
                frag_b.add_node(GraphNode::new(common_kmer.clone(), 4));
                frag_c.add_node(GraphNode::new(common_kmer, 4));
            }

            // Test (A ∪ B) ∪ C
            let mut left_result = frag_a.clone();
            left_result.merge_with(frag_b.clone()).unwrap();
            left_result.merge_with(frag_c.clone()).unwrap();

            // Test A ∪ (B ∪ C)
            let mut temp_bc = frag_b.clone();
            temp_bc.merge_with(frag_c.clone()).unwrap();
            let mut right_result = frag_a.clone();
            right_result.merge_with(temp_bc).unwrap();

            // Results should be equivalent (same nodes and coverage)
            assert_eq!(
                left_result.nodes.len(),
                right_result.nodes.len(),
                "Merge results should have same node count"
            );

            for (hash, left_node) in &left_result.nodes {
                if let Some(right_node) = right_result.nodes.get(hash) {
                    assert_eq!(
                        left_node.coverage, right_node.coverage,
                        "Node coverage should be same in both merge orders"
                    );
                } else {
                    panic!("Node should exist in both merge results");
                }
            }
        }
    }
}
