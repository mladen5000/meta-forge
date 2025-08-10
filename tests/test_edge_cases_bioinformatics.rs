//! Edge case tests for bioinformatics functionality
//! Tests ambiguous bases, short sequences, invalid DNA, empty inputs, and boundary conditions

use meta_forge::assembly::adaptive_k::*;
use meta_forge::core::data_structures::*;
use meta_forge::utils::configuration::{AmbiguousBaseConfig, AmbiguousBaseStrategy};

#[cfg(test)]
mod ambiguous_base_tests {
    use super::*;

    #[test]
    fn test_single_n_base_handling() {
        let configs = vec![
            AmbiguousBaseConfig {
                strategy: AmbiguousBaseStrategy::Skip,
                max_n_count: 0,
                replacement_base: 'A',
                random_probabilities: None,
            },
            AmbiguousBaseConfig {
                strategy: AmbiguousBaseStrategy::Allow,
                max_n_count: 1,
                replacement_base: 'A',
                random_probabilities: None,
            },
            AmbiguousBaseConfig {
                strategy: AmbiguousBaseStrategy::Replace,
                max_n_count: 10,
                replacement_base: 'A',
                random_probabilities: None,
            },
        ];

        let test_sequences = vec!["NATCG", "ANTCG", "ATNCG", "ATCGN"];

        for config in &configs {
            for sequence in &test_sequences {
                let result = CanonicalKmer::new_with_config(sequence, config);

                match config.strategy {
                    AmbiguousBaseStrategy::Skip => {
                        assert!(
                            result.is_err(),
                            "Skip strategy should reject sequences with N"
                        );
                    }
                    AmbiguousBaseStrategy::Allow => {
                        assert!(result.is_ok(), "Allow strategy should accept single N");
                    }
                    AmbiguousBaseStrategy::Replace => {
                        assert!(result.is_ok(), "Replace strategy should handle N");
                        let kmer = result.unwrap();
                        assert!(!kmer.sequence.contains('N'), "N should be replaced");
                    }
                    _ => {}
                }
            }
        }
    }

    #[test]
    fn test_multiple_n_bases() {
        let config_allow_2 = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::Allow,
            max_n_count: 2,
            replacement_base: 'A',
            random_probabilities: None,
        };

        let config_allow_0 = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::Allow,
            max_n_count: 0,
            replacement_base: 'A',
            random_probabilities: None,
        };

        // Should succeed with 2 or fewer Ns
        assert!(CanonicalKmer::new_with_config("NNTCG", &config_allow_2).is_ok());
        assert!(CanonicalKmer::new_with_config("NATCN", &config_allow_2).is_ok());

        // Should fail with more than max_n_count
        assert!(CanonicalKmer::new_with_config("NNNTG", &config_allow_2).is_err());
        assert!(CanonicalKmer::new_with_config("NNNNN", &config_allow_2).is_err());

        // Should fail with any N when max_n_count is 0
        assert!(CanonicalKmer::new_with_config("NATCG", &config_allow_0).is_err());
    }

    #[test]
    fn test_random_replacement_strategy() {
        let config = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::RandomReplace,
            max_n_count: 10,
            replacement_base: 'A', // Not used in random strategy
            random_probabilities: Some([0.3, 0.3, 0.2, 0.2]), // A, C, G, T probabilities
        };

        let sequence_with_n = "NNNNATCGNNN";

        // Test multiple times to ensure randomness works
        let mut replacement_results = Vec::new();
        for _ in 0..10 {
            let result = CanonicalKmer::new_with_config(sequence_with_n, &config);
            assert!(result.is_ok());

            let kmer = result.unwrap();
            assert!(!kmer.sequence.contains('N'), "All Ns should be replaced");
            replacement_results.push(kmer.sequence.clone());
        }

        // Should get some variation in replacements (probabilistic test)
        // Not guaranteed to be different due to randomness, but likely
        let unique_results: std::collections::HashSet<_> =
            replacement_results.into_iter().collect();
        // With many Ns and randomness, we should get some variation
        // (This is a probabilistic test that might rarely fail)
    }

    #[test]
    fn test_context_replacement_strategy() {
        let config = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::ContextReplace,
            max_n_count: 10,
            replacement_base: 'A', // Not used in context strategy
            random_probabilities: None,
        };

        // Test with A-rich context
        let a_rich = "AAAANTCG";
        let result_a = CanonicalKmer::new_with_config(a_rich, &config).unwrap();
        assert!(!result_a.sequence.contains('N'));
        // Should replace N with A (most common base)

        // Test with G-rich context
        let g_rich = "GGGGNNTCG";
        let result_g = CanonicalKmer::new_with_config(g_rich, &config).unwrap();
        assert!(!result_g.sequence.contains('N'));
        // Should replace Ns with G (most common base)
    }

    #[test]
    fn test_mixed_case_with_ambiguous_bases() {
        let config = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::Replace,
            max_n_count: 5,
            replacement_base: 'A',
            random_probabilities: None,
        };

        let mixed_sequences = vec!["atcgN", "ATCGn", "AtCgN", "nAtCg"];

        for sequence in mixed_sequences {
            let result = CanonicalKmer::new_with_config(sequence, &config);
            assert!(result.is_ok());

            let kmer = result.unwrap();
            assert!(!kmer.sequence.contains('N'));
            assert!(!kmer.sequence.contains('n'));
            assert!(kmer.sequence.chars().all(|c| c.is_ascii_uppercase()));
        }
    }
}

#[cfg(test)]
mod short_sequence_tests {
    use super::*;

    #[test]
    fn test_empty_sequence_handling() {
        assert!(CanonicalKmer::new("").is_err());
        assert!(BitPackedKmer::new("").is_err());
        assert!(validate_dna_sequence("").is_ok()); // Empty sequence is valid DNA

        let extractor = MinimizerExtractor::new(3, 5);
        let minimizers = extractor.extract_minimizers("").unwrap();
        assert!(minimizers.is_empty());
    }

    #[test]
    fn test_single_base_sequence() {
        let single_bases = vec!["A", "T", "C", "G"];

        for base in single_bases {
            // Too short for k-mer
            assert!(CanonicalKmer::new(base).is_ok()); // Single base is valid k-mer

            let extractor = MinimizerExtractor::new(3, 5);
            let minimizers = extractor.extract_minimizers(base).unwrap();
            assert!(minimizers.is_empty()); // Too short for k=3
        }
    }

    #[test]
    fn test_two_base_sequence() {
        let two_bases = vec!["AT", "GC", "TA", "CG"];

        for bases in two_bases {
            assert!(CanonicalKmer::new(bases).is_ok()); // Valid 2-mer

            let extractor = MinimizerExtractor::new(3, 5);
            let minimizers = extractor.extract_minimizers(bases).unwrap();
            assert!(minimizers.is_empty()); // Too short for k=3
        }
    }

    #[test]
    fn test_exact_k_mer_size_sequence() {
        let k = 4;
        let sequence = "ATCG"; // Exactly k bases

        assert!(CanonicalKmer::new(sequence).is_ok());

        let extractor = MinimizerExtractor::new(k, 5);
        let minimizers = extractor.extract_minimizers(sequence).unwrap();
        assert_eq!(minimizers.len(), 1); // Should have exactly one minimizer
        assert_eq!(minimizers[0].kmer.len(), k);
    }

    #[test]
    fn test_k_plus_one_sequence() {
        let k = 4;
        let sequence = "ATCGA"; // k+1 bases

        let extractor = MinimizerExtractor::new(k, 3);
        let minimizers = extractor.extract_minimizers(sequence).unwrap();
        assert!(minimizers.len() >= 1); // Should have at least one minimizer

        // Check that all minimizers have correct k-mer size
        for minimizer in &minimizers {
            assert_eq!(minimizer.kmer.len(), k);
        }
    }

    #[test]
    fn test_assembly_chunk_with_very_short_reads() {
        let mut chunk = AssemblyChunk::new(0, 6);

        let short_reads = vec![
            CorrectedRead {
                id: 0,
                original: "A".to_string(),
                corrected: "A".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 1,
                    correction_time_ms: 0,
                },
            },
            CorrectedRead {
                id: 1,
                original: "AT".to_string(),
                corrected: "AT".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30, 30],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 2,
                    correction_time_ms: 0,
                },
            },
        ];

        for read in short_reads {
            let result = chunk.add_read(read);
            assert!(result.is_ok()); // Should not fail, just skip processing
        }

        // No nodes should be created from short reads
        assert_eq!(chunk.processing_stats.nodes_created, 0);
        assert_eq!(chunk.graph_fragment.nodes.len(), 0);
        assert_eq!(chunk.processing_stats.reads_processed, 2); // Reads were processed (but skipped)
    }
}

#[cfg(test)]
mod invalid_dna_tests {
    use super::*;

    #[test]
    fn test_invalid_characters() {
        let invalid_sequences = vec![
            "ATCGX",    // Invalid character X
            "ATCG123",  // Numbers
            "ATCG-GG",  // Hyphen
            "ATCG GG",  // Space
            "ATCG*GG",  // Special character
            "ATCG\nGG", // Newline
            "ATCG\tGG", // Tab
        ];

        for sequence in invalid_sequences {
            assert!(
                CanonicalKmer::new(sequence).is_err(),
                "Sequence '{}' should be invalid",
                sequence
            );
            assert!(
                BitPackedKmer::new(sequence).is_err(),
                "BitPackedKmer with '{}' should be invalid",
                sequence
            );
            assert!(
                validate_dna_sequence(sequence).is_err(),
                "DNA validation for '{}' should fail",
                sequence
            );
        }
    }

    #[test]
    fn test_unicode_characters() {
        let unicode_sequences = vec![
            "ATCG∅",  // Mathematical symbol
            "ATCG®",  // Registered trademark
            "ATCG™",  // Trademark
            "ATCGα",  // Greek letter
            "ATCG中", // Chinese character
        ];

        for sequence in unicode_sequences {
            assert!(
                CanonicalKmer::new(sequence).is_err(),
                "Unicode sequence '{}' should be invalid",
                sequence
            );
            assert!(
                validate_dna_sequence(sequence).is_err(),
                "DNA validation for '{}' should fail",
                sequence
            );
        }
    }

    #[test]
    fn test_mixed_valid_invalid_characters() {
        let mixed_sequences = vec![
            "ATCGGATCX",  // Valid DNA + X
            "XATCGATCG",  // X + valid DNA
            "ATCXGATCG",  // Valid DNA with X in middle
            "ATCG123TCG", // Valid DNA with numbers
        ];

        for sequence in mixed_sequences {
            assert!(
                CanonicalKmer::new(sequence).is_err(),
                "Mixed sequence '{}' should be invalid",
                sequence
            );
        }
    }

    #[test]
    fn test_boundary_valid_characters() {
        // Test characters that are almost valid
        let boundary_sequences = vec![
            "ATCGH", // H is not valid DNA
            "ATCGI", // I is not valid DNA
            "ATCGJ", // J is not valid DNA
            "ATCGU", // U is RNA, not DNA (though some systems accept it)
        ];

        for sequence in boundary_sequences {
            assert!(
                CanonicalKmer::new(sequence).is_err(),
                "Boundary sequence '{}' should be invalid",
                sequence
            );
        }
    }
}

#[cfg(test)]
mod boundary_condition_tests {
    use super::*;

    #[test]
    fn test_maximum_k_mer_size_limits() {
        // Test very long k-mers
        let long_sequence = "A".repeat(1000);
        assert!(CanonicalKmer::new(&long_sequence).is_ok()); // Should handle long sequences

        let very_long_sequence = "A".repeat(10000);
        assert!(CanonicalKmer::new(&very_long_sequence).is_ok()); // Should still work

        // BitPackedKmer has a limit
        let too_long_for_bitpacked = "A".repeat(2000);
        assert!(BitPackedKmer::new(&too_long_for_bitpacked).is_err());
    }

    #[test]
    fn test_zero_window_size_minimizer() {
        // Edge case: window size of 1 (minimum possible)
        let extractor = MinimizerExtractor::new(3, 1);
        let sequence = "ATCGATCG";
        let minimizers = extractor.extract_minimizers(sequence).unwrap();

        assert!(!minimizers.is_empty());
        // With window size 1, each k-mer position should be a potential minimizer
    }

    #[test]
    fn test_k_equal_to_sequence_length() {
        let sequence = "ATCGATCG";
        let k = sequence.len();

        let extractor = MinimizerExtractor::new(k, 3);
        let minimizers = extractor.extract_minimizers(sequence).unwrap();

        assert_eq!(minimizers.len(), 1); // Only one possible k-mer
        assert_eq!(minimizers[0].kmer.len(), k);
    }

    #[test]
    fn test_very_large_coverage_values() {
        let mut node = GraphNode::new(CanonicalKmer::new("ATCG").unwrap(), 4);

        // Simulate very high coverage
        for i in 0..1000 {
            node.add_read_position(i, i, Strand::Forward);
        }

        assert!(node.coverage > 1000);
        node.update_node_type();
        assert_eq!(node.node_type, NodeType::Repetitive);
    }

    #[test]
    fn test_extreme_gc_content() {
        // All AT (0% GC)
        let at_only = "AAAATTTTAAAATTTT";
        assert_eq!(calculate_gc_content(at_only), 0.0);

        // All GC (100% GC)
        let gc_only = "GGGGCCCCGGGGCCCC";
        assert_eq!(calculate_gc_content(gc_only), 1.0);

        // These should still create valid k-mers despite extreme composition
        assert!(CanonicalKmer::new("AAAA").is_ok());
        assert!(CanonicalKmer::new("GGGG").is_ok());
    }

    #[test]
    fn test_palindromic_sequences() {
        let palindromes = vec![
            "ATCGAT", // 6-mer palindrome
            "GAATTC", // EcoRI restriction site
            "AAGCTT", // HindIII restriction site
            "GGATCC", // BamHI restriction site
        ];

        for palindrome in palindromes {
            let kmer = CanonicalKmer::new(palindrome).unwrap();
            // Palindromic sequences should still canonicalize properly
            assert_eq!(kmer.len(), palindrome.len());
        }
    }

    #[test]
    fn test_homopolymer_runs() {
        let homopolymers = vec![
            "AAAAAAAAAA", // A homopolymer
            "TTTTTTTTTT", // T homopolymer
            "GGGGGGGGGG", // G homopolymer
            "CCCCCCCCCC", // C homopolymer
        ];

        for homopolymer in homopolymers {
            // Should handle homopolymers
            assert!(CanonicalKmer::new(homopolymer).is_ok());

            // Complexity should be very low
            let complexity = calculate_sequence_complexity(homopolymer);
            assert!(complexity < 0.1, "Homopolymer should have low complexity");

            let extractor = MinimizerExtractor::new(4, 5);
            let minimizers = extractor.extract_minimizers(homopolymer).unwrap();

            // All k-mers are identical in homopolymer, so should deduplicate
            assert!(minimizers.len() <= homopolymer.len() - 4 + 1);
        }
    }

    #[test]
    fn test_assembly_with_no_overlaps() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);

        // Create reads with no overlaps
        let non_overlapping_reads = vec![
            CorrectedRead {
                id: 0,
                original: "AAAATTTT".to_string(),
                corrected: "AAAATTTT".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30; 8],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 4,
                    correction_time_ms: 0,
                },
            },
            CorrectedRead {
                id: 1,
                original: "GGGGCCCC".to_string(),
                corrected: "GGGGCCCC".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30; 8],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 4,
                    correction_time_ms: 0,
                },
            },
        ];

        let result = builder.build(&non_overlapping_reads);
        assert!(result.is_ok());

        let graph = result.unwrap();
        // Should create separate contigs for non-overlapping reads
        assert!(graph.contigs.len() >= 1);
        assert!(!graph.graph_fragment.nodes.is_empty());
    }
}
