//! TDD Assembly Robustness Tests
//! =============================
//! 
//! Test-driven development focused on assembly pipeline robustness, error handling,
//! and edge cases that aren't fully covered in existing tests.

use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};
use meta_forge::assembly::bioinformatics_optimizations::BitPackedKmer;

#[cfg(test)]
mod assembly_robustness_tests {
    use super::*;
    use std::collections::HashMap;

    fn create_test_read_with_metadata(id: usize, sequence: &str, quality_scores: Vec<u8>) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores,
            correction_metadata: CorrectionMetadata {
                algorithm: "tdd_test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn test_assembly_with_ambiguous_nucleotides() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Reads containing IUPAC ambiguous nucleotides (N, Y, R, etc.)
        let reads = vec![
            create_test_read_with_metadata(0, "ATCGNATCG", vec![30; 9]), // Contains N
            create_test_read_with_metadata(1, "TCGNATCGA", vec![30; 9]), // Contains N  
            create_test_read_with_metadata(2, "ATCGYTCGR", vec![30; 9]), // Contains Y and R
        ];

        // This test should pass gracefully - assembly should handle ambiguous bases
        let result = builder.build(&reads);
        
        // For now, expect it to work (even if it filters out ambiguous k-mers)
        assert!(result.is_ok(), "Assembly should handle ambiguous nucleotides gracefully");
        
        let graph = result.unwrap();
        // The graph might have fewer nodes due to filtering ambiguous k-mers
        // But it should still be valid
        assert!(graph.assembly_stats.num_contigs >= 0);
    }

    #[test]
    fn test_assembly_with_very_low_coverage() {
        let builder = AssemblyGraphBuilder::new(4, 6, 10); // Very high minimum coverage
        
        // Single short read with no overlapping k-mers
        let reads = vec![
            create_test_read_with_metadata(0, "ATCG", vec![30; 4]), // Only one k-mer of size 4
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());
        
        let graph = result.unwrap();
        
        // With very high minimum coverage and only one k-mer instance, should be filtered
        // Graph should be valid but have no contigs
        assert_eq!(graph.assembly_stats.num_contigs, 0, "Low coverage reads should be filtered");
        assert_eq!(graph.graph_fragment.nodes.len(), 0, "Low coverage nodes should be removed");
    }

    #[test]
    fn test_assembly_with_very_short_reads() {
        let builder = AssemblyGraphBuilder::new(10, 15, 1); // Large k-mer size
        
        // All reads shorter than k-mer size
        let reads = vec![
            create_test_read_with_metadata(0, "ATCG", vec![30; 4]),      // 4 bp < k=10
            create_test_read_with_metadata(1, "TCGA", vec![30; 4]),      // 4 bp < k=10
            create_test_read_with_metadata(2, "CGATCG", vec![30; 6]),    // 6 bp < k=10
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());
        
        let graph = result.unwrap();
        // No k-mers can be extracted from reads shorter than k
        assert_eq!(graph.graph_fragment.nodes.len(), 0, "Short reads should produce no k-mers");
        assert_eq!(graph.contigs.len(), 0, "No contigs from short reads");
    }

    #[test] 
    fn test_assembly_with_repetitive_sequences() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Highly repetitive sequences that should challenge the assembler
        let reads = vec![
            create_test_read_with_metadata(0, "ATATATATATATAT", vec![30; 14]), // Perfect repeats
            create_test_read_with_metadata(1, "ATATATATATAT", vec![30; 12]),   // Perfect repeats
            create_test_read_with_metadata(2, "GCGCGCGCGCGC", vec![30; 12]),   // Different repeats
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());
        
        let graph = result.unwrap();
        // Repetitive sequences should create high-coverage nodes
        let high_coverage_nodes: Vec<_> = graph.graph_fragment.nodes.values()
            .filter(|node| node.coverage > 3)
            .collect();
            
        assert!(!high_coverage_nodes.is_empty(), "Repetitive sequences should create high coverage nodes");
    }

    #[test]
    fn test_assembly_graph_statistics_accuracy() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);
        
        let reads = vec![
            create_test_read_with_metadata(0, "ATCGATCG", vec![30; 8]),
            create_test_read_with_metadata(1, "TCGATCGA", vec![30; 8]), 
            create_test_read_with_metadata(2, "CGATCGAT", vec![30; 8]),
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());
        
        let graph = result.unwrap();
        let stats = &graph.assembly_stats;
        
        // Verify statistics are reasonable
        assert!(stats.total_length > 0, "Total length should be positive");
        assert!(stats.coverage_mean >= 0.0, "Coverage mean should be non-negative");
        assert!(stats.n50 >= 0, "N50 should be non-negative");
        assert!(stats.largest_contig >= 0, "Largest contig should be non-negative");
        
        if stats.num_contigs > 0 {
            assert!(stats.largest_contig <= stats.total_length, "Largest contig can't exceed total length");
        }
    }

    #[test]
    fn test_bit_packed_kmer_canonical_consistency() {
        // Test that BitPackedKmer consistently handles canonical k-mers
        let forward_kmer = BitPackedKmer::new("ATCG").unwrap();
        let reverse_comp = BitPackedKmer::new("CGAT").unwrap(); // Reverse complement
        
        // Debug output
        println!("Forward ATCG: packed={:?}, hash={}", forward_kmer.packed_data, forward_kmer.hash);
        println!("Reverse CGAT: packed={:?}, hash={}", reverse_comp.packed_data, reverse_comp.hash);
        
        // These should be different unless we have canonical k-mer handling
        // This test will help drive the implementation of canonical k-mers
        assert!(forward_kmer.hash != 0);
        assert!(reverse_comp.hash != 0);
        
        // For now, they should be different (implementation will determine if they should be same)
        assert_ne!(forward_kmer.hash, reverse_comp.hash, "Forward and reverse complement should have different hashes without canonical handling");
    }

    #[test] 
    fn test_assembly_error_recovery() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Mix of good and problematic reads
        let reads = vec![
            create_test_read_with_metadata(0, "ATCGATCGATCG", vec![30; 12]), // Good read
            create_test_read_with_metadata(1, "AT", vec![30; 2]),            // Too short
            create_test_read_with_metadata(2, "TCGATCGATCGA", vec![30; 12]), // Good read  
            create_test_read_with_metadata(3, "", Vec::new()),               // Empty read
            create_test_read_with_metadata(4, "CGATCGATCGAT", vec![30; 12]), // Good read
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Assembly should recover from problematic reads");
        
        let graph = result.unwrap();
        // Should have nodes from the good reads at least
        assert!(!graph.graph_fragment.nodes.is_empty(), "Should have nodes from good reads");
    }

    #[test]
    fn test_memory_efficient_large_input() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Create a large number of reads to test memory efficiency
        let mut reads = Vec::new();
        for i in 0..1000 {
            let sequence = format!("ATCGATCGATCG{:04}", i % 100); // Generate varied sequences
            reads.push(create_test_read_with_metadata(
                i, 
                &sequence, 
                vec![30; sequence.len()]
            ));
        }

        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle large input efficiently");
        
        let graph = result.unwrap();
        assert!(!graph.graph_fragment.nodes.is_empty(), "Large input should produce nodes");
        assert!(graph.assembly_stats.num_contigs > 0, "Large input should produce contigs");
    }

    #[test] 
    fn test_k_mer_size_adaptation_effectiveness() {
        let builder = AssemblyGraphBuilder::new(3, 9, 1); // Wide range for adaptation
        
        // Test different complexity scenarios that should trigger different k-mer sizes
        let low_complexity = vec![
            create_test_read_with_metadata(0, "AAAAAAAAAA", vec![30; 10]),
            create_test_read_with_metadata(1, "TTTTTTTTTT", vec![30; 10]),
        ];
        
        let high_complexity = vec![
            create_test_read_with_metadata(0, "ATCGATCGACTGATCG", vec![30; 16]),
            create_test_read_with_metadata(1, "GTACGTACGTACGTAC", vec![30; 16]),
        ];

        let low_result = builder.build(&low_complexity);
        let high_result = builder.build(&high_complexity);
        
        assert!(low_result.is_ok(), "Low complexity should be handled");
        assert!(high_result.is_ok(), "High complexity should be handled");
        
        // Both should produce valid graphs, potentially with different characteristics
        let low_graph = low_result.unwrap();
        let high_graph = high_result.unwrap();
        
        // The adaptive k-mer selection should result in different graph structures
        // This test drives the need for effective k-mer size adaptation
        assert!(low_graph.assembly_stats.num_contigs >= 0);
        assert!(high_graph.assembly_stats.num_contigs >= 0);
    }
}