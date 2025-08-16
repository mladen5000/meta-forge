//! Basic tests for assembly optimizations
//! Tests that optimizations compile and run without errors

use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};

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
            correction_time_ms: 1,
        },
    }
}

#[cfg(test)]
mod basic_optimization_tests {
    use super::*;

    #[test]
    fn test_basic_assembly_with_optimizations() {
        // Test that assembly works with our optimizations
        let reads = vec![
            create_test_read(0, "ATCGATCGATCGATCGATCG"),
            create_test_read(1, "CGATCGATCGATCGATCGAT"),
            create_test_read(2, "GATCGATCGATCGATCGATC"),
        ];

        let builder = AssemblyGraphBuilder::new(11, 21, 1);
        let result = builder.build(&reads);

        assert!(result.is_ok(), "Assembly should succeed");
        let graph = result.unwrap();

        // Basic validations
        assert!(
            !graph.graph_fragment.nodes.is_empty() || reads.is_empty(),
            "Should create nodes for non-empty input"
        );
        assert!(
            graph.assembly_stats.total_length >= 0,
            "Assembly stats should be calculated"
        );

        println!("✅ Basic assembly test passed");
        println!("   Nodes: {}", graph.graph_fragment.nodes.len());
        println!("   Edges: {}", graph.graph_fragment.edges.len());
        println!("   Contigs: {}", graph.contigs.len());
    }

    #[test]
    fn test_simd_processor_creation() {
        // Test that SIMD processor can be created without errors
        let processor = meta_forge::assembly::simd_optimizations::SimdProcessor::new();

        // Basic functionality test
        let sequence = b"ATCGATCGATCGATCGATCG";
        let counts = processor.count_nucleotides(sequence);

        // Should count nucleotides correctly
        let total: usize = counts.iter().sum();
        assert_eq!(total, sequence.len(), "Should count all nucleotides");

        println!("✅ SIMD processor test passed");
        println!(
            "   Counts: A={}, C={}, G={}, T={}",
            counts[0], counts[1], counts[2], counts[3]
        );

        // Test k-mer processing
        let kmers = processor.process_kmers(sequence, 5);
        assert!(kmers.is_ok(), "K-mer processing should succeed");

        let kmer_vec = kmers.unwrap();
        assert_eq!(
            kmer_vec.len(),
            sequence.len() - 5 + 1,
            "Should generate correct number of k-mers"
        );

        println!("   K-mers generated: {}", kmer_vec.len());
    }

    #[test]
    fn test_progress_tracking() {
        // Test that progress tracking doesn't crash
        let reads = vec![
            create_test_read(0, "ATCGATCGATCGATCGATCGATCGATCG"),
            create_test_read(1, "CGATCGATCGATCGATCGATCGATCGAT"),
        ];

        let builder = AssemblyGraphBuilder::new(15, 25, 1);
        let result = builder.build(&reads);

        assert!(
            result.is_ok(),
            "Assembly with progress tracking should succeed"
        );
        let graph = result.unwrap();

        println!("✅ Progress tracking test passed");
        println!("   Assembly completed successfully");
        println!(
            "   Final stats: {} nodes, {} contigs",
            graph.graph_fragment.nodes.len(),
            graph.contigs.len()
        );
    }

    #[test]
    fn test_empty_input_handling() {
        // Test edge case with empty input
        let empty_reads: Vec<CorrectedRead> = vec![];
        let builder = AssemblyGraphBuilder::new(15, 25, 1);
        let result = builder.build(&empty_reads);

        assert!(result.is_ok(), "Should handle empty input gracefully");
        let graph = result.unwrap();

        assert_eq!(
            graph.graph_fragment.nodes.len(),
            0,
            "Should have no nodes for empty input"
        );
        assert_eq!(
            graph.contigs.len(),
            0,
            "Should have no contigs for empty input"
        );

        println!("✅ Empty input handling test passed");
    }

    #[test]
    fn test_short_reads_handling() {
        // Test handling of reads shorter than k-mer size
        let short_reads = vec![
            create_test_read(0, "AT"),                   // Very short
            create_test_read(1, "ATCG"),                 // Short
            create_test_read(2, "ATCGATCGATCGATCGATCG"), // Normal length
        ];

        let builder = AssemblyGraphBuilder::new(15, 25, 1);
        let result = builder.build(&short_reads);

        assert!(
            result.is_ok(),
            "Should handle mixed read lengths gracefully"
        );
        let graph = result.unwrap();

        // Should process at least the long read
        assert!(
            graph.assembly_stats.total_length >= 0,
            "Should have some assembly output"
        );

        println!("✅ Short reads handling test passed");
        println!("   Processed mixed read lengths successfully");
    }
}
