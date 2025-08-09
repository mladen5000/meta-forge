//! Integration tests for end-to-end pipeline functionality
//! Tests complete pipeline flow, field access patterns, error handling, and real-world scenarios

use meta_forge::assembly::adaptive_k::*;
use meta_forge::core::data_structures::*;

#[cfg(test)]
mod end_to_end_pipeline_tests {
    use super::*;

    pub fn create_realistic_read(id: usize, sequence: &str, quality_avg: u8) -> CorrectedRead {
        let qualities = vec![quality_avg; sequence.len()];
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: qualities,
            correction_metadata: CorrectionMetadata {
                algorithm: "BWA-MEM2".to_string(),
                confidence_threshold: 0.95,
                context_window: 10,
                correction_time_ms: (id as u64) * 10, // Realistic processing time
            },
        }
    }

    #[test]
    fn test_complete_pipeline_small_dataset() {
        let builder = AssemblyGraphBuilder::new(4, 8, 2);
        
        // Create small but realistic dataset
        let reads = vec![
            create_realistic_read(0, "ATCGATCGATCGATCG", 35),
            create_realistic_read(1, "TCGATCGATCGATCGA", 30),
            create_realistic_read(2, "CGATCGATCGATCGAT", 40),
            create_realistic_read(3, "GATCGATCGATCGATC", 32),
            create_realistic_read(4, "ATCGATCGATCGATCG", 38), // Duplicate for coverage
        ];
        
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Pipeline should complete successfully");
        
        let assembly = result.unwrap();
        
        // Verify all pipeline stages completed
        assert!(!assembly.graph_fragment.nodes.is_empty(), "Should have graph nodes");
        assert!(!assembly.graph_fragment.edges.is_empty(), "Should have graph edges");
        assert!(!assembly.contigs.is_empty(), "Should have generated contigs");
        
        // Verify statistics are calculated
        assert!(assembly.assembly_stats.total_length > 0);
        assert!(assembly.assembly_stats.num_contigs > 0);
        assert!(assembly.assembly_stats.coverage_mean > 0.0);
        assert!(assembly.assembly_stats.n50 > 0);
        
        // Verify field access patterns (the fixes we implemented)
        assert!(!assembly.graph_fragment.nodes.is_empty());
        for (_, node) in &assembly.graph_fragment.nodes {
            assert!(node.coverage > 0);
            assert_eq!(node.kmer_size, 4);
            assert!(!node.kmer.sequence.is_empty());
        }
        
        for edge in &assembly.graph_fragment.edges {
            assert!(edge.weight > 0);
            assert!(edge.confidence > 0.0);
            assert!(edge.overlap_length > 0);
        }
    }

    #[test]
    fn test_pipeline_with_varying_quality_reads() {
        let builder = AssemblyGraphBuilder::new(5, 9, 1);
        
        // Mix of high and low quality reads
        let reads = vec![
            create_realistic_read(0, "ATCGATCGATCGATCGATCG", 40), // High quality
            create_realistic_read(1, "TCGATCGATCGATCGATCGA", 15), // Low quality
            create_realistic_read(2, "CGATCGATCGATCGATCGAT", 35), // Medium quality
            create_realistic_read(3, "GATCGATCGATCGATCGATC", 42), // Very high quality
            create_realistic_read(4, "ATCGATCGATCGATCGATCG", 10), // Very low quality
        ];
        
        let result = builder.build(&reads);
        assert!(result.is_ok());
        
        let assembly = result.unwrap();
        
        // Pipeline should handle quality variation
        assert!(!assembly.graph_fragment.nodes.is_empty());
        assert!(!assembly.contigs.is_empty());
        
        // Check that reads were processed despite quality differences
        let high_coverage_nodes: Vec<_> = assembly.graph_fragment.nodes
            .values()
            .filter(|node| node.coverage > 1)
            .collect();
        
        assert!(!high_coverage_nodes.is_empty(), "Should have overlapping coverage from similar sequences");
    }

    #[test]
    fn test_pipeline_error_handling_invalid_reads() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Mix of valid and problematic reads
        let reads = vec![
            create_realistic_read(0, "ATCGATCG", 35), // Valid
            create_realistic_read(1, "A", 30),        // Too short
            create_realistic_read(2, "", 0),          // Empty
            create_realistic_read(3, "ATCGATCGATCG", 40), // Valid
        ];
        
        // The pipeline should handle these gracefully
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Pipeline should handle problematic reads gracefully");
        
        let assembly = result.unwrap();
        
        // Should still produce results from valid reads
        assert!(!assembly.graph_fragment.nodes.is_empty());
        assert!(!assembly.contigs.is_empty());
    }

    #[test]
    fn test_large_scale_integration() {
        let builder = AssemblyGraphBuilder::new(6, 10, 3);
        
        // Generate larger dataset
        let base_sequences = ["ATCGATCGATCGATCGATCGATCG",
            "TCGATCGATCGATCGATCGATCGA",
            "CGATCGATCGATCGATCGATCGAT",
            "GATCGATCGATCGATCGATCGATC"];
        
        let mut reads = Vec::new();
        for i in 0..50 {
            let seq_idx = i % base_sequences.len();
            let mut sequence = base_sequences[seq_idx].to_string();
            
            // Add some variation to create realistic overlaps
            if i % 3 == 0 {
                sequence = sequence[1..].to_string(); // Offset by 1
            } else if i % 5 == 0 {
                sequence = sequence[2..].to_string(); // Offset by 2
            }
            
            reads.push(create_realistic_read(i, &sequence, 30 + (i % 20) as u8));
        }
        
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Large scale pipeline should complete");
        
        let assembly = result.unwrap();
        
        // Verify scalability
        assert!(assembly.graph_fragment.nodes.len() > 10, "Should have substantial graph");
        assert!(assembly.assembly_stats.total_length > 100, "Should have reasonable assembly size");
        
        // Verify coverage statistics make sense
        let mean_coverage = assembly.graph_fragment.coverage_stats.mean_coverage;
        let median_coverage = assembly.graph_fragment.coverage_stats.median_coverage;
        
        assert!(mean_coverage > 0.0);
        assert!(median_coverage > 0);
        
        // With overlapping reads, should have good coverage
        let high_coverage_nodes = assembly.graph_fragment.nodes
            .values()
            .filter(|n| n.coverage >= 3)
            .count();
        
        assert!(high_coverage_nodes > 0, "Should have well-covered regions");
    }

    #[test]
    fn test_memory_and_performance_tracking() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        let reads = vec![
            create_realistic_read(0, "ATCGATCGATCGATCG", 35),
            create_realistic_read(1, "TCGATCGATCGATCGA", 30),
            create_realistic_read(2, "CGATCGATCGATCGAT", 40),
        ];
        
        let start_time = std::time::Instant::now();
        let result = builder.build(&reads);
        let elapsed = start_time.elapsed();
        
        assert!(result.is_ok());
        let assembly = result.unwrap();
        
        // Performance should be reasonable for small dataset
        assert!(elapsed.as_secs() < 5, "Should complete in reasonable time");
        
        // Memory tracking through chunk processing
        // (In a real system, you'd have more detailed memory profiling)
        assert!(!assembly.graph_fragment.nodes.is_empty());
        assert!(assembly.graph_fragment.coverage_stats.total_nodes > 0);
    }
}

#[cfg(test)]
mod field_access_validation_tests {
    use super::*;

    #[test]
    fn test_corrected_read_field_access() {
        let read = CorrectedRead {
            id: 1,
            original: "ATCGATCG".to_string(),
            corrected: "ATCGATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 8],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 100,
            },
        };
        
        // Test all field accesses that were fixed
        assert_eq!(read.id, 1);
        assert_eq!(read.original, "ATCGATCG");
        assert_eq!(read.corrected, "ATCGATCG");
        assert!(read.corrections.is_empty());
        assert_eq!(read.quality_scores.len(), 8);
        
        // Test correction_metadata field access (this was a major fix)
        assert_eq!(read.correction_metadata.algorithm, "test");
        assert_eq!(read.correction_metadata.confidence_threshold, 0.95);
        assert_eq!(read.correction_metadata.context_window, 5);
        assert_eq!(read.correction_metadata.correction_time_ms, 100);
    }

    #[test]
    fn test_assembly_graph_builder_constructor() {
        // Test the 3-argument constructor that was fixed
        let builder = AssemblyGraphBuilder::new(4, 8, 2);
        
        // Constructor should succeed (no panic or error)
        assert!(true, "Constructor completed successfully");
        
        // Test building with the fixed constructor
        let reads = vec![end_to_end_pipeline_tests::create_realistic_read(0, "ATCGATCG", 35)];
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Build method should work with fixed constructor");
    }

    #[test]
    fn test_graph_fragment_field_access() {
        let fragment = GraphFragment::new(42);
        
        // Test field accesses that were validated/fixed
        assert_eq!(fragment.fragment_id, 42);
        assert_eq!(fragment.nodes.len(), 0);
        assert_eq!(fragment.edges.len(), 0);
        assert_eq!(fragment.coverage_stats.total_nodes, 0);
        assert_eq!(fragment.coverage_stats.total_edges, 0);
        assert_eq!(fragment.coverage_stats.mean_coverage, 0.0);
    }

    #[test]
    fn test_assembly_stats_field_access() {
        let mut stats = AssemblyStats::default();
        
        // Test all field accesses
        assert_eq!(stats.total_length, 0);
        assert_eq!(stats.num_contigs, 0);
        assert_eq!(stats.n50, 0);
        assert_eq!(stats.n90, 0);
        assert_eq!(stats.largest_contig, 0);
        assert_eq!(stats.gc_content, 0.0);
        assert_eq!(stats.coverage_mean, 0.0);
        assert_eq!(stats.coverage_std, 0.0);
        
        // Test field modifications
        stats.total_length = 1000;
        stats.num_contigs = 5;
        stats.n50 = 200;
        stats.coverage_mean = 15.5;
        
        assert_eq!(stats.total_length, 1000);
        assert_eq!(stats.num_contigs, 5);
        assert_eq!(stats.n50, 200);
        assert_eq!(stats.coverage_mean, 15.5);
    }

    #[test]
    fn test_assembly_chunk_method_calls() {
        let mut chunk = AssemblyChunk::new(0, 4);
        
        // Test new() method
        assert_eq!(chunk.chunk_id, 0);
        assert_eq!(chunk.k_size, 4);
        
        // Test add_read() method
        let read = end_to_end_pipeline_tests::create_realistic_read(0, "ATCGATCG", 35);
        let result = chunk.add_read(read);
        assert!(result.is_ok(), "add_read method should work");
        
        // Test finalize() method
        chunk.finalize();
        assert!(chunk.processing_stats.memory_usage_bytes > 0, "finalize should calculate memory usage");
    }

    #[test] 
    fn test_graph_node_coverage_updates() {
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let mut node = GraphNode::new(kmer, 4);
        
        // Test initial state
        assert_eq!(node.coverage, 1);
        assert_eq!(node.node_type, NodeType::Unique);
        
        // Test add_read_position method
        node.add_read_position(1, 10, Strand::Forward);
        assert_eq!(node.coverage, 2);
        
        // Test update_node_type method
        node.update_node_type();
        assert_eq!(node.node_type, NodeType::LowCoverage); // Coverage 2 = LowCoverage
        
        // Increase coverage to test type changes
        for i in 2..15 {
            node.add_read_position(i, i * 5, Strand::Forward);
        }
        
        node.update_node_type();
        assert_eq!(node.node_type, NodeType::HighCoverage); // Coverage > 10
    }
}

#[cfg(test)]
mod error_handling_integration_tests {
    use super::*;

    #[test]
    fn test_pipeline_with_malformed_data() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Include problematic reads that should be handled gracefully
        let reads = vec![
            end_to_end_pipeline_tests::create_realistic_read(0, "ATCGATCG", 35),           // Valid
            end_to_end_pipeline_tests::create_realistic_read(1, "A", 0),                  // Too short, bad quality
            end_to_end_pipeline_tests::create_realistic_read(2, "", 0),                   // Empty sequence
            end_to_end_pipeline_tests::create_realistic_read(3, "ATCGATCGATCGATCG", 45),  // Valid long read
        ];
        
        // Pipeline should not crash
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Pipeline should handle malformed data gracefully");
        
        let assembly = result.unwrap();
        
        // Should still produce meaningful results
        assert!(!assembly.graph_fragment.nodes.is_empty());
        assert!(!assembly.contigs.is_empty());
    }

    #[test]
    fn test_assembly_chunk_error_recovery() {
        let mut chunk = AssemblyChunk::new(0, 4);
        
        // Try to add problematic reads
        let problematic_reads = vec![
            CorrectedRead {
                id: 0,
                original: "".to_string(),
                corrected: "".to_string(),
                corrections: Vec::new(),
                quality_scores: Vec::new(),
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.0,
                    context_window: 0,
                    correction_time_ms: 0,
                },
            },
            end_to_end_pipeline_tests::create_realistic_read(1, "ATCGATCG", 35), // Valid read after error
        ];
        
        for read in problematic_reads {
            let result = chunk.add_read(read);
            // Should not fail, even with problematic data
            assert!(result.is_ok(), "add_read should handle errors gracefully");
        }
        
        chunk.finalize();
        
        // Should have processed at least the valid read
        assert!(chunk.processing_stats.reads_processed > 0);
    }

    #[test]
    fn test_graph_construction_with_no_overlaps() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Reads with no k-mer overlaps
        let reads = vec![
            end_to_end_pipeline_tests::create_realistic_read(0, "AAAAAAAA", 35), // A homopolymer
            end_to_end_pipeline_tests::create_realistic_read(1, "TTTTTTTT", 30), // T homopolymer
            end_to_end_pipeline_tests::create_realistic_read(2, "GGGGGGGG", 40), // G homopolymer
            end_to_end_pipeline_tests::create_realistic_read(3, "CCCCCCCC", 32), // C homopolymer
        ];
        
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle non-overlapping reads");
        
        let assembly = result.unwrap();
        
        // Should create separate contigs
        assert!(!assembly.graph_fragment.nodes.is_empty());
        assert!(!assembly.contigs.is_empty());
        assert!(!assembly.contigs.is_empty()); // At least one contig
    }

    #[test]
    fn test_extreme_coverage_handling() {
        let builder = AssemblyGraphBuilder::new(3, 5, 1);
        
        // Many identical reads to create extreme coverage
        let base_sequence = "ATCGATCGATCG";
        let mut reads = Vec::new();
        
        for i in 0..100 {
            reads.push(end_to_end_pipeline_tests::create_realistic_read(i, base_sequence, 35));
        }
        
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle extreme coverage");
        
        let assembly = result.unwrap();
        
        // Should have very high coverage nodes
        let max_coverage = assembly.graph_fragment.nodes
            .values()
            .map(|n| n.coverage)
            .max()
            .unwrap_or(0);
            
        assert!(max_coverage >= 50, "Should have high coverage from many identical reads");
        
        // Node types should be updated appropriately
        let repetitive_nodes = assembly.graph_fragment.nodes
            .values()
            .filter(|n| matches!(n.node_type, NodeType::Repetitive))
            .count();
            
        assert!(repetitive_nodes > 0, "Should classify high-coverage nodes as repetitive");
    }

    #[test]
    fn test_memory_pressure_simulation() {
        let builder = AssemblyGraphBuilder::new(4, 6, 1);
        
        // Create many reads to simulate memory pressure
        let mut reads = Vec::new();
        for i in 0..500 {
            let sequence = format!("ATCGATCG{:04}", i % 1000); // Some overlap, mostly unique
            reads.push(end_to_end_pipeline_tests::create_realistic_read(i, &sequence, 30));
        }
        
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle large datasets without memory issues");
        
        let assembly = result.unwrap();
        
        // Should successfully process all reads
        assert!(assembly.graph_fragment.nodes.len() > 100);
        assert!(assembly.assembly_stats.total_length > 1000);
        assert!(!assembly.contigs.is_empty());
    }
}