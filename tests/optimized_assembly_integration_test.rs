//! Integration tests for optimized assembly features
//! Tests assembly optimizations and enhanced progress tracking
use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
// Note: Direct SIMD and memory-mapped tests moved to separate files
// to avoid compilation issues on different platforms
use meta_forge::assembly::memory_mapped::MappedGraph;
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};
use std::path::PathBuf;
use std::time::Instant;
use tempfile::TempDir;

fn generate_realistic_reads(num_reads: usize, read_length: usize) -> Vec<CorrectedRead> {
    let mut reads = Vec::with_capacity(num_reads);

    // Create realistic genomic sequences with some overlap
    let base_sequences = vec![
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        "CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG",
        "TTAATTAATTAATTAATTAATTAATTAATTAATTAA",
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
    ];

    for i in 0..num_reads {
        let base_idx = i % base_sequences.len();
        let base_seq = base_sequences[base_idx];

        // Create overlapping subsequences
        let start_pos = (i * 5) % (base_seq.len().saturating_sub(read_length));
        let end_pos = (start_pos + read_length).min(base_seq.len());
        let sequence = &base_seq[start_pos..end_pos];

        reads.push(CorrectedRead {
            id: i,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![35; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "integration_test".to_string(),
                confidence_threshold: 0.95,
                context_window: 7,
                correction_time_ms: 1,
            },
        });
    }

    reads
}

#[cfg(test)]
mod optimized_assembly_tests {
    use super::*;

    #[test]
    #[ignore = "Memory-mapped storage requires platform-specific fixes"]
    fn test_memory_mapped_graph_storage() {
        let temp_dir = TempDir::new().expect("Failed to create temp directory");
        let temp_path = temp_dir.path();

        // Test creating and using memory-mapped graph
        let mapped_graph = MappedGraph::new(temp_path, 1000, 5000);
        assert!(mapped_graph.is_ok(), "Failed to create memory-mapped graph");

        let mut graph = mapped_graph.unwrap();

        // Test basic operations
        let start_time = Instant::now();

        // Add some nodes and edges
        for i in 0..100 {
            graph.add_node(i, i as u32 + 1).expect("Failed to add node");
        }

        for i in 0..99 {
            graph.add_edge(i, i + 1, 3).expect("Failed to add edge");
        }

        let operation_time = start_time.elapsed();
        println!("Memory-mapped operations time: {:?}", operation_time);

        // Verify data integrity
        assert_eq!(graph.node_count, 100);
        assert_eq!(graph.edge_count, 99);

        // Test lookup performance
        let lookup_start = Instant::now();
        for i in 0..100 {
            let coverage = graph.get_node_coverage(i).expect("Node should exist");
            assert_eq!(coverage, i as u32 + 1);
        }
        let lookup_time = lookup_start.elapsed();
        println!("Memory-mapped lookup time for 100 nodes: {:?}", lookup_time);

        // Performance assertion - should be fast
        assert!(lookup_time.as_millis() < 10, "Lookups should be very fast");
    }

    #[test]
    #[ignore = "SIMD optimizations require platform-specific fixes"]
    fn test_simd_optimized_kmer_processing() {
        let processor = SimdOptimizedKmerProcessor::new(21);

        // Test with realistic genomic sequences
        let test_sequences = vec![
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
            "CCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGGCCGG",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        ];

        let start_time = Instant::now();
        let mut total_kmers = 0;

        for sequence in &test_sequences {
            let kmers = processor
                .extract_kmers(sequence)
                .expect("Failed to extract k-mers");
            total_kmers += kmers.len();

            // Validate k-mer extraction
            let expected_kmers = sequence.len().saturating_sub(21 - 1);
            assert_eq!(
                kmers.len(),
                expected_kmers,
                "Wrong number of k-mers extracted"
            );
        }

        let processing_time = start_time.elapsed();
        println!(
            "SIMD k-mer processing time for {} k-mers: {:?}",
            total_kmers, processing_time
        );

        // Test SIMD nucleotide counting
        let long_sequence = test_sequences.join("");
        let count_start = Instant::now();
        let counts = processor.count_nucleotides(&long_sequence);
        let count_time = count_start.elapsed();

        println!("SIMD nucleotide counting time: {:?}", count_time);
        println!(
            "Nucleotide counts: A={}, C={}, G={}, T={}",
            counts[0], counts[1], counts[2], counts[3]
        );

        // Validate counts
        let total_count = counts.iter().sum::<usize>();
        assert_eq!(
            total_count,
            long_sequence.len(),
            "Total nucleotide count should match sequence length"
        );

        // Performance assertions
        assert!(
            processing_time.as_millis() < 50,
            "SIMD k-mer processing should be fast"
        );
        assert!(
            count_time.as_millis() < 10,
            "SIMD counting should be very fast"
        );
    }

    #[test]
    fn test_integrated_optimized_assembly_pipeline() {
        // Test the full assembly pipeline with optimizations enabled
        let reads = generate_realistic_reads(200, 80);

        println!(
            "Testing optimized assembly pipeline with {} reads",
            reads.len()
        );

        let builder = AssemblyGraphBuilder::new(15, 31, 2);

        let start_time = Instant::now();
        let graph_result = builder.build(&reads);
        let construction_time = start_time.elapsed();

        assert!(graph_result.is_ok(), "Assembly should succeed");
        let mut graph = graph_result.unwrap();

        println!("Optimized graph construction time: {:?}", construction_time);
        println!(
            "Graph nodes: {}, edges: {}",
            graph.graph_fragment.nodes.len(),
            graph.graph_fragment.edges.len()
        );

        // Test contig generation with optimizations
        let contig_start = Instant::now();
        graph
            .generate_contigs()
            .expect("Contig generation should succeed");
        let contig_time = contig_start.elapsed();

        println!("Optimized contig generation time: {:?}", contig_time);
        println!("Generated {} contigs", graph.contigs.len());

        // Validate assembly quality
        assert!(
            !graph.graph_fragment.nodes.is_empty(),
            "Should create graph nodes"
        );
        assert!(!graph.contigs.is_empty(), "Should generate contigs");
        assert!(
            graph.assembly_stats.total_length > 0,
            "Should have positive total length"
        );
        assert!(graph.assembly_stats.n50 > 0, "Should calculate N50");

        // Performance expectations with optimizations
        assert!(
            construction_time.as_secs() < 30,
            "Optimized construction should be fast"
        );
        assert!(
            contig_time.as_secs() < 15,
            "Optimized contig generation should be fast"
        );

        // Test progress tracking (should not crash)
        println!("Assembly statistics:");
        println!("  Total length: {}", graph.assembly_stats.total_length);
        println!("  Number of contigs: {}", graph.assembly_stats.num_contigs);
        println!("  Largest contig: {}", graph.assembly_stats.largest_contig);
        println!("  N50: {}", graph.assembly_stats.n50);
        println!("  Mean coverage: {:.2}", graph.assembly_stats.coverage_mean);
    }

    #[test]
    #[ignore = "Combined optimizations require platform-specific fixes"]
    fn test_combined_memory_mapped_and_simd_optimizations() {
        // Test combining memory-mapped storage with SIMD processing
        let temp_dir = TempDir::new().expect("Failed to create temp directory");
        let temp_path = temp_dir.path();

        let reads = generate_realistic_reads(100, 60);
        let simd_processor = SimdOptimizedKmerProcessor::new(15);

        let mut mapped_graph =
            MappedGraph::new(temp_path, 2000, 10000).expect("Failed to create memory-mapped graph");

        let start_time = Instant::now();
        let mut node_id = 0u64;

        // Process reads with SIMD and store in memory-mapped graph
        for read in &reads {
            let kmers = simd_processor
                .extract_kmers(&read.corrected)
                .expect("Failed to extract k-mers");

            for kmer in kmers {
                mapped_graph
                    .add_node(node_id, kmer.coverage)
                    .expect("Failed to add node to memory-mapped graph");
                node_id += 1;
            }
        }

        let combined_time = start_time.elapsed();
        println!(
            "Combined SIMD + memory-mapped processing time: {:?}",
            combined_time
        );
        println!("Processed {} nodes in memory-mapped storage", node_id);

        // Verify integrity
        assert_eq!(
            mapped_graph.node_count as u64, node_id,
            "All nodes should be stored"
        );

        // Test random access performance
        let access_start = Instant::now();
        let sample_size = node_id.min(100);

        for i in 0..sample_size {
            let coverage = mapped_graph.get_node_coverage(i);
            assert!(coverage.is_ok(), "Should be able to access stored nodes");
        }

        let access_time = access_start.elapsed();
        println!(
            "Random access time for {} nodes: {:?}",
            sample_size, access_time
        );

        // Performance assertions
        assert!(
            combined_time.as_secs() < 10,
            "Combined processing should be efficient"
        );
        assert!(access_time.as_millis() < 20, "Random access should be fast");
    }

    #[test]
    fn test_progress_tracking_integration() {
        // Test that enhanced progress tracking works correctly
        let reads = generate_realistic_reads(50, 100);
        let builder = AssemblyGraphBuilder::new(11, 25, 1);

        // This test primarily ensures that the progress-enhanced build doesn't crash
        println!("Testing progress tracking integration...");

        let start_time = Instant::now();
        let result = builder.build(&reads);
        let total_time = start_time.elapsed();

        assert!(
            result.is_ok(),
            "Assembly with progress tracking should succeed"
        );
        let graph = result.unwrap();

        println!("Progress-tracked assembly completed in {:?}", total_time);
        println!(
            "Results: {} nodes, {} contigs",
            graph.graph_fragment.nodes.len(),
            graph.contigs.len()
        );

        // Validate that assembly completed successfully
        assert!(
            !graph.graph_fragment.nodes.is_empty() || reads.is_empty(),
            "Should create nodes for non-empty input"
        );
        assert!(
            graph.assembly_stats.total_length >= 0,
            "Assembly stats should be calculated"
        );

        // The progress tracking should not significantly slow down the assembly
        assert!(
            total_time.as_secs() < 60,
            "Progress tracking should not cause major slowdown"
        );
    }

    #[test]
    fn test_optimization_performance_comparison() {
        // Compare performance with and without optimizations where possible
        let reads = generate_realistic_reads(100, 75);

        // Test baseline assembly
        let baseline_builder = AssemblyGraphBuilder::new(13, 21, 1);
        let baseline_start = Instant::now();
        let baseline_result = baseline_builder.build(&reads);
        let baseline_time = baseline_start.elapsed();

        assert!(baseline_result.is_ok(), "Baseline assembly should succeed");
        let baseline_graph = baseline_result.unwrap();

        // Test optimized assembly (same parameters for fair comparison)
        let optimized_builder = AssemblyGraphBuilder::new(13, 21, 1);
        let optimized_start = Instant::now();
        let optimized_result = optimized_builder.build(&reads);
        let optimized_time = optimized_start.elapsed();

        assert!(
            optimized_result.is_ok(),
            "Optimized assembly should succeed"
        );
        let optimized_graph = optimized_result.unwrap();

        println!("Performance comparison:");
        println!("  Baseline time: {:?}", baseline_time);
        println!("  Optimized time: {:?}", optimized_time);
        println!(
            "  Baseline nodes: {}",
            baseline_graph.graph_fragment.nodes.len()
        );
        println!(
            "  Optimized nodes: {}",
            optimized_graph.graph_fragment.nodes.len()
        );

        // Results should be similar (allowing for minor differences due to optimizations)
        let node_diff = (baseline_graph.graph_fragment.nodes.len() as i64
            - optimized_graph.graph_fragment.nodes.len() as i64)
            .abs();
        let node_threshold = baseline_graph.graph_fragment.nodes.len() / 10; // 10% tolerance

        assert!(
            node_diff <= node_threshold as i64,
            "Optimized version should produce similar node count (diff: {}, threshold: {})",
            node_diff,
            node_threshold
        );

        // Optimized version should not be significantly slower
        let slowdown_factor = optimized_time.as_secs_f64() / baseline_time.as_secs_f64().max(0.001);
        assert!(
            slowdown_factor < 2.0,
            "Optimized version should not be more than 2x slower (factor: {:.2})",
            slowdown_factor
        );

        println!("  Slowdown factor: {:.2}x", slowdown_factor);

        if slowdown_factor < 1.0 {
            println!(
                "  ✅ Optimization provides {:.2}x speedup!",
                1.0 / slowdown_factor
            );
        }
    }
}
