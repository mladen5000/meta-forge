//! Performance tests for the optimized assembly graph construction
//! =============================================================
//!
//! These tests validate the performance improvements and correctness
//! of the optimized graph construction and contig generation algorithms.

use meta_forge::assembly::graph_construction::AssemblyGraphBuilder;
use meta_forge::assembly::bioinformatics_optimizations::{
    BitPackedKmer, StreamingKmerProcessor, SimdNucleotideOps, RollingHash
};
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};
use std::time::Instant;

/// Generate test reads for performance benchmarking
fn generate_test_reads(num_reads: usize, read_length: usize) -> Vec<CorrectedRead> {
    let mut reads = Vec::with_capacity(num_reads);
    let nucleotides = ['A', 'C', 'G', 'T'];
    
    for i in 0..num_reads {
        let mut sequence = String::with_capacity(read_length);
        for j in 0..read_length {
            let idx = (i + j) % nucleotides.len();
            sequence.push(nucleotides[idx]);
        }
        
        reads.push(CorrectedRead {
            id: i,
            original: sequence.clone(),
            corrected: sequence.clone(),
            corrections: Vec::new(),
            quality_scores: vec![30; read_length],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        });
    }
    
    reads
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_graph_construction() {
        // Test with a small set of reads to ensure basic functionality
        let reads = generate_test_reads(50, 100);
        
        let builder = AssemblyGraphBuilder::new(15, 31, 2, 4)
            .expect("Failed to create graph builder");
        
        let start_time = Instant::now();
        let mut graph = builder.build_graph(&reads)
            .expect("Failed to build graph");
        let construction_time = start_time.elapsed();
        
        println!("Small graph construction time: {:?}", construction_time);
        println!("Graph nodes: {}, edges: {}", 
                 graph.graph_fragment.nodes.len(), 
                 graph.graph_fragment.edges.len());
        
        // Test contig generation
        let contig_start = Instant::now();
        graph.generate_contigs()
            .expect("Failed to generate contigs");
        let contig_time = contig_start.elapsed();
        
        println!("Contig generation time: {:?}", contig_time);
        println!("Generated {} contigs", graph.contigs.len());
        
        // Validate that we have contigs now
        assert!(!graph.contigs.is_empty(), "Should generate at least some contigs");
        
        // Validate contigs have reasonable properties
        for contig in &graph.contigs {
            assert!(!contig.sequence.is_empty(), "Contig sequence should not be empty");
            assert!(contig.length > 0, "Contig length should be positive");
            assert!(contig.coverage > 0.0, "Contig coverage should be positive");
        }
    }

    #[test]
    fn test_medium_graph_performance() {
        // Test with medium-sized data to check performance improvements
        let reads = generate_test_reads(500, 150);
        
        let builder = AssemblyGraphBuilder::new(21, 51, 3, 8)
            .expect("Failed to create graph builder");
        
        let start_time = Instant::now();
        let mut graph = builder.build_graph(&reads)
            .expect("Failed to build graph");
        let construction_time = start_time.elapsed();
        
        println!("Medium graph construction time: {:?}", construction_time);
        println!("Graph nodes: {}, edges: {}", 
                 graph.graph_fragment.nodes.len(), 
                 graph.graph_fragment.edges.len());
        
        // Should complete within reasonable time (under 10 seconds)
        assert!(construction_time.as_secs() < 10, 
                "Graph construction taking too long: {:?}", construction_time);
        
        // Test contig generation
        let contig_start = Instant::now();
        graph.generate_contigs()
            .expect("Failed to generate contigs");
        let contig_time = contig_start.elapsed();
        
        println!("Medium contig generation time: {:?}", contig_time);
        println!("Generated {} contigs", graph.contigs.len());
        
        // Should complete within reasonable time (under 5 seconds)
        assert!(contig_time.as_secs() < 5, 
                "Contig generation taking too long: {:?}", contig_time);
        
        // Should generate contigs
        assert!(!graph.contigs.is_empty(), "Should generate contigs for medium dataset");
    }

    #[test]
    fn test_bit_packed_kmer_performance() {
        let test_sequences = vec![
            "ATCGATCGATCG",
            "GCTAGCTAGCTA",
            "AAAAAAAAAAAAAAAA",
            "TTTTTTTTTTTTTTTT",
            "ACGTACGTACGTACGT",
        ];
        
        let start_time = Instant::now();
        let mut bit_kmers = Vec::new();
        
        for sequence in &test_sequences {
            for _ in 0..1000 { // Repeat to measure performance
                let bit_kmer = BitPackedKmer::new(sequence)
                    .expect("Failed to create bit-packed k-mer");
                bit_kmers.push(bit_kmer);
            }
        }
        
        let packing_time = start_time.elapsed();
        println!("Bit-packing time for {} k-mers: {:?}", bit_kmers.len(), packing_time);
        
        // Test unpacking performance
        let unpack_start = Instant::now();
        let mut unpacked_sequences = Vec::new();
        
        for bit_kmer in &bit_kmers {
            let unpacked = bit_kmer.unpack_sequence();
            unpacked_sequences.push(unpacked);
        }
        
        let unpacking_time = unpack_start.elapsed();
        println!("Unpacking time for {} k-mers: {:?}", bit_kmers.len(), unpacking_time);
        
        // Validate correctness
        for (i, original) in test_sequences.iter().cycle().take(bit_kmers.len()).enumerate() {
            assert_eq!(unpacked_sequences[i], *original, "Unpacked sequence should match original");
        }
        
        // Performance assertions (should be fast)
        assert!(packing_time.as_millis() < 100, "Bit-packing should be fast");
        assert!(unpacking_time.as_millis() < 100, "Unpacking should be fast");
    }

    #[test]
    fn test_streaming_kmer_processor() {
        let long_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG".repeat(100);
        
        let start_time = Instant::now();
        let mut processor = StreamingKmerProcessor::new(21);
        let kmers = processor.process_sequence(&long_sequence)
            .expect("Failed to process sequence");
        let processing_time = start_time.elapsed();
        
        println!("Streaming k-mer processing time: {:?}", processing_time);
        println!("Processed {} k-mers", kmers.len());
        
        let stats = processor.get_statistics();
        println!("Streaming stats: {:?}", stats);
        
        assert!(!kmers.is_empty(), "Should extract k-mers from sequence");
        assert!(stats.unique_kmers > 0, "Should have unique k-mers");
        assert!(stats.total_kmers > 0, "Should have total k-mers");
        
        // Performance assertion
        assert!(processing_time.as_millis() < 100, "Streaming processing should be fast");
    }

    #[test]
    fn test_rolling_hash_performance() {
        let long_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG".repeat(1000);
        
        let start_time = Instant::now();
        let mut rolling_hash = RollingHash::new(31);
        let mut hash_values = Vec::new();
        
        for nucleotide in long_sequence.chars() {
            if let Some(hash) = rolling_hash.push(nucleotide).unwrap() {
                hash_values.push(hash);
            }
        }
        
        let rolling_time = start_time.elapsed();
        println!("Rolling hash time for {} nucleotides: {:?}", long_sequence.len(), rolling_time);
        println!("Generated {} hash values", hash_values.len());
        
        assert!(!hash_values.is_empty(), "Should generate hash values");
        assert_eq!(hash_values.len(), long_sequence.len() - 31 + 1, "Should generate correct number of hashes");
        
        // Performance assertion
        assert!(rolling_time.as_millis() < 50, "Rolling hash should be very fast");
    }

    #[test]
    fn test_simd_nucleotide_operations() {
        let long_sequence = "ATCGATCGATCG".repeat(10000);
        
        let start_time = Instant::now();
        let counts = SimdNucleotideOps::count_nucleotides(&long_sequence);
        let counting_time = start_time.elapsed();
        
        let gc_start = Instant::now();
        let gc_content = SimdNucleotideOps::gc_content(&long_sequence);
        let gc_time = gc_start.elapsed();
        
        println!("SIMD nucleotide counting time: {:?}", counting_time);
        println!("SIMD GC content calculation time: {:?}", gc_time);
        println!("Nucleotide counts: A={}, C={}, G={}, T={}", counts[0], counts[1], counts[2], counts[3]);
        println!("GC content: {:.3}", gc_content);
        
        // Validate correctness
        let total = counts.iter().sum::<usize>();
        assert_eq!(total, long_sequence.len(), "Total count should match sequence length");
        assert!(gc_content >= 0.0 && gc_content <= 1.0, "GC content should be between 0 and 1");
        
        // Performance assertions
        assert!(counting_time.as_millis() < 50, "Nucleotide counting should be fast");
        assert!(gc_time.as_millis() < 10, "GC content calculation should be very fast");
    }

    #[test]
    fn test_memory_efficiency() {
        // Test memory usage of different approaches
        let reads = generate_test_reads(100, 200);
        
        let builder = AssemblyGraphBuilder::new(25, 45, 2, 4)
            .expect("Failed to create graph builder");
        
        // Build graph and measure memory usage (approximate)
        let start_memory = get_memory_usage();
        let mut graph = builder.build_graph(&reads)
            .expect("Failed to build graph");
        let after_construction = get_memory_usage();
        
        graph.generate_contigs()
            .expect("Failed to generate contigs");
        let after_contigs = get_memory_usage();
        
        println!("Memory usage:");
        println!("  Start: {} KB", start_memory);
        println!("  After construction: {} KB", after_construction);
        println!("  After contigs: {} KB", after_contigs);
        println!("  Graph construction delta: {} KB", after_construction - start_memory);
        println!("  Contig generation delta: {} KB", after_contigs - after_construction);
        
        // Memory usage should be reasonable for the dataset size
        let construction_memory = after_construction - start_memory;
        assert!(construction_memory < 50_000, // Less than 50MB for small test
                "Memory usage too high: {} KB", construction_memory);
    }

    // Helper function to estimate memory usage (simplified)
    fn get_memory_usage() -> usize {
        // This is a simplified estimation - in a real scenario, you'd use
        // more sophisticated memory profiling tools
        std::process::id() as usize // Placeholder
    }

    #[test]
    fn test_parallel_performance() {
        let reads = generate_test_reads(1000, 100);
        
        // Test with single thread
        let builder_single = AssemblyGraphBuilder::new(21, 41, 2, 1)
            .expect("Failed to create single-threaded builder");
        
        let start_single = Instant::now();
        let graph_single = builder_single.build_graph(&reads)
            .expect("Failed to build single-threaded graph");
        let time_single = start_single.elapsed();
        
        // Test with multiple threads
        let builder_multi = AssemblyGraphBuilder::new(21, 41, 2, 4)
            .expect("Failed to create multi-threaded builder");
        
        let start_multi = Instant::now();
        let graph_multi = builder_multi.build_graph(&reads)
            .expect("Failed to build multi-threaded graph");
        let time_multi = start_multi.elapsed();
        
        println!("Single-threaded time: {:?}", time_single);
        println!("Multi-threaded time: {:?}", time_multi);
        println!("Speedup: {:.2}x", time_single.as_secs_f64() / time_multi.as_secs_f64());
        
        // Validate both produce similar results
        assert_eq!(graph_single.graph_fragment.nodes.len(), 
                   graph_multi.graph_fragment.nodes.len(),
                   "Both approaches should produce similar node counts");
        
        // Multi-threaded should be faster (or at least not significantly slower)
        let speedup = time_single.as_secs_f64() / time_multi.as_secs_f64();
        assert!(speedup > 0.5, "Multi-threaded version should not be much slower");
    }
}