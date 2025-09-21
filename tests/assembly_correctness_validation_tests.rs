//! Assembly Correctness Validation Tests
//! =====================================
//!
//! Comprehensive correctness validation for the assembly implementation.
//! Focuses on edge cases, boundary conditions, error handling, and consistency.
//! Validates that optimizations don't break correctness.

use meta_forge::assembly::laptop_assembly::{
    LaptopAssembler, LaptopConfig, LaptopAssemblyGraph, CompactKmer, BoundedKmerCounter
};
use meta_forge::core::data_structures::{CorrectedRead, Contig, ContigType, CorrectionMetadata};
use anyhow::Result;
use std::time::Duration;
use std::collections::HashSet;

/// Test compact k-mer encoding correctness
#[test]
fn test_compact_kmer_encoding_correctness() {
    println!("🧬 Testing CompactKmer encoding correctness...");

    // Test basic encoding/decoding
    let test_sequences = vec![
        "A", "T", "C", "G",
        "AT", "CG", "TA", "GC",
        "ATCG", "GCTA", "AAAA", "TTTT",
        "ATCGATCGATCGATCGATCGATCGATCG", // Max length (32)
    ];

    for seq in &test_sequences {
        let kmer = CompactKmer::new(seq).unwrap();
        let decoded = kmer.to_string();
        assert_eq!(seq, &decoded, "Encoding/decoding mismatch for: {}", seq);

        // Test hash consistency
        let kmer2 = CompactKmer::new(seq).unwrap();
        assert_eq!(kmer.hash(), kmer2.hash(), "Hash inconsistency for: {}", seq);
    }

    println!("✅ CompactKmer encoding correctness validated");
}

#[test]
fn test_compact_kmer_boundary_conditions() {
    println!("🔍 Testing CompactKmer boundary conditions...");

    // Test empty sequence
    assert!(CompactKmer::new("").is_err(), "Empty sequence should fail");

    // Test invalid characters
    let invalid_sequences = vec!["X", "ATCGX", "N", "ATCN"];
    for seq in &invalid_sequences {
        assert!(CompactKmer::new(seq).is_err(), "Invalid sequence should fail: {}", seq);
    }

    // Test maximum length (32 nucleotides - fits in u64 with 2 bits per base)
    let max_seq = "A".repeat(32);
    assert!(CompactKmer::new(&max_seq).is_ok(), "Max length sequence should work");

    // Test over maximum length
    let over_max = "A".repeat(33);
    assert!(CompactKmer::new(&over_max).is_err(), "Over-max length should fail");

    // Test case insensitivity
    let lower_case = CompactKmer::new("atcg").unwrap();
    let upper_case = CompactKmer::new("ATCG").unwrap();
    assert_eq!(lower_case.to_string(), upper_case.to_string(), "Case should be normalized");

    println!("✅ CompactKmer boundary conditions validated");
}

#[test]
fn test_bounded_kmer_counter_memory_management() {
    println!("💾 Testing BoundedKmerCounter memory management...");

    // Test with very small memory budget
    let mut counter = BoundedKmerCounter::new(1); // 1MB limit

    // Add many k-mers to test memory bounding
    let test_kmers: Vec<u64> = (0..100000).collect();
    for kmer_hash in &test_kmers {
        counter.add_kmer(*kmer_hash);
    }

    let (unique, total, dropped, memory) = counter.get_stats();

    // Verify memory constraints are respected
    assert!(memory <= 2 * 1024 * 1024, "Memory usage exceeds budget: {} bytes", memory); // Allow some overhead
    assert_eq!(total, 100000, "Total k-mers seen should be accurate");
    assert!(memory <= 10 * 1024 * 1024, "Memory usage should be bounded: {} bytes", memory); // Use memory as proxy for size check

    // Test that frequent k-mers are prioritized
    counter.add_kmer(12345); // Add same k-mer multiple times
    counter.add_kmer(12345);
    counter.add_kmer(12345);

    let frequent = counter.get_frequent_kmers(3);
    assert!(frequent.iter().any(|(hash, _)| *hash == 12345), "Frequent k-mer should be retained");

    println!("✅ BoundedKmerCounter memory management validated");
}

#[test]
fn test_assembly_graph_malformed_data() {
    println!("🚨 Testing assembly graph with malformed data...");

    let config = LaptopConfig::low_memory();
    let mut graph = LaptopAssemblyGraph::new(config);

    // Test with empty reads
    let empty_reads = vec![];
    assert!(graph.build_from_reads(&empty_reads, 21).is_ok(), "Empty reads should be handled gracefully");

    // Test with reads shorter than k
    let short_reads = vec![
        create_test_read(0, "AT"),      // Length 2, k=21
        create_test_read(1, "ATCGATCG"), // Length 8, k=21
    ];
    assert!(graph.build_from_reads(&short_reads, 21).is_ok(), "Short reads should be handled gracefully");

    // Test with reads containing only one type of nucleotide
    let homopolymer_reads = vec![
        create_test_read(0, "AAAAAAAAAAAAAAAAAAAAAAAAA"),
        create_test_read(1, "TTTTTTTTTTTTTTTTTTTTTTTTT"),
    ];
    assert!(graph.build_from_reads(&homopolymer_reads, 21).is_ok(), "Homopolymer reads should be handled");

    // Test with extremely repetitive reads
    let repetitive_reads = vec![
        create_test_read(0, "ATCGATCGATCGATCGATCGATCGATCGATCG"),
        create_test_read(1, "TCGATCGATCGATCGATCGATCGATCGATCGA"),
        create_test_read(2, "CGATCGATCGATCGATCGATCGATCGATCGAT"),
    ];
    assert!(graph.build_from_reads(&repetitive_reads, 21).is_ok(), "Repetitive reads should be handled");

    println!("✅ Assembly graph malformed data handling validated");
}

#[test]
fn test_assembly_timeout_and_error_handling() {
    println!("⏱️ Testing assembly timeout and error handling...");

    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);

    // Test with very short timeout
    let reads = generate_large_test_dataset(1000);
    let result = assembler.assemble_with_timeout(&reads, Duration::from_millis(1));

    // Should either succeed quickly or timeout gracefully
    match result {
        Ok(contigs) => {
            println!("✅ Assembly completed within 1ms: {} contigs", contigs.len());
        }
        Err(e) => {
            if e.to_string().contains("timeout") {
                println!("✅ Assembly timed out gracefully as expected");
            } else {
                panic!("Unexpected error during timeout test: {}", e);
            }
        }
    }

    // Test with reasonable timeout
    let result = assembler.assemble_with_timeout(&reads, Duration::from_secs(30));
    assert!(result.is_ok(), "Assembly should succeed with reasonable timeout");

    println!("✅ Timeout and error handling validated");
}

#[test]
fn test_assembly_output_consistency() {
    println!("🔄 Testing assembly output consistency across runs...");

    let config = LaptopConfig::medium_memory();
    let assembler = LaptopAssembler::new(config);

    // Create deterministic test data
    let reads = create_deterministic_test_reads();

    // Run assembly multiple times
    let mut results = Vec::new();
    for i in 0..3 {
        println!("  Run {}: Assembling...", i + 1);
        let contigs = assembler.assemble(&reads).unwrap();
        results.push(contigs);
    }

    // Compare results for consistency
    let first_result = &results[0];
    for (i, result) in results.iter().enumerate().skip(1) {
        println!("  Comparing run {} with run 1...", i + 1);

        // Results should be deterministic (same number of contigs)
        assert_eq!(
            first_result.len(),
            result.len(),
            "Inconsistent number of contigs between runs: {} vs {}",
            first_result.len(),
            result.len()
        );

        // Total assembled length should be consistent
        let first_total: usize = first_result.iter().map(|c| c.length).sum();
        let current_total: usize = result.iter().map(|c| c.length).sum();

        let length_diff = (first_total as i64 - current_total as i64).abs();
        let tolerance = (first_total as f64 * 0.05) as i64; // 5% tolerance

        assert!(
            length_diff <= tolerance,
            "Total assembly length inconsistent: {} vs {} (diff: {})",
            first_total, current_total, length_diff
        );
    }

    println!("✅ Assembly output consistency validated");
}

#[test]
fn test_memory_cleanup_and_emergency_recovery() {
    println!("🧹 Testing memory cleanup and emergency recovery...");

    let config = LaptopConfig::low_memory(); // Force low memory conditions
    let mut graph = LaptopAssemblyGraph::new(config);

    // Generate large dataset to trigger memory pressure
    let reads = generate_large_test_dataset(2000);

    // Build graph which should trigger emergency cleanup
    let result = graph.build_from_reads(&reads, 21);

    match result {
        Ok(()) => {
            println!("✅ Assembly completed despite memory pressure");

            // Test emergency cleanup function
            let memory_before = graph.memory_usage_mb();
            graph.emergency_cleanup().unwrap();
            let memory_after = graph.memory_usage_mb();

            assert!(
                memory_after <= memory_before,
                "Emergency cleanup should reduce or maintain memory usage"
            );

            // Graph should still be functional after cleanup
            let contigs = graph.generate_contigs().unwrap();
            println!("  Generated {} contigs after cleanup", contigs.len());
        }
        Err(e) => {
            // If assembly fails due to memory, check error message is informative
            let error_msg = e.to_string();
            assert!(
                error_msg.contains("memory") || error_msg.contains("timeout"),
                "Error should be informative about memory/timeout: {}",
                error_msg
            );
            println!("✅ Assembly gracefully handled memory constraints: {}", error_msg);
        }
    }

    println!("✅ Memory cleanup and emergency recovery validated");
}

#[test]
fn test_contig_generation_correctness() {
    println!("🧬 Testing contig generation correctness...");

    let config = LaptopConfig::medium_memory();
    let assembler = LaptopAssembler::new(config);

    // Test with known overlapping reads
    let overlapping_reads = vec![
        create_test_read(0, "ATCGATCGATCGATCGATCGATCG"),
        create_test_read(1, "TCGATCGATCGATCGATCGATCGA"),
        create_test_read(2, "CGATCGATCGATCGATCGATCGAT"),
        create_test_read(3, "GATCGATCGATCGATCGATCGATC"),
    ];

    let contigs = assembler.assemble(&overlapping_reads).unwrap();

    // Validate contig properties
    for contig in &contigs {
        // Basic sanity checks
        assert!(contig.length > 0, "Contig length should be positive");
        assert_eq!(contig.sequence.len(), contig.length, "Sequence length should match reported length");
        assert!(contig.coverage > 0.0, "Coverage should be positive");

        // Sequence should contain only valid DNA bases
        for c in contig.sequence.chars() {
            assert!(
                matches!(c, 'A' | 'T' | 'C' | 'G'),
                "Invalid character in contig sequence: {}",
                c
            );
        }

        // Check for reasonable coverage values
        assert!(
            contig.coverage <= 100.0,
            "Coverage seems unreasonably high: {}",
            contig.coverage
        );
    }

    // Test that contigs are non-empty and reasonable
    assert!(!contigs.is_empty(), "Should generate at least one contig from overlapping reads");

    let total_length: usize = contigs.iter().map(|c| c.length).sum();
    assert!(total_length >= 20, "Total assembly length should be reasonable");

    println!("✅ Contig generation correctness validated");
}

#[test]
fn test_concurrent_kmer_processing() {
    println!("⚡ Testing concurrent k-mer processing...");

    let config = LaptopConfig::high_memory(); // Enable parallel processing
    let mut graph = LaptopAssemblyGraph::new(config);

    // Create large dataset to trigger parallel processing
    let reads = generate_large_test_dataset(5000);

    let start_time = std::time::Instant::now();
    let result = graph.build_from_reads(&reads, 21);
    let elapsed = start_time.elapsed();

    assert!(result.is_ok(), "Concurrent processing should complete successfully");

    println!("  Processed {} reads in {:.2}s", reads.len(), elapsed.as_secs_f64());

    // Verify graph state is consistent after concurrent processing
    let contigs = graph.generate_contigs().unwrap();

    // Basic consistency checks
    assert!(!contigs.is_empty(), "Concurrent processing should generate contigs");

    let total_length: usize = contigs.iter().map(|c| c.length).sum();
    assert!(total_length > 0, "Total assembly length should be positive");

    // Check for data races by verifying coverage values are reasonable
    for contig in &contigs {
        assert!(
            contig.coverage > 0.0 && contig.coverage < 10000.0,
            "Coverage value suggests potential data race: {}",
            contig.coverage
        );
    }

    println!("✅ Concurrent k-mer processing validated");
}

#[test]
fn test_adaptive_k_selection_correctness() {
    println!("📏 Testing adaptive k-mer selection correctness...");

    use meta_forge::assembly::adaptive_k::{AdaptiveKSelector, AdaptiveKConfig};

    // Test with different read characteristics
    let test_cases = vec![
        ("short_reads", create_short_reads(100, 25)),
        ("long_reads", create_long_reads(100, 100)),
        ("variable_reads", create_variable_length_reads(100)),
    ];

    for (test_name, reads) in test_cases {
        println!("  Testing {}", test_name);

        let config = AdaptiveKConfig {
            min_k: 15,
            max_k: 31,
            sample_size: 50,
            memory_budget_mb: 1024,
        };

        let selector = AdaptiveKSelector::new(config);
        let k = selector.select_optimal_k(&reads).unwrap();

        // K should be within bounds
        assert!(k >= 15 && k <= 31, "K should be within configured bounds: {}", k);

        // K should be odd (conventional for assembly)
        assert!(k % 2 == 1, "K should be odd: {}", k);

        // K should be reasonable for the read lengths
        let avg_read_len: f64 = reads.iter().map(|r| r.corrected.len()).sum::<usize>() as f64 / reads.len() as f64;
        assert!(
            k as f64 <= avg_read_len * 0.8,
            "K should not be too large relative to read length: k={}, avg_len={:.1}",
            k, avg_read_len
        );

        println!("    Selected k={} for avg read length {:.1}", k, avg_read_len);
    }

    println!("✅ Adaptive k-mer selection correctness validated");
}

// Helper functions

fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
    CorrectedRead {
        id,
        original: sequence.to_string(),
        corrected: sequence.to_string(),
        corrections: Vec::new(),
        quality_scores: vec![30; sequence.len()],
        correction_metadata: CorrectionMetadata {
            algorithm: "test".to_string(),
            confidence_threshold: 0.95,
            context_window: 5,
            correction_time_ms: 0,
        },
    }
}

fn generate_large_test_dataset(num_reads: usize) -> Vec<CorrectedRead> {
    let base_sequence = "ATCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGT";
    let mut reads = Vec::new();

    for i in 0..num_reads {
        let start = i % (base_sequence.len() - 30);
        let sequence = &base_sequence[start..start + 30];
        reads.push(create_test_read(i, sequence));
    }

    reads
}

fn create_deterministic_test_reads() -> Vec<CorrectedRead> {
    // Create reads with known, deterministic overlap pattern
    vec![
        create_test_read(0, "ATCGATCGATCGATCGATCGATCG"),
        create_test_read(1, "TCGATCGATCGATCGATCGATCGA"),
        create_test_read(2, "CGATCGATCGATCGATCGATCGAT"),
        create_test_read(3, "GATCGATCGATCGATCGATCGATC"),
        create_test_read(4, "ATCGATCGATCGATCGATCGATCG"), // Duplicate for coverage
    ]
}

fn create_short_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    let bases = ['A', 'T', 'C', 'G'];
    (0..count)
        .map(|i| {
            let sequence: String = (0..length)
                .map(|j| bases[(i + j) % 4])
                .collect();
            create_test_read(i, &sequence)
        })
        .collect()
}

fn create_long_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    let pattern = "ATCGATCGATCG";
    (0..count)
        .map(|i| {
            let sequence = pattern.repeat(length / pattern.len() + 1);
            let sequence = &sequence[..length];
            create_test_read(i, sequence)
        })
        .collect()
}

fn create_variable_length_reads(count: usize) -> Vec<CorrectedRead> {
    let base_pattern = "ATCGATCGATCGAAGTCGATCGATCGAAGT";
    (0..count)
        .map(|i| {
            let length = 20 + (i % 50); // Lengths from 20 to 69
            let sequence = &base_pattern[..length.min(base_pattern.len())];
            create_test_read(i, sequence)
        })
        .collect()
}

/// Integration test for complete assembly pipeline correctness
#[test]
fn test_complete_assembly_pipeline_correctness() {
    println!("🏁 Testing complete assembly pipeline correctness...");

    // Test different configurations
    let configs = vec![
        ("low_memory", LaptopConfig::low_memory()),
        ("medium_memory", LaptopConfig::medium_memory()),
        ("high_memory", LaptopConfig::high_memory()),
    ];

    for (config_name, config) in configs {
        println!("  Testing {} configuration...", config_name);

        let assembler = LaptopAssembler::new(config);
        let reads = generate_realistic_test_dataset();

        let start_time = std::time::Instant::now();
        let result = assembler.assemble(&reads);
        let elapsed = start_time.elapsed();

        match result {
            Ok(contigs) => {
                println!("    ✅ Assembly completed in {:.2}s with {} contigs",
                        elapsed.as_secs_f64(), contigs.len());

                // Validate assembly quality
                validate_assembly_quality(&contigs, &reads);
            }
            Err(e) => {
                // Allow timeout or memory errors for low-end configurations
                if config_name == "low_memory" && (e.to_string().contains("timeout") || e.to_string().contains("memory")) {
                    println!("    ⚠️ {} configuration hit resource limits (acceptable): {}", config_name, e);
                } else {
                    panic!("Unexpected assembly failure for {}: {}", config_name, e);
                }
            }
        }
    }

    println!("✅ Complete assembly pipeline correctness validated");
}

fn generate_realistic_test_dataset() -> Vec<CorrectedRead> {
    // Generate reads that simulate real metagenomic data
    let reference_segments = vec![
        "ATGCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGT",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG",
        "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG",
    ];

    let mut reads = Vec::new();
    let mut read_id = 0;

    for segment in &reference_segments {
        // Generate overlapping reads from each segment
        for i in 0..=(segment.len().saturating_sub(40)) {
            if i % 5 == 0 { // Every 5th position to create overlaps
                let read_seq = &segment[i..i + 40];
                reads.push(create_test_read(read_id, read_seq));
                read_id += 1;
            }
        }
    }

    reads
}

fn validate_assembly_quality(contigs: &[Contig], original_reads: &[CorrectedRead]) {
    // Basic quality checks
    assert!(!contigs.is_empty(), "Assembly should produce contigs");

    // Check total assembly length is reasonable
    let total_assembled: usize = contigs.iter().map(|c| c.length).sum();
    let total_read_length: usize = original_reads.iter().map(|r| r.corrected.len()).sum();

    assert!(
        total_assembled >= total_read_length / 10, // Should recover at least 10% of total read length
        "Assembly length too small: {} vs {} read length",
        total_assembled, total_read_length
    );

    // Check contigs are properly formed
    for contig in contigs {
        assert!(contig.length >= 15, "Contig too short: {}", contig.length);
        assert!(contig.coverage > 0.0, "Contig should have positive coverage");

        // Sequence should be valid DNA
        for c in contig.sequence.chars() {
            assert!(
                matches!(c, 'A' | 'T' | 'C' | 'G' | 'N'),
                "Invalid character in contig: {}",
                c
            );
        }
    }

    // Check for reasonable N50 (simplified)
    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
    lengths.sort_by(|a, b| b.cmp(a));

    if !lengths.is_empty() {
        let largest = lengths[0];
        assert!(largest >= 20, "Largest contig should be reasonable: {}", largest);
    }
}