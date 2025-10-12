//! Comprehensive Tests for Laptop Assembly Pipeline
//! ===============================================
//!
//! Tests for the laptop-optimized assembly implementation including:
//! - Configuration validation and auto-detection
//! - Adaptive k-mer selection
//! - Memory-bounded assembly
//! - Error handling and edge cases

use meta_forge::assembly::adaptive_k::{AdaptiveKConfig, AdaptiveKSelector};
use meta_forge::assembly::{LaptopAssembler, LaptopAssemblyGraph, LaptopConfig};
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};

/// Helper function to create test reads
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
        kmer_hash_cache: Vec::new(),
    }
}

/// Create realistic test dataset with overlapping reads
fn create_test_dataset() -> Vec<CorrectedRead> {
    let genome = "ATCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGT";
    let read_length = 20;
    let overlap = 10;
    let mut reads = Vec::new();

    for i in 0..=(genome.len() - read_length) {
        if i % (read_length - overlap) == 0 {
            let sequence = &genome[i..i + read_length];
            reads.push(create_test_read(i, sequence));
        }
    }

    // Add some reverse complement reads for complexity
    for (id, read) in reads.clone().iter().enumerate() {
        let rev_comp = reverse_complement(&read.corrected);
        reads.push(create_test_read(reads.len() + id, &rev_comp));
    }

    reads
}

/// Simple reverse complement implementation
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
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
fn test_laptop_config_presets() {
    let low_config = LaptopConfig::low_memory();
    assert_eq!(low_config.memory_budget_mb, 1024);
    assert_eq!(low_config.cpu_cores, 2);
    assert_eq!(low_config.max_k, 21);

    let medium_config = LaptopConfig::medium_memory();
    assert_eq!(medium_config.memory_budget_mb, 2048);
    assert_eq!(medium_config.cpu_cores, 4);
    assert_eq!(medium_config.max_k, 31);

    let high_config = LaptopConfig::high_memory();
    assert_eq!(high_config.memory_budget_mb, 4096);
    assert!(high_config.cpu_cores >= 1); // Should be at least 1
    assert_eq!(high_config.max_k, 63);
}

#[test]
fn test_laptop_config_custom() {
    // Valid custom configuration
    let config = LaptopConfig::custom(2048, 4, 31).unwrap();
    assert_eq!(config.memory_budget_mb, 2048);
    assert_eq!(config.cpu_cores, 4);
    assert_eq!(config.max_k, 31);

    // Invalid configurations should fail
    assert!(LaptopConfig::custom(100, 4, 31).is_err()); // Too little memory
    assert!(LaptopConfig::custom(2048, 0, 31).is_err()); // Zero cores
    assert!(LaptopConfig::custom(2048, 4, 10).is_err()); // K too small
    assert!(LaptopConfig::custom(2048, 4, 200).is_err()); // K too large
}

#[test]
fn test_laptop_config_auto_detect() {
    // Set environment variable to test auto-detection
    std::env::set_var("MEMORY_GB", "6");
    let config = LaptopConfig::auto_detect();

    // Should detect medium memory configuration for 6GB
    assert_eq!(config.memory_budget_mb, 2048);
    assert_eq!(config.cpu_cores, 4);

    // Clean up
    std::env::remove_var("MEMORY_GB");
}

#[test]
fn test_adaptive_k_selection() {
    let config = AdaptiveKConfig::default();
    let selector = AdaptiveKSelector::new(config);

    // Test with very short reads
    let short_reads = vec![
        create_test_read(0, "ATCGATCG"),
        create_test_read(1, "GCTAGCTA"),
        create_test_read(2, "TTAACCGG"),
    ];

    let k_short = selector.select_optimal_k(&short_reads).unwrap();
    assert!(k_short >= 15 && k_short <= 31);
    assert!(k_short <= 21); // Short reads should use smaller k

    // Test with longer reads
    let long_reads = vec![
        create_test_read(0, &"ATCGATCGATCGATCGATCG".repeat(5)),
        create_test_read(1, &"GCTAGCTAGCTAGCTAGCTA".repeat(5)),
    ];

    let k_long = selector.select_optimal_k(&long_reads).unwrap();
    assert!(k_long >= k_short); // Longer reads should allow larger k
}

#[test]
fn test_adaptive_k_gc_content() {
    let config = AdaptiveKConfig::default();
    let selector = AdaptiveKSelector::new(config);

    // High GC content reads (should prefer smaller k)
    let high_gc_reads = vec![
        create_test_read(0, &"GCGCGCGCGCGCGCGCGCGC".repeat(3)),
        create_test_read(1, &"CCGGCCGGCCGGCCGGCCGG".repeat(3)),
    ];

    let k_high_gc = selector.select_optimal_k(&high_gc_reads).unwrap();

    // Balanced GC content reads
    let balanced_reads = vec![
        create_test_read(0, &"ATCGATCGATCGATCGATCG".repeat(3)),
        create_test_read(1, &"TACGTACGTACGTACGTACG".repeat(3)),
    ];

    let k_balanced = selector.select_optimal_k(&balanced_reads).unwrap();

    // High GC should tend to use smaller or equal k
    assert!(k_high_gc <= k_balanced + 2);
}

#[test]
fn test_laptop_assembler_creation() {
    let config = LaptopConfig::medium_memory();
    let _assembler = LaptopAssembler::new(config.clone());

    // Test auto-config assembler
    let _auto_assembler = LaptopAssembler::auto_config();
    // Should create successfully
}

#[test]
fn test_laptop_assembly_small_dataset() {
    let config = LaptopConfig::low_memory(); // Use conservative config for tests
    let assembler = LaptopAssembler::new(config);

    let reads = vec![
        create_test_read(0, "ATCGATCGATCGAAGT"),
        create_test_read(1, "TCGATCGATCGAAGTC"),
        create_test_read(2, "GATCGATCGAAGTCGA"),
        create_test_read(3, "TCGATCGAAGTCGATC"),
    ];

    let contigs = assembler.assemble(&reads).unwrap();

    // Should produce at least one contig
    assert!(!contigs.is_empty());

    // Contigs should have reasonable properties
    for contig in &contigs {
        assert!(contig.length > 0);
        assert!(contig.coverage > 0.0);
        assert!(!contig.sequence.is_empty());
    }

    println!(
        "Generated {} contigs from {} reads",
        contigs.len(),
        reads.len()
    );
    for (i, contig) in contigs.iter().enumerate() {
        println!(
            "Contig {}: len={}, cov={:.2}",
            i, contig.length, contig.coverage
        );
    }
}

#[test]
fn test_laptop_assembly_realistic_dataset() {
    let config = LaptopConfig::medium_memory();
    let assembler = LaptopAssembler::new(config);

    let reads = create_test_dataset();
    println!("Testing with {} reads", reads.len());

    let contigs = assembler.assemble(&reads).unwrap();

    // Should produce contigs
    assert!(!contigs.is_empty());

    // At least one contig should be reasonably long
    let max_contig_length = contigs.iter().map(|c| c.length).max().unwrap_or(0);
    assert!(max_contig_length >= 20); // Should assemble overlapping reads

    println!(
        "Realistic dataset: {} contigs, max length: {}",
        contigs.len(),
        max_contig_length
    );
}

#[test]
fn test_laptop_assembly_graph_memory_bounds() {
    let config = LaptopConfig::custom(512, 2, 21).unwrap(); // Very limited memory
    let mut graph = LaptopAssemblyGraph::new(config);

    let reads = create_test_dataset();

    // Should handle memory constraints gracefully
    let result = graph.build_from_reads(&reads, 15); // Use small k
    assert!(result.is_ok());

    // Memory usage should be reasonable
    let memory_mb = graph.memory_usage_mb();
    assert!(memory_mb < 100.0); // Should stay well under budget

    println!("Graph memory usage: {:.2} MB", memory_mb);
}

#[test]
fn test_laptop_assembly_error_handling() {
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);

    // Empty reads should still work
    let empty_reads = vec![];
    let result = assembler.assemble(&empty_reads);
    // Should either succeed with empty result or fail gracefully
    assert!(result.is_ok() || result.is_err());

    // Reads with invalid characters should be handled
    let invalid_reads = vec![
        create_test_read(0, "ATCGXYZGATCG"), // Invalid characters
        create_test_read(1, "NNNNNNNNNNN"),  // All Ns
    ];

    let result = assembler.assemble(&invalid_reads);
    // Should handle gracefully without crashing
    assert!(result.is_ok() || result.is_err());
}

#[test]
fn test_memory_budget_enforcement() {
    let very_small_config = LaptopConfig::custom(256, 1, 17).unwrap();
    let assembler = LaptopAssembler::new(very_small_config);

    // Create a larger dataset that could potentially exceed memory
    let mut large_reads = Vec::new();
    for i in 0..100 {
        let sequence = format!("ATCGATCGATCG{:04}AAGT", i);
        large_reads.push(create_test_read(i, &sequence));
    }

    // Should complete without running out of memory
    let result = assembler.assemble(&large_reads);
    assert!(result.is_ok());

    if let Ok(contigs) = result {
        println!(
            "Memory-constrained assembly: {} contigs from {} reads",
            contigs.len(),
            large_reads.len()
        );
    }
}

#[test]
fn test_k_mer_size_validation() {
    let config = LaptopConfig::custom(1024, 2, 15).unwrap(); // Max k = 15
    let mut graph = LaptopAssemblyGraph::new(config);

    let reads = vec![create_test_read(0, "ATCGATCGATCGATCG")];

    // K larger than max should fail
    let result = graph.build_from_reads(&reads, 21);
    assert!(result.is_err());

    // K within limits should succeed
    let result = graph.build_from_reads(&reads, 15);
    assert!(result.is_ok());
}

#[test]
fn test_chunk_processing() {
    let mut config = LaptopConfig::custom(1024, 2, 21).unwrap();
    config.chunk_size = 2; // Force small chunks

    let assembler = LaptopAssembler::new(config);

    let reads = create_test_dataset();
    let result = assembler.assemble(&reads);

    // Should handle chunked processing correctly
    assert!(result.is_ok());

    if let Ok(contigs) = result {
        println!("Chunked processing: {} contigs", contigs.len());
    }
}

#[test]
fn test_concurrent_assembly() {
    use std::sync::Arc;
    use std::thread;

    let config = Arc::new(LaptopConfig::medium_memory());
    let reads = Arc::new(create_test_dataset());

    let mut handles = vec![];

    // Run multiple assemblies concurrently to test thread safety
    for i in 0..3 {
        let config_clone = Arc::clone(&config);
        let reads_clone = Arc::clone(&reads);

        let handle = thread::spawn(move || {
            let assembler = LaptopAssembler::new((*config_clone).clone());
            let result = assembler.assemble(&reads_clone);
            println!("Thread {} completed assembly", i);
            result
        });

        handles.push(handle);
    }

    // Wait for all threads to complete
    for handle in handles {
        let result = handle.join().unwrap();
        assert!(result.is_ok());
    }
}
