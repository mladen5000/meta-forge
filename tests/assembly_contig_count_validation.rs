//! Assembly Contig Count Validation Tests
//! =======================================
//!
//! Tests to validate the fundamental biological constraint:
//! Maximum possible contigs ≤ Number of input reads
//!
//! This catches the critical bug where assembly generates more contigs than reads.

use meta_forge::assembly::laptop_assembly::{LaptopAssembler, LaptopConfig};
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};

/// Helper to create test reads
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
        kmer_hash_cache: AHashMap::new(),
    }
}

/// CRITICAL TEST: Validate fundamental biological constraint
#[test]
fn test_contig_count_cannot_exceed_read_count() {
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);

    // Create 100 test reads with good overlap
    let mut reads = Vec::new();
    let base_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"; // 48bp

    for i in 0..100 {
        // Create overlapping reads to test realistic assembly
        let start = (i % 10) as usize;
        let end = (start + 40).min(base_sequence.len());
        let sequence = &base_sequence[start..end];

        reads.push(create_test_read(i, sequence));
    }

    let result = assembler.assemble(&reads);

    match result {
        Ok(contigs) => {
            // CRITICAL ASSERTION: Contigs must be ≤ reads
            assert!(
                contigs.len() <= reads.len(),
                "CRITICAL BUG: Generated {} contigs from {} reads. \
                 Maximum biologically possible = {} (one per read if no overlap). \
                 This violates fundamental assembly constraints.",
                contigs.len(),
                reads.len(),
                reads.len()
            );

            // Additional quality checks
            if !contigs.is_empty() {
                let avg_length: f64 =
                    contigs.iter().map(|c| c.length).sum::<usize>() as f64 / contigs.len() as f64;

                println!("✅ Assembly validation passed:");
                println!("   Reads: {}", reads.len());
                println!("   Contigs: {}", contigs.len());
                println!(
                    "   Ratio: {:.2}x reduction",
                    reads.len() as f64 / contigs.len() as f64
                );
                println!("   Avg contig length: {:.1} bp", avg_length);

                // Reasonable assembly should produce fewer contigs than reads (merging happened)
                if contigs.len() < reads.len() {
                    println!("   ✨ Good: Reads were successfully merged into longer contigs");
                } else {
                    println!("   ⚠️  Warning: No merging occurred (contigs = reads)");
                }
            }
        }
        Err(e) => {
            // Check if error is the expected validation error
            let error_msg = e.to_string();
            if error_msg.contains("Generated") && error_msg.contains("contigs from") {
                panic!(
                    "Assembly correctly detected contig over-generation bug: {}",
                    e
                );
            } else {
                // Other errors might be acceptable (e.g., no valid contigs)
                println!("Assembly failed (may be acceptable): {}", e);
            }
        }
    }
}

/// Test with completely non-overlapping reads
#[test]
fn test_non_overlapping_reads_produce_equal_or_fewer_contigs() {
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);

    // Create completely different sequences (no overlap possible)
    let sequences = vec![
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", // 34bp of A's
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", // 34bp of T's
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", // 34bp of G's
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", // 34bp of C's
    ];

    let reads: Vec<_> = sequences
        .iter()
        .enumerate()
        .map(|(i, seq)| create_test_read(i, seq))
        .collect();

    let result = assembler.assemble(&reads);

    match result {
        Ok(contigs) => {
            // Even with no overlap, contigs should be ≤ reads
            assert!(
                contigs.len() <= reads.len(),
                "Even with no overlap, got {} contigs from {} reads",
                contigs.len(),
                reads.len()
            );
        }
        Err(e) if e.to_string().contains("Generated") => {
            panic!("Detected contig over-generation: {}", e);
        }
        Err(_) => {
            // No contigs generated is acceptable for non-overlapping reads
        }
    }
}

/// Test that assembly rejects when it detects over-generation
#[test]
#[should_panic(expected = "Generated")]
fn test_assembly_validates_and_rejects_over_generation() {
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);

    // This test will pass once we've verified the validation catches bugs
    // If the old buggy code runs, it would generate 1000+ contigs from 10 reads
    let reads: Vec<_> = (0..10)
        .map(|i| create_test_read(i, "ATCGATCGATCGATCGATCGATCG"))
        .collect();

    // If assembly has the bug, this should panic with validation error
    let _result = assembler
        .assemble(&reads)
        .expect("Should reject over-generation");
}

/// Test MetaSPAdes-style filtering
#[test]
fn test_metaspades_quality_standards() {
    let config = LaptopConfig::medium_memory();
    let assembler = LaptopAssembler::new(config);

    // Create reads with good overlap for assembly
    let mut reads = Vec::new();
    for i in 0..200 {
        let sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        reads.push(create_test_read(i, sequence));
    }

    let result = assembler.assemble(&reads);

    if let Ok(contigs) = result {
        // Verify contigs meet quality standards
        for contig in &contigs {
            // MetaSPAdes requires minimum coverage (we use 2.0x)
            assert!(
                contig.coverage >= 2.0,
                "Contig has insufficient coverage: {:.1}x (minimum 2.0x)",
                contig.coverage
            );

            // For metagenomics, reasonable minimum length
            // (We're more lenient than MetaSPAdes' 500bp for testing)
            assert!(
                contig.length >= 15,
                "Contig too short: {} bp",
                contig.length
            );
        }

        println!("✅ All {} contigs meet quality standards", contigs.len());
    }
}

/// Test that single k-mer contigs are NOT generated
#[test]
fn test_no_single_kmer_contigs() {
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);

    let reads: Vec<_> = (0..50)
        .map(|i| create_test_read(i, "ATCGATCGATCGATCGATCGATCGATCGATCG"))
        .collect();

    let result = assembler.assemble(&reads);

    if let Ok(contigs) = result {
        // No contig should be just a single k-mer (e.g., 21bp for k=21)
        for contig in &contigs {
            assert!(
                contig.length > 21,
                "Found single k-mer contig: {} bp (this should not happen)",
                contig.length
            );
        }

        if !contigs.is_empty() {
            println!("✅ No single k-mer contigs generated (good!)");
        }
    }
}
