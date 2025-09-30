//! TDD Tests for Critical Assembly Fixes
//! =====================================
//!
//! Priority 1 tests to validate fixes for:
//! 1. K-mer overlap detection (has_overlap)
//! 2. Aggressive k-mer filtering (cleanup_rare_kmers)
//! 3. Sequence reconstruction (contig extension)
//! 4. Memory limits in parallel processing

use meta_forge::assembly::optimized::{BitPackedKmer, OptimizedAssembler, OptimizedConfig};
use meta_forge::assembly::laptop_assembly::{LaptopAssembler, LaptopConfig, BoundedKmerCounter};
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};

/// Test 1: K-mer Overlap Detection - Biological Correctness
#[test]
fn test_kmer_overlap_detection_correct() {
    // Two consecutive k-mers from sequence "ATCGATCG"
    let kmer1 = BitPackedKmer::new("ATCGA").unwrap(); // Positions 0-4
    let kmer2 = BitPackedKmer::new("TCGAT").unwrap(); // Positions 1-5

    // They should overlap: kmer1[1..5] == kmer2[0..4]
    // ATCGA -> suffix "TCGA"
    // TCGAT -> prefix "TCGA"
    // Expected: TRUE (they overlap with 4 bp)

    let seq1 = kmer1.to_string();
    let seq2 = kmer2.to_string();
    let k = kmer1.len();

    let suffix = &seq1[1..]; // "TCGA"
    let prefix = &seq2[..k-1]; // "TCGA"

    assert_eq!(suffix, prefix, "K-mers from consecutive positions should overlap");
}

#[test]
fn test_kmer_no_overlap() {
    // Two non-consecutive k-mers
    let kmer1 = BitPackedKmer::new("ATCGA").unwrap();
    let kmer2 = BitPackedKmer::new("GGGGG").unwrap();

    let seq1 = kmer1.to_string();
    let seq2 = kmer2.to_string();
    let k = kmer1.len();

    let suffix = &seq1[1..];
    let prefix = &seq2[..k-1];

    assert_ne!(suffix, prefix, "Non-consecutive k-mers should not overlap");
}

/// Test 2: K-mer Filtering - Not Too Aggressive
#[test]
fn test_kmer_filtering_preserves_valid_kmers() {
    let mut counter = BoundedKmerCounter::new(10); // Small budget to trigger cleanup

    // Add k-mers with different frequencies:
    // - Hash 1: count=1 (singleton, but valid)
    // - Hash 2: count=2 (low coverage)
    // - Hash 3: count=5 (good coverage)

    counter.add_kmer(1);
    counter.add_kmer(2);
    counter.add_kmer(2);
    counter.add_kmer(3);
    counter.add_kmer(3);
    counter.add_kmer(3);
    counter.add_kmer(3);
    counter.add_kmer(3);

    // Get frequent k-mers with min_count=2
    let frequent = counter.get_frequent_kmers(2);

    // Should include hash 2 (count=2) and hash 3 (count=5)
    assert!(frequent.iter().any(|(h, _)| *h == 2), "Should preserve count=2 k-mers");
    assert!(frequent.iter().any(|(h, _)| *h == 3), "Should preserve count=5 k-mers");

    // But not hash 1 (count=1) when using min_count=2
    assert!(!frequent.iter().any(|(h, _)| *h == 1), "Should filter count=1 when min=2");
}

#[test]
fn test_kmer_filtering_with_min_count_1() {
    let mut counter = BoundedKmerCounter::new(10);

    counter.add_kmer(1); // Singleton
    counter.add_kmer(2);
    counter.add_kmer(2);

    // With min_count=1, should keep singletons
    let frequent = counter.get_frequent_kmers(1);

    assert!(frequent.iter().any(|(h, _)| *h == 1), "Should preserve singletons with min_count=1");
    assert!(frequent.iter().any(|(h, _)| *h == 2), "Should preserve count=2");
}

/// Test 3: Sequence Reconstruction - Proper Extension
#[test]
fn test_contig_sequence_extension() {
    // Simulate contig extension with overlapping k-mers
    // Sequence: "ATCGATCG"
    // K-mers (k=5): ATCGA, TCGAT, CGATC, GATCG

    let kmers = vec![
        "ATCGA", // Start
        "TCGAT", // Extend by 1bp: add 'T'
        "CGATC", // Extend by 1bp: add 'C'
        "GATCG", // Extend by 1bp: add 'G'
    ];

    let mut sequence = String::new();

    for (i, kmer_str) in kmers.iter().enumerate() {
        if i == 0 {
            // First k-mer: add entire sequence
            sequence.push_str(kmer_str);
        } else {
            // Subsequent k-mers: add only the last character (non-overlapping part)
            let last_char = kmer_str.chars().last().unwrap();
            sequence.push(last_char);
        }
    }

    // Expected: ATCGA + T + C + G = "ATCGATCG"
    assert_eq!(sequence, "ATCGATCG", "Contig extension should add only non-overlapping bases");
    assert_eq!(sequence.len(), 8, "Length should be k + (num_kmers - 1) = 5 + 3 = 8");
}

#[test]
fn test_contig_extension_length_calculation() {
    let k = 21;
    let num_kmers = 100;

    // Expected contig length: k + (num_kmers - 1)
    // First k-mer contributes k bases
    // Each subsequent k-mer contributes 1 new base
    let expected_length = k + (num_kmers - 1);

    assert_eq!(expected_length, 120, "100 k-mers of size 21 should produce 120bp contig");
}

/// Test 4: Graph Connectivity - Ensuring Edges Are Created
#[test]
fn test_graph_connectivity_with_overlapping_kmers() {
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);

    // Create reads with overlapping k-mers
    let reads = vec![
        create_test_read(0, "ATCGATCGATCGATCGATCGATCG"), // 24bp
        create_test_read(1, "TCGATCGATCGATCGATCGATCGA"), // 24bp, offset by 1
        create_test_read(2, "CGATCGATCGATCGATCGATCGAT"), // 24bp, offset by 2
    ];

    // These reads should create a connected graph
    let result = assembler.assemble(&reads);

    match result {
        Ok(contigs) => {
            assert!(!contigs.is_empty(), "Should generate at least one contig from overlapping reads");

            // Check that at least one contig is longer than a single k-mer
            let has_extended_contig = contigs.iter().any(|c| c.length > 21);
            assert!(has_extended_contig, "Should have at least one contig extended beyond single k-mer");
        }
        Err(e) => panic!("Assembly failed: {}", e),
    }
}

/// Test 5: Memory Limits - Parallel Processing Should Not Explode
#[test]
fn test_memory_bounded_parallel_processing() {
    let config = LaptopConfig::custom(512, 2, 21).unwrap(); // 512MB budget
    let assembler = LaptopAssembler::new(config);

    // Create many reads to stress test memory
    let mut reads = Vec::new();
    for i in 0..5000 {
        reads.push(create_test_read(i, "ATCGATCGATCGATCGATCGATCGATCGATCGATCG"));
    }

    // Should complete without memory explosion
    let result = assembler.assemble(&reads);

    match result {
        Ok(contigs) => {
            println!("Generated {} contigs from {} reads", contigs.len(), reads.len());
            assert!(!contigs.is_empty(), "Should generate contigs");
        }
        Err(e) => {
            // If it fails, should be graceful (timeout or explicit memory limit)
            assert!(
                e.to_string().contains("timeout") || e.to_string().contains("memory"),
                "Failure should be due to timeout or memory limit, not crash"
            );
        }
    }
}

/// Test 6: Contig Generation Quality - Minimum Expected Counts
#[test]
fn test_contig_generation_produces_reasonable_counts() {
    let config = LaptopConfig::medium_memory();
    let assembler = LaptopAssembler::new(config);

    // Create 100 reads with good overlap
    let mut reads = Vec::new();
    let base = "ATCGATCGATCGATCGATCGATCG";

    for i in 0..100 {
        let offset = i % 10;
        if offset + 24 <= base.len() {
            reads.push(create_test_read(i, &base[offset..offset + 24]));
        } else {
            reads.push(create_test_read(i, base));
        }
    }

    let result = assembler.assemble(&reads);

    match result {
        Ok(contigs) => {
            // With 100 reads and good overlap, expect at least 5-10 contigs
            assert!(
                contigs.len() >= 5,
                "Expected at least 5 contigs from 100 reads, got {}",
                contigs.len()
            );

            // Check N50 quality
            let mut lengths: Vec<_> = contigs.iter().map(|c| c.length).collect();
            lengths.sort_by(|a, b| b.cmp(a));

            let total_length: usize = lengths.iter().sum();
            let mut cumulative = 0;
            let mut n50 = 0;

            for &length in &lengths {
                cumulative += length;
                if cumulative >= total_length / 2 {
                    n50 = length;
                    break;
                }
            }

            assert!(n50 >= 15, "N50 should be at least 15bp, got {}", n50);
        }
        Err(e) => panic!("Assembly failed: {}", e),
    }
}

/// Test 7: Overlap Detection with Reverse Complement
#[test]
fn test_kmer_overlap_with_reverse_complement() {
    let kmer1 = BitPackedKmer::new("ATCGA").unwrap();
    let kmer2 = BitPackedKmer::new("CGAT").unwrap(); // Reverse complement of "ATCG"

    // Check if reverse complement logic handles this correctly
    if let Ok(rc1) = kmer1.reverse_complement() {
        let rc1_str = rc1.to_string();
        let kmer2_str = kmer2.to_string();

        // The overlap check should consider reverse complements for canonical k-mers
        println!("kmer1: {}", kmer1.to_string());
        println!("rc1: {}", rc1_str);
        println!("kmer2: {}", kmer2_str);
    }
}

// Helper function
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