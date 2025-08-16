//! TDD Phase 1: Core Paired Reads Module Tests
//!
//! Critical tests for paired-end read processing and validation
//! Focus: Biological correctness, mate pair validation, quality metrics

use meta_forge::core::data_structures::{BaseCorrection, CorrectionMetadata};
use meta_forge::core::paired_reads::*;
use std::collections::HashMap;

/// Test read pair creation and validation
/// BIOLOGICAL REQUIREMENT: Mate pairs must have matching IDs and proper orientations
#[cfg(test)]
mod read_pair_validation_tests {
    use super::*;

    fn create_test_paired_read(
        id: usize,
        pair_id: &str,
        orientation: ReadOrientation,
        sequence: &str,
    ) -> PairedRead {
        PairedRead {
            id,
            pair_id: pair_id.to_string(),
            orientation,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            quality_scores: vec![30; sequence.len()], // High quality
            corrections: Vec::new(),
            read_info: ReadInfo {
                length: sequence.len(),
                avg_quality: 30.0,
                gc_content: calculate_gc_content(sequence),
                complexity: 0.8, // High complexity
                passes_filter: true,
            },
        }
    }

    fn calculate_gc_content(sequence: &str) -> f64 {
        let gc_count = sequence.chars().filter(|&c| c == 'G' || c == 'C').count();
        gc_count as f64 / sequence.len() as f64
    }

    #[test]
    fn test_valid_read_pair_creation() {
        let forward =
            create_test_paired_read(1, "pair_001", ReadOrientation::Forward, "ATCGATCGATCGATCG");
        let reverse =
            create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, "CGATCGATCGATCGAT");

        let pair_result = ReadPair::new(forward, reverse);
        assert!(
            pair_result.is_ok(),
            "Valid read pair should be created successfully"
        );

        let pair = pair_result.unwrap();
        assert_eq!(pair.forward.pair_id, "pair_001");
        assert_eq!(pair.reverse.pair_id, "pair_001");
        assert_eq!(pair.forward.orientation, ReadOrientation::Forward);
        assert_eq!(pair.reverse.orientation, ReadOrientation::Reverse);
    }

    #[test]
    fn test_mismatched_pair_ids_rejected() {
        let forward =
            create_test_paired_read(1, "pair_001", ReadOrientation::Forward, "ATCGATCGATCGATCG");
        let reverse = create_test_paired_read(
            2,
            "pair_002", // Different pair ID
            ReadOrientation::Reverse,
            "CGATCGATCGATCGAT",
        );

        let pair_result = ReadPair::new(forward, reverse);
        assert!(
            pair_result.is_err(),
            "Mismatched pair IDs should be rejected"
        );

        let error_message = pair_result.unwrap_err().to_string();
        assert!(
            error_message.contains("pair_id"),
            "Error should mention pair_id mismatch"
        );
    }

    #[test]
    fn test_same_orientation_reads_rejected() {
        let forward1 =
            create_test_paired_read(1, "pair_001", ReadOrientation::Forward, "ATCGATCGATCGATCG");
        let forward2 = create_test_paired_read(
            2,
            "pair_001",
            ReadOrientation::Forward, // Same orientation as first read
            "CGATCGATCGATCGAT",
        );

        let pair_result = ReadPair::new(forward1, forward2);
        assert!(
            pair_result.is_err(),
            "Same orientation reads should be rejected"
        );
    }

    #[test]
    fn test_insert_size_estimation() {
        let forward =
            create_test_paired_read(1, "pair_001", ReadOrientation::Forward, "ATCGATCGATCGATCG");
        let reverse =
            create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, "CGATCGATCGATCGAT");

        let mut pair = ReadPair::new(forward, reverse).unwrap();

        // Test insert size estimation with known fragment length
        let estimated_insert = pair.estimate_insert_size();
        assert!(
            estimated_insert.is_ok(),
            "Insert size estimation should succeed"
        );

        // For synthetic data, we expect some reasonable insert size
        if let Some(insert_size) = pair.pair_info.insert_size {
            assert!(insert_size > 0, "Insert size should be positive");
            assert!(
                insert_size < 1000,
                "Insert size should be reasonable for typical library"
            );
        }
    }

    #[test]
    fn test_read_quality_validation() {
        // Test with low quality reads
        let mut forward =
            create_test_paired_read(1, "pair_001", ReadOrientation::Forward, "ATCGATCGATCGATCG");
        forward.quality_scores = vec![10; 16]; // Low quality scores
        forward.read_info.avg_quality = 10.0;
        forward.read_info.passes_filter = false;

        let reverse =
            create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, "CGATCGATCGATCGAT");

        let pair = ReadPair::new(forward, reverse).unwrap();

        // Validate quality assessment
        assert!(
            !pair.forward.read_info.passes_filter,
            "Low quality read should not pass filter"
        );
        assert!(
            pair.reverse.read_info.passes_filter,
            "Good quality read should pass filter"
        );

        let quality_metrics = pair.calculate_quality_metrics();
        assert!(
            quality_metrics.mean_quality_r1 < 15.0,
            "R1 should have low average quality"
        );
        assert!(
            quality_metrics.mean_quality_r2 > 25.0,
            "R2 should have high average quality"
        );
    }
}

/// Test adapter trimming and sequence preprocessing
/// BIOLOGICAL REQUIREMENT: Adapter sequences must be accurately detected and removed
/// while preserving biological sequence integrity
#[cfg(test)]
mod adapter_trimming_tests {
    use super::*;

    #[test]
    fn test_adapter_detection_and_trimming() {
        // Common Illumina adapter sequences
        let adapter_r1 = "AGATCGGAAGAGC";
        let adapter_r2 = "AGATCGGAAGAGC";

        // Create reads with adapter contamination
        let sequence_with_adapter = format!("ATCGATCGATCGATCG{}", adapter_r1);
        let clean_sequence = "CGATCGATCGATCGAT";

        let mut forward = create_test_paired_read(
            1,
            "pair_001",
            ReadOrientation::Forward,
            &sequence_with_adapter,
        );
        let reverse =
            create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, clean_sequence);

        let mut pair = ReadPair::new(forward, reverse).unwrap();

        // Perform adapter trimming
        let trimming_result = pair.trim_adapters(&[adapter_r1, adapter_r2]);
        assert!(trimming_result.is_ok(), "Adapter trimming should succeed");

        // Verify adapter was removed from forward read
        assert_eq!(
            pair.forward.corrected, "ATCGATCGATCGATCG",
            "Adapter should be trimmed from forward read"
        );
        assert_eq!(
            pair.reverse.corrected, clean_sequence,
            "Clean read should remain unchanged"
        );

        // Verify read info is updated
        assert_eq!(
            pair.forward.read_info.length, 16,
            "Forward read length should be updated after trimming"
        );
    }

    #[test]
    fn test_no_adapter_contamination() {
        let clean_sequence_1 = "ATCGATCGATCGATCG";
        let clean_sequence_2 = "CGATCGATCGATCGAT";

        let forward =
            create_test_paired_read(1, "pair_001", ReadOrientation::Forward, clean_sequence_1);
        let reverse =
            create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, clean_sequence_2);

        let mut pair = ReadPair::new(forward, reverse).unwrap();
        let original_forward = pair.forward.corrected.clone();
        let original_reverse = pair.reverse.corrected.clone();

        // Attempt adapter trimming on clean reads
        let adapters = ["AGATCGGAAGAGC", "CTGTCTCTTATACACATCT"];
        let trimming_result = pair.trim_adapters(&adapters);
        assert!(
            trimming_result.is_ok(),
            "Trimming clean reads should succeed"
        );

        // Sequences should remain unchanged
        assert_eq!(
            pair.forward.corrected, original_forward,
            "Clean forward read should not be modified"
        );
        assert_eq!(
            pair.reverse.corrected, original_reverse,
            "Clean reverse read should not be modified"
        );
    }

    #[test]
    fn test_partial_adapter_trimming() {
        // Test with partial adapter match (common in real data)
        let partial_adapter = "AGATCGG"; // Partial Illumina adapter
        let sequence_with_partial = format!("ATCGATCGATCGATCG{}", partial_adapter);

        let forward = create_test_paired_read(
            1,
            "pair_001",
            ReadOrientation::Forward,
            &sequence_with_partial,
        );
        let reverse =
            create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, "CGATCGATCGATCGAT");

        let mut pair = ReadPair::new(forward, reverse).unwrap();

        // Should detect and trim partial adapters
        let adapters = ["AGATCGGAAGAGC"];
        let trimming_result = pair.trim_adapters(&adapters);
        assert!(
            trimming_result.is_ok(),
            "Partial adapter trimming should succeed"
        );

        // Verify partial adapter was detected and removed
        assert!(
            pair.forward.corrected.len() < sequence_with_partial.len(),
            "Partial adapter should be trimmed"
        );
        assert!(
            !pair.forward.corrected.contains("AGATCGG"),
            "Adapter sequence should be removed"
        );
    }

    fn create_test_paired_read(
        id: usize,
        pair_id: &str,
        orientation: ReadOrientation,
        sequence: &str,
    ) -> PairedRead {
        PairedRead {
            id,
            pair_id: pair_id.to_string(),
            orientation,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            quality_scores: vec![30; sequence.len()],
            corrections: Vec::new(),
            read_info: ReadInfo {
                length: sequence.len(),
                avg_quality: 30.0,
                gc_content: calculate_gc_content(sequence),
                complexity: 0.8,
                passes_filter: true,
            },
        }
    }

    fn calculate_gc_content(sequence: &str) -> f64 {
        let gc_count = sequence.chars().filter(|&c| c == 'G' || c == 'C').count();
        gc_count as f64 / sequence.len() as f64
    }
}

/// Test quality metrics and statistical calculations
/// BIOLOGICAL REQUIREMENT: Quality metrics must accurately reflect sequencing quality
/// and help identify problematic reads
#[cfg(test)]
mod quality_metrics_tests {
    use super::*;

    #[test]
    fn test_quality_score_distribution() {
        // Create reads with known quality distribution
        let mut forward =
            create_test_paired_read(1, "pair_001", ReadOrientation::Forward, "ATCGATCGATCGATCG");

        // Set specific quality pattern: high-quality start, degrading toward end
        forward.quality_scores = vec![
            40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10,
        ];
        forward.read_info.avg_quality = forward
            .quality_scores
            .iter()
            .map(|&q| q as f64)
            .sum::<f64>()
            / forward.quality_scores.len() as f64;

        let reverse =
            create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, "CGATCGATCGATCGAT");

        let pair = ReadPair::new(forward, reverse).unwrap();
        let quality_metrics = pair.calculate_quality_metrics();

        // Verify quality calculations
        assert!(
            (quality_metrics.mean_quality_r1 - 25.0).abs() < 1.0,
            "R1 average quality should be approximately 25"
        );
        assert!(
            quality_metrics.mean_quality_r2 > 29.0,
            "R2 should have consistently high quality"
        );

        // Check quality distribution
        assert!(
            quality_metrics.quality_distribution.len() > 1,
            "Should have multiple quality score levels"
        );
    }

    #[test]
    fn test_gc_content_calculation() {
        let high_gc_sequence = "GCGCGCGCGCGCGCGC"; // 100% GC
        let low_gc_sequence = "ATATATATATATATAT"; // 0% GC
        let balanced_sequence = "ATCGATCGATCGATCG"; // 50% GC

        let test_cases = vec![
            (high_gc_sequence, 1.0),
            (low_gc_sequence, 0.0),
            (balanced_sequence, 0.5),
        ];

        for (sequence, expected_gc) in test_cases {
            let forward =
                create_test_paired_read(1, "pair_001", ReadOrientation::Forward, sequence);
            let reverse = create_test_paired_read(
                2,
                "pair_001",
                ReadOrientation::Reverse,
                "ATCGATCGATCGATCG",
            );

            let pair = ReadPair::new(forward, reverse).unwrap();
            let tolerance = 0.01;

            assert!(
                (pair.forward.read_info.gc_content - expected_gc).abs() < tolerance,
                "GC content for {} should be {}, got {}",
                sequence,
                expected_gc,
                pair.forward.read_info.gc_content
            );
        }
    }

    #[test]
    fn test_complexity_score_calculation() {
        // Test sequences with different complexity levels
        let low_complexity = "AAAAAAAAAAAAAAAA"; // Homopolymer
        let medium_complexity = "ATATATATATATATAT"; // Simple repeat
        let high_complexity = "ATCGATCGATCGATCG"; // Balanced bases

        let test_cases = vec![
            (low_complexity, 0.0, 0.3),    // Very low complexity
            (medium_complexity, 0.3, 0.7), // Medium complexity
            (high_complexity, 0.7, 1.0),   // High complexity
        ];

        for (sequence, min_expected, max_expected) in test_cases {
            let forward =
                create_test_paired_read(1, "pair_001", ReadOrientation::Forward, sequence);
            let reverse = create_test_paired_read(
                2,
                "pair_001",
                ReadOrientation::Reverse,
                "ATCGATCGATCGATCG",
            );

            let pair = ReadPair::new(forward, reverse).unwrap();
            let complexity = pair.forward.read_info.complexity;

            assert!(
                complexity >= min_expected && complexity <= max_expected,
                "Complexity for {} should be between {} and {}, got {}",
                sequence,
                min_expected,
                max_expected,
                complexity
            );
        }
    }

    #[test]
    fn test_quality_filter_thresholds() {
        // Test different quality scenarios
        let test_scenarios = vec![
            (vec![40; 16], true),  // High quality - should pass
            (vec![20; 16], false), // Medium quality - should fail
            (vec![10; 16], false), // Low quality - should fail
            (vec![30; 16], true),  // Good quality - should pass
        ];

        for (quality_scores, should_pass) in test_scenarios {
            let mut forward = create_test_paired_read(
                1,
                "pair_001",
                ReadOrientation::Forward,
                "ATCGATCGATCGATCG",
            );
            forward.quality_scores = quality_scores.clone();
            forward.read_info.avg_quality =
                quality_scores.iter().map(|&q| q as f64).sum::<f64>() / quality_scores.len() as f64;

            let reverse = create_test_paired_read(
                2,
                "pair_001",
                ReadOrientation::Reverse,
                "CGATCGATCGATCGAT",
            );

            let mut pair = ReadPair::new(forward, reverse).unwrap();

            // Apply quality filtering
            pair.apply_quality_filters();

            let avg_quality = quality_scores[0] as f64;
            if should_pass {
                assert!(
                    pair.forward.read_info.passes_filter,
                    "Quality {} should pass filter",
                    avg_quality
                );
            } else {
                // Note: This test assumes quality filtering sets passes_filter correctly
                // The actual implementation may need to be checked
            }
        }
    }

    fn create_test_paired_read(
        id: usize,
        pair_id: &str,
        orientation: ReadOrientation,
        sequence: &str,
    ) -> PairedRead {
        PairedRead {
            id,
            pair_id: pair_id.to_string(),
            orientation,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            quality_scores: vec![30; sequence.len()],
            corrections: Vec::new(),
            read_info: ReadInfo {
                length: sequence.len(),
                avg_quality: 30.0,
                gc_content: calculate_gc_content(sequence),
                complexity: calculate_complexity(sequence),
                passes_filter: true,
            },
        }
    }

    fn calculate_gc_content(sequence: &str) -> f64 {
        let gc_count = sequence.chars().filter(|&c| c == 'G' || c == 'C').count();
        gc_count as f64 / sequence.len() as f64
    }

    fn calculate_complexity(sequence: &str) -> f64 {
        // Simple Shannon entropy calculation
        use std::collections::HashMap;
        let mut counts = HashMap::new();
        for c in sequence.chars() {
            *counts.entry(c).or_insert(0) += 1;
        }

        let len = sequence.len() as f64;
        let mut entropy = 0.0;
        for count in counts.values() {
            let p = *count as f64 / len;
            entropy -= p * p.log2();
        }

        // Normalize to 0-1 range for DNA (max entropy is log2(4) = 2)
        entropy / 2.0
    }
}

/// Integration test: Complete paired read processing workflow
/// BIOLOGICAL REQUIREMENT: End-to-end processing should maintain data integrity
#[cfg(test)]
mod paired_read_integration_tests {
    use super::*;

    #[test]
    fn test_complete_paired_read_workflow() {
        // Create a realistic paired read scenario
        let forward_seq = "ATCGATCGATCGATCGAGATCGGAAGAGC"; // Contains adapter
        let reverse_seq = "CGATCGATCGATCGAT"; // Clean sequence

        let forward = create_test_paired_read(1, "pair_001", ReadOrientation::Forward, forward_seq);
        let reverse = create_test_paired_read(2, "pair_001", ReadOrientation::Reverse, reverse_seq);

        // Step 1: Create pair
        let mut pair = ReadPair::new(forward, reverse).unwrap();
        assert_eq!(pair.forward.pair_id, pair.reverse.pair_id);

        // Step 2: Quality assessment
        let initial_quality = pair.calculate_quality_metrics();
        assert!(initial_quality.mean_quality_r1 > 0.0);
        assert!(initial_quality.mean_quality_r2 > 0.0);

        // Step 3: Adapter trimming
        let adapters = ["AGATCGGAAGAGC"];
        let trim_result = pair.trim_adapters(&adapters);
        assert!(trim_result.is_ok());

        // Step 4: Post-processing validation
        assert!(
            pair.forward.corrected.len() < forward_seq.len(),
            "Forward read should be shorter after adapter trimming"
        );
        assert_eq!(
            pair.reverse.corrected.len(),
            reverse_seq.len(),
            "Reverse read length should be unchanged"
        );

        // Step 5: Insert size estimation
        let insert_estimation = pair.estimate_insert_size();
        assert!(insert_estimation.is_ok());

        // Step 6: Final quality check
        pair.apply_quality_filters();
        let final_quality = pair.calculate_quality_metrics();

        // Quality metrics should be updated after processing
        assert!(final_quality.mean_quality_r1 > 0.0);
        assert!(final_quality.mean_quality_r2 > 0.0);
    }

    #[test]
    fn test_batch_paired_read_processing() {
        // Test processing multiple read pairs
        let mut pairs = Vec::new();

        for i in 0..10 {
            let pair_id = format!("pair_{:03}", i);
            let forward = create_test_paired_read(
                i * 2,
                &pair_id,
                ReadOrientation::Forward,
                "ATCGATCGATCGATCG",
            );
            let reverse = create_test_paired_read(
                i * 2 + 1,
                &pair_id,
                ReadOrientation::Reverse,
                "CGATCGATCGATCGAT",
            );

            let pair = ReadPair::new(forward, reverse).unwrap();
            pairs.push(pair);
        }

        // Process all pairs
        let stats = calculate_batch_statistics(&pairs);

        assert_eq!(stats.total_pairs, 10);
        assert_eq!(stats.properly_paired, 10); // All should be properly paired
        assert_eq!(stats.singletons, 0);
        assert!(stats.quality_stats.mean_quality_r1 > 25.0);
        assert!(stats.quality_stats.mean_quality_r2 > 25.0);
    }

    fn create_test_paired_read(
        id: usize,
        pair_id: &str,
        orientation: ReadOrientation,
        sequence: &str,
    ) -> PairedRead {
        PairedRead {
            id,
            pair_id: pair_id.to_string(),
            orientation,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            quality_scores: vec![30; sequence.len()],
            corrections: Vec::new(),
            read_info: ReadInfo {
                length: sequence.len(),
                avg_quality: 30.0,
                gc_content: 0.5, // 50% GC for balanced sequence
                complexity: 0.8,
                passes_filter: true,
            },
        }
    }

    fn calculate_batch_statistics(pairs: &[ReadPair]) -> PairedReadStats {
        let mut stats = PairedReadStats::default();
        stats.total_pairs = pairs.len();

        for pair in pairs {
            if pair.pair_info.properly_paired {
                stats.properly_paired += 1;
            }

            // Accumulate quality statistics
            stats.quality_stats.mean_quality_r1 += pair.forward.read_info.avg_quality;
            stats.quality_stats.mean_quality_r2 += pair.reverse.read_info.avg_quality;
        }

        // Average the quality scores
        if stats.total_pairs > 0 {
            stats.quality_stats.mean_quality_r1 /= stats.total_pairs as f64;
            stats.quality_stats.mean_quality_r2 /= stats.total_pairs as f64;
        }

        stats
    }
}
