//! Assembly Benchmark Validation Suite
//!
//! Validates that assembly optimizations maintain quality while improving performance.

use ahash::AHashMap;
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};
use meta_forge::utils::benchmark_validator::{
    BenchmarkConfig, BenchmarkValidator, PerformanceRequirements, QualityThresholds,
    ValidationOptions,
};
use std::collections::HashMap;

/// Create test dataset for benchmarking
fn create_benchmark_dataset(name: &str, size: usize) -> Vec<CorrectedRead> {
    (0..size)
        .map(|i| {
            // Generate synthetic metagenomic reads
            let seq = generate_synthetic_sequence(100 + (i % 50));
            CorrectedRead {
                id: i,
                original: seq.clone(),
                corrected: seq.clone(),
                corrections: Vec::new(),
                quality_scores: vec![35; seq.len()],
                correction_metadata: CorrectionMetadata {
                    algorithm: "benchmark".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: Vec::new(),
            }
        })
        .collect()
}

/// Generate synthetic metagenomic sequence
fn generate_synthetic_sequence(length: usize) -> String {
    const BASES: &[u8] = b"ACGT";
    (0..length)
        .map(|_| {
            let idx = fastrand::usize(..4);
            BASES[idx] as char
        })
        .collect()
}

#[test]
fn test_small_dataset_validation() {
    let mut validator = BenchmarkValidator::new();

    // Configure benchmark for small dataset
    validator.add_benchmark(BenchmarkConfig {
        dataset_name: "small_synthetic".to_string(),
        quality_thresholds: QualityThresholds {
            min_n50: 500,
            min_total_length: 2000,
            max_contigs: 100,
            min_avg_coverage: 2.0,
            max_gap_percentage: 10.0,
        },
        performance_requirements: PerformanceRequirements {
            max_execution_time_secs: 30,
            max_memory_usage_mb: 512,
            max_cpu_utilization: 95.0,
            min_speedup_factor: 1.0,
        },
        validation_options: ValidationOptions {
            iterations: 3,
            validate_against_reference: false,
            detailed_quality_analysis: true,
            compare_with_baseline: false,
        },
    });

    // Create test dataset
    let reads = create_benchmark_dataset("small_synthetic", 100);

    // Run benchmark
    let result = validator.run_benchmark("small_synthetic", &reads);
    assert!(result.is_ok(), "Benchmark should complete successfully");

    let benchmark_result = result.unwrap();
    println!("Benchmark Status: {:?}", benchmark_result.pass_status);
    println!(
        "Quality Passed: {}",
        benchmark_result.quality_validation.overall_quality_passed
    );
    println!(
        "Performance Passed: {}",
        benchmark_result
            .performance_validation
            .overall_performance_passed
    );
}

#[test]
fn test_medium_dataset_validation() {
    let mut validator = BenchmarkValidator::new();

    validator.add_benchmark(BenchmarkConfig {
        dataset_name: "medium_synthetic".to_string(),
        quality_thresholds: QualityThresholds {
            min_n50: 1000,
            min_total_length: 10000,
            max_contigs: 200,
            min_avg_coverage: 3.0,
            max_gap_percentage: 5.0,
        },
        performance_requirements: PerformanceRequirements {
            max_execution_time_secs: 120,
            max_memory_usage_mb: 1024,
            max_cpu_utilization: 90.0,
            min_speedup_factor: 1.0,
        },
        validation_options: ValidationOptions {
            iterations: 2,
            validate_against_reference: false,
            detailed_quality_analysis: true,
            compare_with_baseline: false,
        },
    });

    let reads = create_benchmark_dataset("medium_synthetic", 500);
    let result = validator.run_benchmark("medium_synthetic", &reads);
    assert!(result.is_ok());
}

#[test]
fn test_benchmark_suite() {
    let validator = BenchmarkValidator::default();

    let mut datasets = HashMap::new();
    datasets.insert(
        "small_test".to_string(),
        create_benchmark_dataset("small_test", 50),
    );

    let results = validator.run_all_benchmarks(&datasets);
    assert!(results.is_ok(), "Benchmark suite should complete");

    let benchmark_results = results.unwrap();
    assert!(
        !benchmark_results.is_empty(),
        "Should have benchmark results"
    );

    for (name, result) in benchmark_results.iter() {
        println!("\nBenchmark: {}", name);
        println!("  Status: {:?}", result.pass_status);
        println!(
            "  Speedup: {:.2}x",
            result.optimization_impact.speedup_factor
        );
        println!(
            "  Memory: {:.1}% change",
            result.optimization_impact.memory_reduction_percent
        );
    }
}
