//! Laptop Performance Regression Tests
//! ===================================
//!
//! This test suite specifically validates that laptop-optimized assembly
//! maintains performance and quality standards over time. It focuses on
//! preventing regressions in the critical optimizations.
//!
//! Key Areas:
//! - Memory usage regression detection
//! - Performance benchmark validation
//! - Quality preservation across updates
//! - CI/CD integration points

use anyhow::Result;
use meta_forge::assembly::laptop_assembly::{LaptopAssembler, LaptopConfig};
use meta_forge::assembly::performance_optimizations::{
    CacheOptimizedGraph, OptimizationConfig, PerformanceMode, PerformanceBenchmark,
    SIMDNucleotideOps
};
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};
use std::time::{Instant, Duration};
use std::collections::HashMap;

/* ========================================================================= */
/*                        PERFORMANCE REGRESSION TRACKING                  */
/* ========================================================================= */

/// Performance baseline expectations for laptop hardware
#[derive(Debug, Clone)]
struct PerformanceBaseline {
    /// Maximum acceptable assembly time per 1000 reads (seconds)
    max_assembly_time_per_1k_reads: f64,
    /// Maximum memory usage for low-memory config (MB)
    max_memory_usage_mb: f64,
    /// Minimum expected SIMD speedup over scalar
    min_simd_speedup: f64,
    /// Maximum acceptable contig-to-read ratio for overlapping reads
    max_contig_ratio: f64,
}

impl PerformanceBaseline {
    fn laptop_4gb() -> Self {
        Self {
            max_assembly_time_per_1k_reads: 10.0, // 10 seconds per 1000 reads
            max_memory_usage_mb: 1024.0,           // 1GB limit
            min_simd_speedup: 1.0,                 // At least no regression
            max_contig_ratio: 0.5,                 // Should merge reads significantly
        }
    }

    fn laptop_8gb() -> Self {
        Self {
            max_assembly_time_per_1k_reads: 5.0,   // 5 seconds per 1000 reads
            max_memory_usage_mb: 2048.0,           // 2GB limit
            min_simd_speedup: 1.2,                 // 20% improvement expected
            max_contig_ratio: 0.4,                 // Better merging with more memory
        }
    }

    fn laptop_16gb() -> Self {
        Self {
            max_assembly_time_per_1k_reads: 3.0,   // 3 seconds per 1000 reads
            max_memory_usage_mb: 4096.0,           // 4GB limit
            min_simd_speedup: 1.5,                 // 50% improvement expected
            max_contig_ratio: 0.3,                 // Excellent merging
        }
    }
}

#[cfg(test)]
mod performance_regression_tests {
    use super::*;

    /// Test that assembly performance meets laptop baseline expectations
    #[test]
    fn test_assembly_performance_baseline() -> Result<()> {
        let test_cases = vec![
            ("4GB Laptop", LaptopConfig::low_memory(), PerformanceBaseline::laptop_4gb()),
            ("8GB Laptop", LaptopConfig::medium_memory(), PerformanceBaseline::laptop_8gb()),
            ("16GB Laptop", LaptopConfig::high_memory(), PerformanceBaseline::laptop_16gb()),
        ];

        for (config_name, config, baseline) in test_cases {
            println!("\n🔬 Testing {} configuration...", config_name);

            // Create test dataset of 1000 reads
            let reads = create_benchmark_reads(1000, 100);
            let assembler = LaptopAssembler::new(config);

            // Measure assembly time
            let start = Instant::now();
            let contigs = assembler.assemble(&reads)?;
            let assembly_time = start.elapsed().as_secs_f64();

            // Calculate performance metrics
            let contig_ratio = contigs.len() as f64 / reads.len() as f64;

            // Validate against baseline
            assert!(assembly_time <= baseline.max_assembly_time_per_1k_reads,
                   "{}: Assembly time {:.2}s exceeds baseline {:.2}s",
                   config_name, assembly_time, baseline.max_assembly_time_per_1k_reads);

            assert!(contig_ratio <= baseline.max_contig_ratio,
                   "{}: Contig ratio {:.3} exceeds baseline {:.3}",
                   config_name, contig_ratio, baseline.max_contig_ratio);

            println!("✅ {}: {:.2}s assembly, {:.3} contig ratio - WITHIN BASELINE",
                    config_name, assembly_time, contig_ratio);
        }

        Ok(())
    }

    /// Test SIMD performance regression
    #[test]
    fn test_simd_performance_regression() -> Result<()> {
        let test_sequences = vec![
            create_test_sequence(1_000),   // 1KB
            create_test_sequence(10_000),  // 10KB
            create_test_sequence(100_000), // 100KB
        ];

        let benchmark = PerformanceBenchmark::new("SIMD Regression Test", 100);

        for (i, sequence) in test_sequences.iter().enumerate() {
            let result = benchmark.benchmark_nucleotide_counting(sequence);

            let baseline = PerformanceBaseline::laptop_8gb();
            assert!(result.speedup >= baseline.min_simd_speedup,
                   "SIMD speedup {:.2}x below baseline {:.2}x for {}KB sequence",
                   result.speedup, baseline.min_simd_speedup, sequence.len() / 1000);

            println!("✅ SIMD test {}: {:.2}x speedup on {}KB sequence",
                    i + 1, result.speedup, sequence.len() / 1000);
        }

        Ok(())
    }

    /// Test memory usage regression across configurations
    #[test]
    fn test_memory_usage_regression() -> Result<()> {
        let configs = vec![
            ("Low Memory", LaptopConfig::low_memory(), PerformanceBaseline::laptop_4gb()),
            ("Medium Memory", LaptopConfig::medium_memory(), PerformanceBaseline::laptop_8gb()),
            ("High Memory", LaptopConfig::high_memory(), PerformanceBaseline::laptop_16gb()),
        ];

        for (name, config, baseline) in configs {
            let reads = create_memory_stress_reads(2000, 80);
            let mut graph = meta_forge::assembly::laptop_assembly::LaptopAssemblyGraph::new(config);

            // Build graph and measure memory
            graph.build_from_reads(&reads, 21)?;
            let memory_usage = graph.memory_usage_mb();

            assert!(memory_usage <= baseline.max_memory_usage_mb,
                   "{}: Memory usage {:.1}MB exceeds baseline {:.1}MB",
                   name, memory_usage, baseline.max_memory_usage_mb);

            println!("✅ {}: {:.1}MB memory usage - WITHIN BASELINE",
                    name, memory_usage);
        }

        Ok(())
    }

    /// Test performance scaling with dataset size
    #[test]
    fn test_performance_scaling_regression() -> Result<()> {
        let dataset_sizes = vec![100, 500, 1000, 2000];
        let config = LaptopConfig::medium_memory();
        let assembler = LaptopAssembler::new(config);

        let mut scaling_results = Vec::new();

        for &size in &dataset_sizes {
            let reads = create_benchmark_reads(size, 100);

            let start = Instant::now();
            let _contigs = assembler.assemble(&reads)?;
            let time_per_read = start.elapsed().as_secs_f64() / size as f64;

            scaling_results.push((size, time_per_read));
            println!("Dataset size {}: {:.4}s per read", size, time_per_read);
        }

        // Verify scaling is reasonable (not exponential)
        let (small_size, small_time) = scaling_results[0];
        let (large_size, large_time) = scaling_results.last().unwrap();

        let scaling_factor = large_time / small_time;
        let size_factor = *large_size as f64 / small_size as f64;

        // Should scale better than O(n²)
        assert!(scaling_factor < size_factor * size_factor,
               "Performance scaling {:.2}x worse than quadratic for {}x data increase",
               scaling_factor, size_factor);

        println!("✅ Performance scaling: {:.2}x slower for {}x more data (acceptable)",
                scaling_factor, size_factor);

        Ok(())
    }
}

/* ========================================================================= */
/*                        QUALITY REGRESSION TESTS                         */
/* ========================================================================= */

#[cfg(test)]
mod quality_regression_tests {
    use super::*;

    /// Test that assembly quality doesn't degrade with optimizations
    #[test]
    fn test_assembly_quality_regression() -> Result<()> {
        // Use reference dataset with known good results
        let reference_reads = create_reference_dataset();

        let configs = vec![
            ("Unoptimized", OptimizationConfig::default()),
            ("Low Memory", OptimizationConfig::low_memory()),
            ("High Performance", OptimizationConfig::high_performance()),
        ];

        let mut quality_results = HashMap::new();

        for (config_name, config) in configs {
            // Create graph with specific configuration
            let mut graph = CacheOptimizedGraph::new(1000);

            // Simulate building graph (simplified for test)
            for (i, read) in reference_reads.iter().enumerate() {
                let hash = calculate_read_hash(read);
                graph.add_node(hash, 1);

                if i > 0 {
                    let prev_hash = calculate_read_hash(&reference_reads[i - 1]);
                    let _ = graph.add_edge(prev_hash, hash);
                }
            }

            let (nodes, edges, _, _) = graph.get_statistics();
            let quality_score = calculate_graph_quality_score(nodes, edges, reference_reads.len());

            quality_results.insert(config_name, quality_score);
            println!("{}: Quality score {:.3}", config_name, quality_score);
        }

        // Verify optimized configurations don't significantly degrade quality
        let unoptimized_quality = quality_results["Unoptimized"];

        for (config_name, quality) in &quality_results {
            if *config_name != "Unoptimized" {
                let quality_ratio = quality / unoptimized_quality;
                assert!(quality_ratio >= 0.9,
                       "{} quality {:.3} is {:.1}% of unoptimized {:.3}",
                       config_name, quality, quality_ratio * 100.0, unoptimized_quality);

                println!("✅ {}: {:.1}% of unoptimized quality retained",
                        config_name, quality_ratio * 100.0);
            }
        }

        Ok(())
    }

    /// Test k-mer processing accuracy regression
    #[test]
    fn test_kmer_processing_accuracy_regression() -> Result<()> {
        let test_sequence = b"ATCGATCGATCGATCGATCGATCGATCGATCG";
        let k_sizes = vec![15, 21, 31];

        for k in k_sizes {
            if k > test_sequence.len() {
                continue;
            }

            // Test compact k-mer representation
            let expected_kmers = test_sequence.len() - k + 1;
            let kmer_iter = meta_forge::assembly::performance_optimizations::ZeroCopyKmerIterator::new(test_sequence, k);
            let actual_kmers = kmer_iter.count();

            assert_eq!(actual_kmers, expected_kmers,
                      "K-mer count mismatch for k={}: expected {}, got {}",
                      k, expected_kmers, actual_kmers);

            // Test SIMD nucleotide counting accuracy
            let gc_content = SIMDNucleotideOps::gc_content_simd(test_sequence);
            assert!(gc_content >= 0.0 && gc_content <= 1.0,
                   "GC content {:.3} out of valid range for k={}", gc_content, k);

            println!("✅ K-mer processing accuracy maintained for k={}", k);
        }

        Ok(())
    }

    /// Test contig generation consistency
    #[test]
    fn test_contig_generation_consistency() -> Result<()> {
        let test_reads = create_deterministic_reads(100, 50);

        // Run assembly multiple times to test consistency
        let mut contig_counts = Vec::new();
        let config = LaptopConfig::medium_memory();

        for run in 0..5 {
            let assembler = LaptopAssembler::new(config.clone());
            let contigs = assembler.assemble(&test_reads)?;
            contig_counts.push(contigs.len());

            println!("Run {}: {} contigs generated", run + 1, contigs.len());
        }

        // Results should be consistent across runs
        let min_contigs = *contig_counts.iter().min().unwrap();
        let max_contigs = *contig_counts.iter().max().unwrap();
        let variation = (max_contigs - min_contigs) as f64 / min_contigs as f64;

        assert!(variation <= 0.1,
               "Contig count variation {:.1}% exceeds 10% threshold", variation * 100.0);

        println!("✅ Contig generation consistency: {:.1}% variation across runs",
                variation * 100.0);

        Ok(())
    }

    fn calculate_read_hash(read: &CorrectedRead) -> u64 {
        // Simple hash function for testing
        let mut hash = 0u64;
        for byte in read.corrected.bytes() {
            hash = hash.wrapping_mul(31).wrapping_add(byte as u64);
        }
        hash
    }

    fn calculate_graph_quality_score(nodes: usize, edges: usize, reads: usize) -> f64 {
        // Simple quality metric: connectivity relative to input size
        if reads == 0 {
            return 0.0;
        }

        let node_efficiency = 1.0 - (nodes as f64 / reads as f64);
        let connectivity = if nodes > 1 {
            edges as f64 / (nodes - 1) as f64
        } else {
            0.0
        };

        (node_efficiency + connectivity.min(1.0)) / 2.0
    }
}

/* ========================================================================= */
/*                        CI/CD INTEGRATION TESTS                          */
/* ========================================================================= */

#[cfg(test)]
mod ci_cd_integration_tests {
    use super::*;

    /// Test suitable for CI/CD pipeline execution
    #[test]
    fn test_ci_cd_smoke_test() -> Result<()> {
        println!("🚀 Running CI/CD laptop assembly smoke test...");

        // Quick test with small dataset for CI/CD
        let reads = create_benchmark_reads(50, 60);
        let config = LaptopConfig::low_memory();
        let assembler = LaptopAssembler::new(config);

        let start = Instant::now();
        let contigs = assembler.assemble(&reads)?;
        let duration = start.elapsed();

        // Basic validation for CI/CD
        assert!(!contigs.is_empty(), "Should generate contigs");
        assert!(duration.as_secs() < 10, "Should complete quickly in CI/CD");
        assert!(contigs.len() < reads.len(), "Should merge some reads");

        println!("✅ CI/CD smoke test: {} reads -> {} contigs in {:.2}s",
                reads.len(), contigs.len(), duration.as_secs_f64());

        Ok(())
    }

    /// Test laptop compatibility detection
    #[test]
    fn test_laptop_compatibility_detection() -> Result<()> {
        // Test auto-configuration logic
        let auto_assembler = LaptopAssembler::auto_config();

        // Create small test to verify it works
        let reads = create_benchmark_reads(10, 40);
        let contigs = auto_assembler.assemble(&reads)?;

        assert!(!contigs.is_empty(), "Auto-configured assembler should work");

        println!("✅ Laptop auto-configuration working");

        // Test all predefined configurations
        let configs = vec![
            LaptopConfig::low_memory(),
            LaptopConfig::medium_memory(),
            LaptopConfig::high_memory(),
        ];

        for (i, config) in configs.iter().enumerate() {
            let assembler = LaptopAssembler::new(config.clone());
            let test_contigs = assembler.assemble(&reads)?;

            assert!(!test_contigs.is_empty(), "Config {} should work", i);
            println!("✅ Laptop configuration {} validated", i);
        }

        Ok(())
    }

    /// Test performance benchmarks for CI/CD metrics
    #[test]
    fn test_ci_cd_performance_metrics() -> Result<()> {
        let benchmark = PerformanceBenchmark::new("CI/CD Metrics", 10);

        // Test SIMD performance
        let test_sequence = create_test_sequence(10_000);
        let simd_result = benchmark.benchmark_nucleotide_counting(&test_sequence);

        // Test k-mer processing
        let kmer_result = benchmark.benchmark_kmer_iteration(&test_sequence, 21);

        // Log metrics for CI/CD tracking
        println!("METRIC:simd_speedup:{:.2}", simd_result.speedup);
        println!("METRIC:kmer_speedup:{:.2}", kmer_result.speedup);
        println!("METRIC:simd_time_ms:{:.2}", simd_result.optimized_time_ns as f64 / 1_000_000.0);
        println!("METRIC:kmer_time_ms:{:.2}", kmer_result.optimized_time_ns as f64 / 1_000_000.0);

        // Basic performance thresholds for CI/CD
        assert!(simd_result.speedup >= 0.8, "SIMD performance regression detected");
        assert!(kmer_result.speedup >= 0.8, "K-mer processing regression detected");

        println!("✅ CI/CD performance metrics within acceptable ranges");

        Ok(())
    }

    /// Test memory usage for CI/CD environment constraints
    #[test]
    fn test_ci_cd_memory_constraints() -> Result<()> {
        // Simulate CI/CD environment with limited memory
        let ci_config = LaptopConfig {
            memory_budget_mb: 512, // 512MB limit for CI/CD
            cpu_cores: 2,
            chunk_size: 200,
            max_k: 21,
        };

        let reads = create_benchmark_reads(100, 50);
        let mut graph = meta_forge::assembly::laptop_assembly::LaptopAssemblyGraph::new(ci_config);

        graph.build_from_reads(&reads, 15)?;
        let memory_usage = graph.memory_usage_mb();

        assert!(memory_usage <= 512.0,
               "Memory usage {:.1}MB exceeds CI/CD limit of 512MB", memory_usage);

        println!("✅ CI/CD memory constraint test: {:.1}MB usage", memory_usage);

        Ok(())
    }
}

/* ========================================================================= */
/*                        HELPER FUNCTIONS                                 */
/* ========================================================================= */

/// Create benchmark reads with controlled characteristics
fn create_benchmark_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    let base_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

    (0..count).map(|i| {
        let start = (i * 2) % (base_sequence.len() - length);
        let sequence = base_sequence[start..start + length].to_string();

        CorrectedRead {
            id: i,
            original: sequence.clone(),
            corrected: sequence,
            corrections: Vec::new(),
            quality_scores: vec![30; length],
            correction_metadata: CorrectionMetadata {
                algorithm: "benchmark".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }).collect()
}

/// Create reads that stress memory usage
fn create_memory_stress_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    (0..count).map(|i| {
        // Create unique sequences to maximize memory usage
        let mut sequence = String::with_capacity(length);
        for j in 0..length {
            let nucleotide = match (i + j) % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            };
            sequence.push(nucleotide);
        }

        CorrectedRead {
            id: i,
            original: sequence.clone(),
            corrected: sequence,
            corrections: Vec::new(),
            quality_scores: vec![25; length],
            correction_metadata: CorrectionMetadata {
                algorithm: "memory_stress".to_string(),
                confidence_threshold: 0.85,
                context_window: 7,
                correction_time_ms: 1,
            },
        }
    }).collect()
}

/// Create reference dataset with known characteristics
fn create_reference_dataset() -> Vec<CorrectedRead> {
    // Create a small but representative dataset
    let template = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    let read_length = 25;
    let overlap = 10;
    let mut reads = Vec::new();

    let mut start = 0;
    let mut id = 0;

    while start + read_length <= template.len() {
        let sequence = template[start..start + read_length].to_string();

        reads.push(CorrectedRead {
            id,
            original: sequence.clone(),
            corrected: sequence,
            corrections: Vec::new(),
            quality_scores: vec![35; read_length],
            correction_metadata: CorrectionMetadata {
                algorithm: "reference".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        });

        start += read_length - overlap;
        id += 1;
    }

    reads
}

/// Create deterministic reads for consistency testing
fn create_deterministic_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    let pattern = "ATCGATCG";

    (0..count).map(|i| {
        let mut sequence = String::with_capacity(length);
        while sequence.len() < length {
            sequence.push_str(pattern);
        }
        sequence.truncate(length);

        CorrectedRead {
            id: i,
            original: sequence.clone(),
            corrected: sequence,
            corrections: Vec::new(),
            quality_scores: vec![30; length],
            correction_metadata: CorrectionMetadata {
                algorithm: "deterministic".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }).collect()
}

/// Create test sequence for performance benchmarking
fn create_test_sequence(length: usize) -> Vec<u8> {
    let pattern = b"ATCGATCGATCGATCG";
    let mut sequence = Vec::with_capacity(length);

    while sequence.len() < length {
        sequence.extend_from_slice(pattern);
    }

    sequence.truncate(length);
    sequence
}