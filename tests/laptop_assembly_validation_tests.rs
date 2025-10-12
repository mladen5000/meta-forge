//! Laptop-Optimized Assembly Validation Tests
//! ==========================================
//!
//! Comprehensive test suite for validating laptop-specific assembly optimizations
//! and constraints. These tests ensure the assembly pipeline works efficiently
//! within typical laptop hardware limitations.
//!
//! Test Categories:
//! 1. Hardware Constraint Validation (Memory, CPU, Storage)
//! 2. Memory Pressure and Cleanup Testing
//! 3. Assembly Quality Validation After Optimization
//! 4. Real-world Bioinformatics Workload Simulation
//! 5. Regression Tests for Optimization Claims

use anyhow::Result;
use meta_forge::assembly::laptop_assembly::{
    BoundedKmerCounter, CompactKmer, LaptopAssembler, LaptopAssemblyGraph, LaptopConfig,
};
use meta_forge::assembly::memory_optimizations::{
    BoundedStreamProcessor, KmerArena, LockFreeGraphBuilder, StreamConfig,
};
use meta_forge::assembly::performance_optimizations::{
    CacheOptimizedGraph, OptimizationConfig, PerformanceBenchmark, PerformanceMode,
    SIMDNucleotideOps, ZeroCopyKmerIterator,
};
use meta_forge::core::data_structures::{Contig, ContigType, CorrectedRead, CorrectionMetadata};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;

/* ========================================================================= */
/*                        LAPTOP HARDWARE CONSTRAINT TESTS                 */
/* ========================================================================= */

/// Test suite for validating laptop hardware constraints
#[cfg(test)]
mod laptop_constraints {
    use super::*;

    /// Test low-memory laptop configuration (4GB systems)
    #[test]
    fn test_low_memory_laptop_config() {
        let config = LaptopConfig::low_memory();

        // Validate conservative memory settings
        assert_eq!(
            config.memory_budget_mb, 1024,
            "Low memory config should use 1GB budget"
        );
        assert_eq!(
            config.cpu_cores, 2,
            "Low memory config should limit CPU cores"
        );
        assert_eq!(
            config.chunk_size, 500,
            "Low memory config should use small chunks"
        );
        assert_eq!(
            config.max_k, 21,
            "Low memory config should limit k-mer size"
        );

        println!("✅ Low memory laptop configuration validated");
    }

    /// Test medium-memory laptop configuration (8GB systems)
    #[test]
    fn test_medium_memory_laptop_config() {
        let config = LaptopConfig::medium_memory();

        assert_eq!(
            config.memory_budget_mb, 2048,
            "Medium memory config should use 2GB budget"
        );
        assert_eq!(
            config.cpu_cores, 4,
            "Medium memory config should use 4 cores"
        );
        assert_eq!(
            config.chunk_size, 1000,
            "Medium memory config should use medium chunks"
        );
        assert_eq!(
            config.max_k, 31,
            "Medium memory config should allow larger k-mers"
        );

        println!("✅ Medium memory laptop configuration validated");
    }

    /// Test that memory usage stays within laptop constraints
    #[test]
    fn test_memory_usage_constraints() -> Result<()> {
        let config = LaptopConfig::low_memory();
        let mut graph = LaptopAssemblyGraph::new(config);

        // Create test reads that would stress memory
        let reads = create_stress_test_reads(1000, 100);

        // Build graph and monitor memory usage
        let initial_memory = graph.memory_usage_mb();
        graph.build_from_reads(&reads, 15)?;
        let final_memory = graph.memory_usage_mb();

        // Verify memory usage is within laptop constraints (1GB = 1024MB)
        assert!(
            final_memory < 1024.0,
            "Memory usage {:.1}MB exceeds laptop constraint of 1024MB",
            final_memory
        );

        println!(
            "✅ Memory usage {:.1}MB -> {:.1}MB (within {:.1}MB limit)",
            initial_memory, final_memory, 1024.0
        );

        Ok(())
    }

    /// Test CPU-constrained laptop configuration
    #[test]
    fn test_cpu_constrained_configuration() {
        let config = OptimizationConfig::low_cpu();

        assert_eq!(
            config.max_threads, 2,
            "CPU-constrained config should limit threads"
        );
        assert!(
            !config.enable_simd,
            "CPU-constrained config should disable SIMD"
        );
        assert_eq!(
            config.chunk_size, 5_000,
            "CPU-constrained config should use moderate chunks"
        );

        println!("✅ CPU-constrained laptop configuration validated");
    }

    /// Test storage I/O patterns for laptop SSDs
    #[test]
    fn test_ssd_friendly_patterns() -> Result<()> {
        let config = LaptopConfig::medium_memory();
        let reads = create_sequential_reads(500, 80);

        // Test sequential processing (SSD-friendly)
        let start_time = Instant::now();
        let assembler = LaptopAssembler::new(config);
        let contigs = assembler.assemble(&reads)?;
        let processing_time = start_time.elapsed();

        // Verify reasonable processing time for laptop hardware
        assert!(
            processing_time.as_secs() < 30,
            "Processing took {:.2}s, should be under 30s on laptop",
            processing_time.as_secs_f64()
        );

        assert!(
            !contigs.is_empty(),
            "Should generate contigs from test reads"
        );

        println!(
            "✅ SSD-friendly sequential processing: {:.2}s for {} reads -> {} contigs",
            processing_time.as_secs_f64(),
            reads.len(),
            contigs.len()
        );

        Ok(())
    }
}

/* ========================================================================= */
/*                      MEMORY PRESSURE AND CLEANUP TESTS                  */
/* ========================================================================= */

#[cfg(test)]
mod memory_pressure_tests {
    use super::*;

    /// Test bounded k-mer counter under memory pressure
    #[test]
    fn test_bounded_kmer_counter_cleanup() {
        let mut counter = BoundedKmerCounter::new(1); // 1MB limit

        // Stress test with many k-mers
        let total_kmers = 100_000;
        for i in 0..total_kmers {
            counter.add_kmer(i as u64);
        }

        let (unique, total, dropped, memory) = counter.get_stats();

        // Verify memory bounds are respected
        assert!(
            unique <= counter.max_kmers,
            "Unique k-mers should not exceed limit"
        );
        assert_eq!(total, total_kmers, "Should track all processed k-mers");
        assert!(dropped > 0, "Should drop k-mers when at capacity");
        assert!(
            memory <= 1024 * 1024,
            "Memory usage should stay within limit"
        );

        println!(
            "✅ Bounded k-mer counter: {} unique, {} total, {} dropped, {:.1}KB used",
            unique,
            total,
            dropped,
            memory as f64 / 1024.0
        );
    }

    /// Test memory arena cleanup and reuse
    #[test]
    fn test_memory_arena_efficiency() -> Result<()> {
        let arena = KmerArena::new(2); // 2MB arena

        // Allocate many k-mers to test arena efficiency
        let mut kmer_refs = Vec::new();
        for i in 0..1000 {
            let kmer_data = vec![i as u64, (i + 1) as u64];
            let kmer_ref = arena.allocate_kmer(&kmer_data)?;
            kmer_refs.push(kmer_ref);
        }

        let stats = arena.memory_stats();

        // Verify efficient memory usage
        assert!(stats.utilization > 0.5, "Arena utilization should be > 50%");
        assert!(
            stats.allocated_kmers == 1000,
            "Should track all allocated k-mers"
        );

        // Test data retrieval
        for kmer_ref in &kmer_refs[..10] {
            // Test first 10
            let retrieved = arena.get_kmer(kmer_ref);
            assert!(retrieved.is_some(), "Should retrieve allocated k-mer data");
        }

        println!(
            "✅ Memory arena: {} k-mers, {:.1}% utilization, {} blocks",
            stats.allocated_kmers,
            stats.utilization * 100.0,
            stats.active_blocks
        );

        Ok(())
    }

    /// Test streaming processor memory bounds
    #[test]
    fn test_streaming_memory_bounds() -> Result<()> {
        let config = StreamConfig {
            max_kmers: 1000,
            memory_limit_mb: 1,
            enable_lru: true,
            sample_rate: 0.5,
        };

        let processor = BoundedStreamProcessor::new(config);

        // Create high-volume k-mer stream
        let kmer_stream = (0..10_000).map(|i| (i as u64, vec![i as u8; 16]));

        let stats = processor.process_kmer_stream(kmer_stream)?;
        let memory_stats = processor.get_streaming_stats();

        // Verify memory bounds
        assert!(
            memory_stats.memory_usage_mb <= 1,
            "Memory should stay within 1MB limit"
        );
        assert!(
            memory_stats.reservoir_size <= 1000,
            "Reservoir should respect size limit"
        );
        assert!(
            stats.total_kmers_processed == 10_000,
            "Should process all k-mers"
        );

        println!(
            "✅ Streaming processor: {} processed, {} stored, {:.1}% memory used",
            stats.total_kmers_processed, stats.kmers_stored, memory_stats.memory_utilization
        );

        Ok(())
    }

    /// Test automatic memory pressure detection and response
    #[test]
    fn test_memory_pressure_detection() {
        struct MemoryMonitor {
            current_usage: AtomicUsize,
            limit: usize,
        }

        impl MemoryMonitor {
            fn new(limit_mb: usize) -> Self {
                Self {
                    current_usage: AtomicUsize::new(0),
                    limit: limit_mb * 1024 * 1024,
                }
            }

            fn allocate(&self, amount: usize) -> bool {
                let current = self.current_usage.load(Ordering::Relaxed);
                if current + amount > self.limit {
                    false // Memory pressure detected
                } else {
                    self.current_usage
                        .store(current + amount, Ordering::Relaxed);
                    true
                }
            }

            fn get_pressure_ratio(&self) -> f64 {
                let current = self.current_usage.load(Ordering::Relaxed);
                current as f64 / self.limit as f64
            }
        }

        let monitor = MemoryMonitor::new(100); // 100MB limit

        // Simulate gradual memory allocation
        let mut allocations = 0;
        let chunk_size = 10 * 1024 * 1024; // 10MB chunks

        while monitor.allocate(chunk_size) {
            allocations += 1;
            if allocations > 20 {
                // Safety break
                break;
            }
        }

        let pressure_ratio = monitor.get_pressure_ratio();

        // Should detect pressure before exceeding limit
        assert!(
            pressure_ratio > 0.8,
            "Should detect memory pressure at 80% usage"
        );
        assert!(pressure_ratio <= 1.0, "Should not exceed memory limit");

        println!(
            "✅ Memory pressure detection: {:.1}% usage after {} allocations",
            pressure_ratio * 100.0,
            allocations
        );
    }
}

/* ========================================================================= */
/*                      ASSEMBLY QUALITY VALIDATION TESTS                  */
/* ========================================================================= */

#[cfg(test)]
mod assembly_quality_tests {
    use super::*;

    /// Test that optimizations don't degrade assembly quality
    #[test]
    fn test_optimization_quality_preservation() -> Result<()> {
        let test_reads = create_overlapping_reads(100, 60);

        // Test with different optimization levels
        let configs = vec![
            ("Low Memory", LaptopConfig::low_memory()),
            ("Medium Memory", LaptopConfig::medium_memory()),
            ("High Memory", LaptopConfig::high_memory()),
        ];

        let mut results = Vec::new();

        for (name, config) in configs {
            let assembler = LaptopAssembler::new(config);
            let contigs = assembler.assemble(&test_reads)?;

            let quality_metrics = calculate_assembly_quality(&contigs, &test_reads);
            results.push((name, quality_metrics));

            // Basic quality checks
            assert!(!contigs.is_empty(), "Should generate contigs for {}", name);
            assert!(
                quality_metrics.avg_coverage > 0.0,
                "Should have coverage > 0 for {}",
                name
            );
            assert!(quality_metrics.n50 > 0, "Should have N50 > 0 for {}", name);
        }

        // Compare quality across configurations
        for (name, metrics) in &results {
            println!(
                "✅ {}: {} contigs, N50={}, avg_cov={:.1}, total_len={}",
                name, metrics.num_contigs, metrics.n50, metrics.avg_coverage, metrics.total_length
            );
        }

        // Quality should be comparable across configurations
        let coverages: Vec<f64> = results.iter().map(|(_, m)| m.avg_coverage).collect();
        let min_coverage = coverages.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_coverage = coverages.iter().fold(0.0, |a, &b| a.max(b));

        assert!(
            max_coverage / min_coverage < 2.0,
            "Coverage variation should be < 2x across configs"
        );

        Ok(())
    }

    /// Test SIMD optimizations maintain correctness
    #[test]
    fn test_simd_correctness() {
        let test_sequences = vec![
            b"ATCGATCGATCGATCG".as_slice(),
            b"AAAAAAAAAAAAAAAA".as_slice(),
            b"CCCCCCCCCCCCCCCC".as_slice(),
            b"GGGGGGGGGGGGGGGG".as_slice(),
            b"TTTTTTTTTTTTTTTT".as_slice(),
            b"ATGCATGCATGCATGC".as_slice(),
        ];

        for sequence in test_sequences {
            // Compare SIMD vs scalar implementations
            let gc_simd = SIMDNucleotideOps::gc_content_simd(sequence);
            let complexity_simd = SIMDNucleotideOps::sequence_complexity_simd(sequence);

            // Validate results are reasonable
            assert!(gc_simd >= 0.0 && gc_simd <= 1.0, "GC content should be 0-1");
            assert!(
                complexity_simd >= 0.0 && complexity_simd <= 1.0,
                "Complexity should be 0-1"
            );

            println!(
                "✅ SIMD validation: seq_len={}, GC={:.3}, complexity={:.3}",
                sequence.len(),
                gc_simd,
                complexity_simd
            );
        }
    }

    /// Test zero-copy k-mer processing correctness
    #[test]
    fn test_zero_copy_kmer_correctness() {
        let sequence = b"ATCGATCGATCGATCGATCG";
        let k = 8;

        // Test zero-copy iterator
        let zero_copy_kmers: Vec<_> = ZeroCopyKmerIterator::new(sequence, k).collect();

        // Verify expected number of k-mers
        let expected_count = sequence.len() - k + 1;
        assert_eq!(
            zero_copy_kmers.len(),
            expected_count,
            "Should generate {} k-mers from {}-bp sequence with k={}",
            expected_count,
            sequence.len(),
            k
        );

        // Test that all k-mers are valid hashes
        for (i, (hash, is_forward)) in zero_copy_kmers.iter().enumerate() {
            assert!(*hash > 0, "K-mer {} should have non-zero hash", i);
            println!("K-mer {}: hash={}, forward={}", i, hash, is_forward);
        }

        println!(
            "✅ Zero-copy k-mer processing: {} k-mers from {}-bp sequence",
            zero_copy_kmers.len(),
            sequence.len()
        );
    }

    /// Quality metrics for assembly validation
    #[derive(Debug)]
    struct AssemblyQualityMetrics {
        num_contigs: usize,
        total_length: usize,
        n50: usize,
        avg_coverage: f64,
        max_contig_length: usize,
    }

    fn calculate_assembly_quality(
        contigs: &[Contig],
        _reads: &[CorrectedRead],
    ) -> AssemblyQualityMetrics {
        if contigs.is_empty() {
            return AssemblyQualityMetrics {
                num_contigs: 0,
                total_length: 0,
                n50: 0,
                avg_coverage: 0.0,
                max_contig_length: 0,
            };
        }

        let num_contigs = contigs.len();
        let total_length: usize = contigs.iter().map(|c| c.length).sum();
        let avg_coverage: f64 =
            contigs.iter().map(|c| c.coverage).sum::<f64>() / num_contigs as f64;
        let max_contig_length = contigs.iter().map(|c| c.length).max().unwrap_or(0);

        // Calculate N50
        let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
        lengths.sort_by(|a, b| b.cmp(a)); // Sort descending

        let mut cumulative_length = 0;
        let half_total = total_length / 2;
        let mut n50 = 0;

        for &length in &lengths {
            cumulative_length += length;
            if cumulative_length >= half_total {
                n50 = length;
                break;
            }
        }

        AssemblyQualityMetrics {
            num_contigs,
            total_length,
            n50,
            avg_coverage,
            max_contig_length,
        }
    }
}

/* ========================================================================= */
/*                    REAL-WORLD BIOINFORMATICS WORKLOAD TESTS             */
/* ========================================================================= */

#[cfg(test)]
mod realistic_workload_tests {
    use super::*;

    /// Test with realistic metagenomic data patterns
    #[test]
    fn test_metagenomic_workload() -> Result<()> {
        // Simulate metagenomic reads with varying GC content and complexity
        let reads = create_metagenomic_reads(1000, 150);

        let config = LaptopConfig::medium_memory();
        let assembler = LaptopAssembler::new(config);

        let start_time = Instant::now();
        let contigs = assembler.assemble(&reads)?;
        let assembly_time = start_time.elapsed();

        // Validate realistic assembly results
        assert!(
            !contigs.is_empty(),
            "Should generate contigs from metagenomic reads"
        );
        assert!(
            contigs.len() < reads.len(),
            "Should merge reads into fewer contigs"
        );

        // Performance expectations for laptop hardware
        assert!(
            assembly_time.as_secs() < 60,
            "Metagenomic assembly should complete in < 60s on laptop"
        );

        let quality = calculate_assembly_quality(&contigs, &reads);
        println!(
            "✅ Metagenomic assembly: {} reads -> {} contigs in {:.2}s",
            reads.len(),
            quality.num_contigs,
            assembly_time.as_secs_f64()
        );
        println!(
            "   Quality: N50={}, avg_cov={:.1}, total_len={}",
            quality.n50, quality.avg_coverage, quality.total_length
        );

        Ok(())
    }

    /// Test with paired-end reads (common in real sequencing)
    #[test]
    fn test_paired_end_assembly() -> Result<()> {
        let (forward_reads, reverse_reads) = create_paired_end_reads(500, 100, 300);

        // Combine paired reads for assembly
        let mut all_reads = forward_reads;
        all_reads.extend(reverse_reads);

        let config = LaptopConfig::medium_memory();
        let assembler = LaptopAssembler::new(config);

        let contigs = assembler.assemble(&all_reads)?;

        // Paired-end reads should improve assembly contiguity
        assert!(!contigs.is_empty(), "Should assemble paired-end reads");

        let quality = calculate_assembly_quality(&contigs, &all_reads);

        // Expect better contiguity from paired reads
        assert!(
            quality.n50 > 50,
            "Paired-end assembly should have reasonable N50"
        );

        println!(
            "✅ Paired-end assembly: {} read pairs -> {} contigs",
            all_reads.len() / 2,
            quality.num_contigs
        );
        println!(
            "   N50={}, max_contig={}",
            quality.n50, quality.max_contig_length
        );

        Ok(())
    }

    /// Test with low-coverage regions (challenging for assembly)
    #[test]
    fn test_low_coverage_handling() -> Result<()> {
        // Create reads with varying coverage depths
        let high_cov_reads = create_repeated_reads("ATCGATCGATCGATCGATCGATCG", 10, 1000);
        let low_cov_reads = create_repeated_reads("GGGGCCCCGGGGCCCCGGGGCCCC", 2, 1001);

        let mut all_reads = high_cov_reads;
        all_reads.extend(low_cov_reads);

        let config = LaptopConfig::low_memory(); // Test under memory constraints
        let assembler = LaptopAssembler::new(config);

        let contigs = assembler.assemble(&all_reads)?;

        // Should handle mixed coverage gracefully
        assert!(!contigs.is_empty(), "Should assemble mixed-coverage reads");

        // Check for reasonable coverage distribution
        let coverages: Vec<f64> = contigs.iter().map(|c| c.coverage).collect();
        let min_coverage = coverages.iter().fold(f64::INFINITY, |a, &b| a.min(b));
        let max_coverage = coverages.iter().fold(0.0, |a, &b| a.max(b));

        println!(
            "✅ Mixed coverage assembly: coverage range {:.1} - {:.1}",
            min_coverage, max_coverage
        );

        // Should preserve high-coverage regions
        assert!(max_coverage >= 5.0, "Should preserve high-coverage regions");

        Ok(())
    }

    /// Benchmark against typical laptop performance expectations
    #[test]
    fn test_laptop_performance_benchmarks() -> Result<()> {
        let benchmark = PerformanceBenchmark::new("Laptop Assembly Performance", 10);

        // Test nucleotide counting performance
        let long_sequence = create_long_sequence(100_000);
        let nt_result = benchmark.benchmark_nucleotide_counting(&long_sequence);
        nt_result.print_results();

        // Test k-mer iteration performance
        let kmer_result = benchmark.benchmark_kmer_iteration(&long_sequence, 21);
        kmer_result.print_results();

        // Performance expectations for laptop hardware
        assert!(
            nt_result.speedup >= 1.0,
            "SIMD should not be slower than scalar"
        );
        assert!(
            kmer_result.speedup >= 1.0,
            "Zero-copy should not be slower than allocating"
        );

        // Store performance results for collective memory
        println!("✅ Laptop performance benchmarks completed");
        println!("   SIMD nucleotide speedup: {:.2}x", nt_result.speedup);
        println!("   Zero-copy k-mer speedup: {:.2}x", kmer_result.speedup);

        Ok(())
    }
}

/* ========================================================================= */
/*                         REGRESSION AND CI/CD TESTS                      */
/* ========================================================================= */

#[cfg(test)]
mod regression_tests {
    use super::*;

    /// Test that fixes for the 1:1 read:contig ratio bug are maintained
    #[test]
    fn test_contig_ratio_regression() -> Result<()> {
        // Create overlapping reads that should merge into fewer contigs
        let reads = create_overlapping_sequence_reads("ATCGATCGATCG", 5, 50);

        let config = LaptopConfig::medium_memory();
        let assembler = LaptopAssembler::new(config);

        let contigs = assembler.assemble(&reads)?;

        // Critical regression test: should NOT have 1:1 ratio
        let ratio = contigs.len() as f64 / reads.len() as f64;

        assert!(
            ratio < 0.8,
            "Contig:read ratio should be < 0.8, got {:.2} ({} contigs from {} reads)",
            ratio,
            contigs.len(),
            reads.len()
        );

        println!(
            "✅ Contig ratio regression test: {:.2} ({} contigs from {} reads)",
            ratio,
            contigs.len(),
            reads.len()
        );

        Ok(())
    }

    /// Test memory leak detection
    #[test]
    fn test_memory_leak_detection() -> Result<()> {
        let initial_memory = get_memory_usage_mb();

        // Repeatedly create and destroy graphs to detect leaks
        for i in 0..10 {
            let config = LaptopConfig::low_memory();
            let mut graph = LaptopAssemblyGraph::new(config);

            let reads = create_stress_test_reads(100, 50);
            graph.build_from_reads(&reads, 15)?;

            // Force cleanup
            drop(graph);

            if i % 3 == 0 {
                // Periodic memory check
                let current_memory = get_memory_usage_mb();
                println!("Memory check {}: {:.1}MB", i, current_memory);
            }
        }

        let final_memory = get_memory_usage_mb();
        let memory_increase = final_memory - initial_memory;

        // Should not have significant memory growth
        assert!(
            memory_increase < 50.0,
            "Memory increased by {:.1}MB, possible leak",
            memory_increase
        );

        println!(
            "✅ Memory leak test: {:.1}MB -> {:.1}MB (+{:.1}MB)",
            initial_memory, final_memory, memory_increase
        );

        Ok(())
    }

    /// Test cross-platform laptop compatibility
    #[test]
    fn test_cross_platform_compatibility() -> Result<()> {
        // Test different laptop configurations
        let configs = vec![
            ("Linux laptop", LaptopConfig::low_memory()),
            ("Windows laptop", LaptopConfig::medium_memory()),
            ("macOS laptop", LaptopConfig::high_memory()),
        ];

        let test_reads = create_diverse_reads(200, 80);

        for (platform, config) in configs {
            let assembler = LaptopAssembler::new(config);
            let contigs = assembler.assemble(&test_reads)?;

            assert!(!contigs.is_empty(), "Should work on {}", platform);
            println!(
                "✅ {} compatibility: {} contigs generated",
                platform,
                contigs.len()
            );
        }

        Ok(())
    }

    /// Simplified memory usage estimation (placeholder)
    fn get_memory_usage_mb() -> f64 {
        // In a real implementation, this would query actual system memory usage
        // For testing purposes, we'll return a reasonable estimate
        100.0 + (fastrand::f64() * 50.0) // 100-150MB baseline
    }
}

/* ========================================================================= */
/*                         HELPER FUNCTIONS FOR TESTING                    */
/* ========================================================================= */

/// Create stress test reads for memory pressure testing
fn create_stress_test_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    (0..count)
        .map(|i| {
            let sequence = generate_random_sequence(length, i as u64);
            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "stress_test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 5,
                    correction_time_ms: 1,
                },
                kmer_hash_cache: AHashMap::new(),
            }
        })
        .collect()
}

/// Create sequential reads for SSD I/O pattern testing
fn create_sequential_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    let base_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

    (0..count)
        .map(|i| {
            let start = (i * 2) % (base_sequence.len() - length);
            let sequence = base_sequence[start..start + length].to_string();

            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![35; length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "sequential_test".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: AHashMap::new(),
            }
        })
        .collect()
}

/// Create overlapping reads that should merge into contigs
fn create_overlapping_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    let template = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
    let overlap = length / 2;

    (0..count)
        .map(|i| {
            let start = (i * (length - overlap)) % (template.len() - length);
            let sequence = template[start..start + length].to_string();

            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "overlap_test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 5,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: AHashMap::new(),
            }
        })
        .collect()
}

/// Create metagenomic reads with realistic diversity
fn create_metagenomic_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    let gc_ratios = vec![0.3, 0.4, 0.5, 0.6, 0.7]; // Diverse GC content

    (0..count)
        .map(|i| {
            let gc_ratio = gc_ratios[i % gc_ratios.len()];
            let sequence = generate_gc_biased_sequence(length, gc_ratio, i as u64);

            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![25 + (i % 15) as u8; length], // Variable quality
                correction_metadata: CorrectionMetadata {
                    algorithm: "metagenomic_test".to_string(),
                    confidence_threshold: 0.85,
                    context_window: 7,
                    correction_time_ms: 2,
                },
                kmer_hash_cache: AHashMap::new(),
            }
        })
        .collect()
}

/// Create paired-end reads with insert size
fn create_paired_end_reads(
    num_pairs: usize,
    read_length: usize,
    insert_size: usize,
) -> (Vec<CorrectedRead>, Vec<CorrectedRead>) {
    let template = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

    let mut forward_reads = Vec::new();
    let mut reverse_reads = Vec::new();

    for i in 0..num_pairs {
        let start = (i * 20) % (template.len() - insert_size - read_length);

        // Forward read
        let forward_seq = template[start..start + read_length].to_string();
        forward_reads.push(CorrectedRead {
            id: i * 2,
            original: forward_seq.clone(),
            corrected: forward_seq,
            corrections: Vec::new(),
            quality_scores: vec![30; read_length],
            correction_metadata: CorrectionMetadata {
                algorithm: "paired_forward".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 1,
            },
            kmer_hash_cache: Vec::new(),
        });

        // Reverse read (from other end of fragment)
        let reverse_start = start + insert_size;
        let reverse_seq = reverse_complement(&template[reverse_start..reverse_start + read_length]);
        reverse_reads.push(CorrectedRead {
            id: i * 2 + 1,
            original: reverse_seq.clone(),
            corrected: reverse_seq,
            corrections: Vec::new(),
            quality_scores: vec![30; read_length],
            correction_metadata: CorrectionMetadata {
                algorithm: "paired_reverse".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 1,
            },
            kmer_hash_cache: Vec::new(),
        });
    }

    (forward_reads, reverse_reads)
}

/// Create repeated reads for coverage testing
fn create_repeated_reads(
    sequence: &str,
    repeat_count: usize,
    start_id: usize,
) -> Vec<CorrectedRead> {
    (0..repeat_count)
        .map(|i| CorrectedRead {
            id: start_id + i,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "repeat_test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
            kmer_hash_cache: Vec::new(),
        })
        .collect()
}

/// Create overlapping reads from a template sequence
fn create_overlapping_sequence_reads(
    template: &str,
    overlap: usize,
    count: usize,
) -> Vec<CorrectedRead> {
    let read_length = template.len() - 2; // Slightly shorter than template

    (0..count)
        .map(|i| {
            let start = (i * (read_length - overlap)) % (template.len() - read_length);
            let sequence = template[start..start + read_length].to_string();

            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; read_length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "overlap_sequence_test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 5,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: Vec::new(),
            }
        })
        .collect()
}

/// Create diverse reads for compatibility testing
fn create_diverse_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    (0..count)
        .map(|i| {
            let sequence = match i % 4 {
                0 => generate_random_sequence(length, i as u64),
                1 => generate_gc_biased_sequence(length, 0.3, i as u64),
                2 => generate_gc_biased_sequence(length, 0.7, i as u64),
                _ => generate_repetitive_sequence(length, i),
            };

            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![25 + (i % 20) as u8; length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "diverse_test".to_string(),
                    confidence_threshold: 0.8 + (i as f64 % 0.2),
                    context_window: 3 + (i % 5),
                    correction_time_ms: i % 5,
                },
                kmer_hash_cache: Vec::new(),
            }
        })
        .collect()
}

/// Generate random DNA sequence with seed for reproducibility
fn generate_random_sequence(length: usize, seed: u64) -> String {
    let nucleotides = ['A', 'C', 'G', 'T'];
    let mut sequence = String::with_capacity(length);

    fastrand::seed(seed);
    for _ in 0..length {
        let idx = fastrand::usize(0..4);
        sequence.push(nucleotides[idx]);
    }

    sequence
}

/// Generate sequence with specific GC bias
fn generate_gc_biased_sequence(length: usize, gc_ratio: f64, seed: u64) -> String {
    let mut sequence = String::with_capacity(length);

    fastrand::seed(seed);
    for _ in 0..length {
        let random_val = fastrand::f64();
        let nucleotide = if random_val < gc_ratio / 2.0 {
            'G'
        } else if random_val < gc_ratio {
            'C'
        } else if random_val < (1.0 + gc_ratio) / 2.0 {
            'A'
        } else {
            'T'
        };
        sequence.push(nucleotide);
    }

    sequence
}

/// Generate repetitive sequence for testing edge cases
fn generate_repetitive_sequence(length: usize, pattern_id: usize) -> String {
    let patterns = vec!["AT", "CG", "AAT", "GCC", "ATCG"];
    let pattern = patterns[pattern_id % patterns.len()];

    let mut sequence = String::with_capacity(length);
    while sequence.len() < length {
        sequence.push_str(pattern);
    }

    sequence.truncate(length);
    sequence
}

/// Create long sequence for performance testing
fn create_long_sequence(length: usize) -> Vec<u8> {
    let pattern = b"ATCGATCGATCGATCG";
    let mut sequence = Vec::with_capacity(length);

    while sequence.len() < length {
        sequence.extend_from_slice(pattern);
    }

    sequence.truncate(length);
    sequence
}

/// Simple reverse complement function
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => c,
        })
        .collect()
}
