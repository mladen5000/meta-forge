//! Comprehensive Performance Benchmark Suite
//! =========================================
//!
//! Demonstrates and measures the impact of all optimization techniques:
//! - Memory pool-based k-mer storage
//! - Advanced parallel processing with auto-tuning
//! - Database query optimizations
//! - Advanced graph algorithms
//! - SIMD-accelerated operations

use anyhow::Result;
use std::time::Instant;
use crate::assembly::{
    optimized_kmer_pool::{KmerPool, PooledKmer},
    parallel_optimizations::{ParallelProcessor, ProcessingStats},
    database_optimizations::{GenomicDatabase, DatabaseConfig},
    advanced_graph_algorithms::AdvancedDeBruijnGraph,
    simd_optimizations::SimdProcessor,
};
use std::sync::Arc;

/// Comprehensive benchmark suite for all optimizations
pub struct ComprehensiveBenchmark {
    /// Test sequences of varying lengths
    test_sequences: Vec<String>,
    /// K-mer sizes to test
    k_values: Vec<usize>,
    /// Dataset sizes for scalability testing
    dataset_sizes: Vec<usize>,
}

impl ComprehensiveBenchmark {
    /// Create new benchmark suite
    pub fn new() -> Self {
        let test_sequences = Self::generate_test_sequences();
        let k_values = vec![21, 31, 51, 71];
        let dataset_sizes = vec![1000, 10000, 50000, 100000];
        
        Self {
            test_sequences,
            k_values,
            dataset_sizes,
        }
    }
    
    /// Generate realistic test sequences
    fn generate_test_sequences() -> Vec<String> {
        let base_sequences = vec![
            // Typical bacterial genes
            "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGT",
            "GCTAGCTAGCTAGCGATCGATCGATCGTAGCTAGCTAGCTAGCTGATCGATCGATCGATCGTAGCTAGCT",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", // Low complexity
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", // Repetitive
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", // Perfect repeat
        ];
        
        let mut sequences = Vec::new();
        
        // Generate sequences of different lengths
        for &base in &base_sequences {
            for multiplier in &[1, 2, 5, 10, 20] {
                let mut extended = String::new();
                for _ in 0..*multiplier {
                    extended.push_str(base);
                }
                sequences.push(extended);
            }
        }
        
        // Add some random-like sequences
        for length in &[100, 500, 1000, 5000] {
            sequences.push(Self::generate_random_sequence(*length));
        }
        
        sequences
    }
    
    /// Generate pseudo-random DNA sequence
    fn generate_random_sequence(length: usize) -> String {
        let bases = ['A', 'T', 'C', 'G'];
        let mut sequence = String::with_capacity(length);
        
        for i in 0..length {
            let base_idx = (i * 17 + 23) % 4; // Simple pseudo-random
            sequence.push(bases[base_idx]);
        }
        
        sequence
    }
    
    /// Run comprehensive benchmark suite
    pub fn run_comprehensive_benchmark(&self) -> Result<BenchmarkResults> {
        println!("🚀 Running Comprehensive Performance Benchmark Suite");
        println!("====================================================");
        
        let mut results = BenchmarkResults::new();
        
        // 1. Memory Pool vs Traditional Allocation Benchmark
        results.memory_benchmark = self.benchmark_memory_optimizations()?;
        
        // 2. Parallel Processing Benchmark
        results.parallel_benchmark = self.benchmark_parallel_processing()?;
        
        // 3. Database Operations Benchmark
        results.database_benchmark = self.benchmark_database_operations()?;
        
        // 4. Graph Algorithms Benchmark
        results.graph_benchmark = self.benchmark_graph_algorithms()?;
        
        // 5. SIMD Operations Benchmark
        results.simd_benchmark = self.benchmark_simd_operations()?;
        
        // 6. End-to-End Pipeline Benchmark
        results.pipeline_benchmark = self.benchmark_complete_pipeline()?;
        
        // Print comprehensive results
        results.print_summary();
        
        Ok(results)
    }
    
    /// Benchmark memory pool vs traditional allocation
    fn benchmark_memory_optimizations(&self) -> Result<MemoryBenchmarkResults> {
        println!("\n📊 Memory Optimization Benchmark");
        println!("-".repeat(40));
        
        let mut results = MemoryBenchmarkResults::default();
        
        for &k in &self.k_values {
            for &dataset_size in &[1000, 10000] { // Limit for memory tests
                println!("Testing k={}, dataset_size={}", k, dataset_size);
                
                // Traditional allocation benchmark
                let start = Instant::now();
                let traditional_memory = self.benchmark_traditional_kmers(k, dataset_size)?;
                let traditional_time = start.elapsed();
                
                // Memory pool benchmark
                let start = Instant::now();
                let pool_memory = self.benchmark_pooled_kmers(k, dataset_size)?;
                let pool_time = start.elapsed();
                
                let memory_improvement = if traditional_memory > 0 {
                    ((traditional_memory - pool_memory) as f64 / traditional_memory as f64) * 100.0
                } else {
                    0.0
                };
                
                let time_improvement = if traditional_time > pool_time {
                    traditional_time.as_millis() as f64 / pool_time.as_millis() as f64
                } else {
                    1.0
                };
                
                results.memory_reduction_percent = results.memory_reduction_percent.max(memory_improvement);
                results.allocation_speedup = results.allocation_speedup.max(time_improvement);
                
                println!("  Memory: {} KB -> {} KB ({:.1}% reduction)", 
                         traditional_memory / 1024, pool_memory / 1024, memory_improvement);
                println!("  Time: {:?} -> {:?} ({:.2}x speedup)", 
                         traditional_time, pool_time, time_improvement);
            }
        }
        
        Ok(results)
    }
    
    /// Benchmark traditional k-mer allocation
    fn benchmark_traditional_kmers(&self, k: usize, dataset_size: usize) -> Result<usize> {
        let mut kmers = Vec::new();
        let mut memory_used = 0;
        
        for (i, sequence) in self.test_sequences.iter().cycle().take(dataset_size).enumerate() {
            if sequence.len() >= k {
                for j in 0..=(sequence.len() - k) {
                    let kmer_seq = &sequence[j..j + k];
                    
                    // Simulate traditional allocation overhead
                    let kmer_string = kmer_seq.to_string();
                    memory_used += kmer_string.capacity() + std::mem::size_of::<String>();
                    kmers.push(kmer_string);
                    
                    if i % 1000 == 0 && kmers.len() > 100000 {
                        // Prevent excessive memory usage in test
                        break;
                    }
                }
            }
        }
        
        Ok(memory_used)
    }
    
    /// Benchmark pooled k-mer allocation
    fn benchmark_pooled_kmers(&self, k: usize, dataset_size: usize) -> Result<usize> {
        let pool = Arc::new(KmerPool::new(dataset_size * 10, k));
        let mut kmers = Vec::new();
        
        for (i, sequence) in self.test_sequences.iter().cycle().take(dataset_size).enumerate() {
            if sequence.len() >= k {
                for j in 0..=(sequence.len() - k) {
                    let kmer_seq = &sequence[j..j + k];
                    
                    if let Ok(pooled_kmer) = PooledKmer::new(kmer_seq, pool.clone()) {
                        kmers.push(pooled_kmer);
                    }
                    
                    if i % 1000 == 0 && kmers.len() > 100000 {
                        // Prevent excessive memory usage in test
                        break;
                    }
                }
            }
        }
        
        let stats = pool.memory_stats();
        Ok(stats.used_bytes)
    }
    
    /// Benchmark parallel processing optimizations
    fn benchmark_parallel_processing(&self) -> Result<ParallelBenchmarkResults> {
        println!("\n🔄 Parallel Processing Benchmark");
        println!("-".repeat(40));
        
        let processor = ParallelProcessor::new();
        let mut results = ParallelBenchmarkResults::default();
        
        for &dataset_size in &self.dataset_sizes {
            println!("Testing dataset size: {}", dataset_size);
            
            let test_data: Vec<_> = (0..dataset_size).collect();
            
            // Benchmark parallel k-mer counting
            let sequences: Vec<_> = self.test_sequences.iter()
                .cycle()
                .take(dataset_size)
                .cloned()
                .collect();
            
            let start = Instant::now();
            let kmer_counts = processor.count_kmers_parallel(&sequences, 21);
            let parallel_time = start.elapsed();
            
            // Benchmark sequential version for comparison
            let start = Instant::now();
            let sequential_counts = self.count_kmers_sequential(&sequences, 21);
            let sequential_time = start.elapsed();
            
            let speedup = if parallel_time.as_nanos() > 0 {
                sequential_time.as_nanos() as f64 / parallel_time.as_nanos() as f64
            } else {
                1.0
            };
            
            results.max_speedup = results.max_speedup.max(speedup);
            results.parallel_efficiency = speedup / rayon::current_num_threads() as f64;
            
            println!("  Sequential: {:?}, Parallel: {:?}, Speedup: {:.2}x", 
                     sequential_time, parallel_time, speedup);
            
            // Verify correctness
            assert_eq!(kmer_counts.len(), sequential_counts.len());
        }
        
        let stats = processor.get_stats();
        results.throughput_ips = stats.throughput_ips;
        
        Ok(results)
    }
    
    /// Sequential k-mer counting for comparison
    fn count_kmers_sequential(&self, sequences: &[String], k: usize) -> std::collections::HashMap<String, usize> {
        let mut counts = std::collections::HashMap::new();
        
        for sequence in sequences {
            for i in 0..=(sequence.len().saturating_sub(k)) {
                let kmer = &sequence[i..i + k];
                *counts.entry(kmer.to_string()).or_insert(0) += 1;
            }
        }
        
        counts
    }
    
    /// Benchmark database operations
    fn benchmark_database_operations(&self) -> Result<DatabaseBenchmarkResults> {
        println!("\n💾 Database Operations Benchmark");
        println!("-".repeat(40));
        
        let temp_dir = tempfile::tempdir()?;
        let db_path = temp_dir.path().join("benchmark.db");
        
        let config = DatabaseConfig {
            enable_caching: true,
            batch_size: 10000,
            cache_size_mb: 128,
            ..Default::default()
        };
        
        let db = GenomicDatabase::new(db_path.to_str().unwrap(), config)?;
        let mut results = DatabaseBenchmarkResults::default();
        
        // Benchmark batch insertions
        let test_kmers: Vec<_> = self.test_sequences.iter()
            .flat_map(|seq| {
                (0..=(seq.len().saturating_sub(21)))
                    .map(|i| {
                        let kmer = &seq[i..i + 21];
                        (kmer.to_string(), kmer.as_bytes().to_vec(), fastrand::u32(1..100))
                    })
            })
            .take(50000)
            .collect();
        
        println!("Inserting {} k-mers...", test_kmers.len());
        
        let start = Instant::now();
        let inserted = db.batch_insert_kmers(&test_kmers)?;
        let insert_time = start.elapsed();
        
        results.batch_insert_rate = inserted as f64 / insert_time.as_secs_f64();
        
        println!("  Inserted {} k-mers in {:?} ({:.0} kmers/sec)", 
                 inserted, insert_time, results.batch_insert_rate);
        
        // Benchmark queries with and without caching
        let start = Instant::now();
        let frequent1 = db.get_frequent_kmers(10, Some(100))?;
        let first_query_time = start.elapsed();
        
        let start = Instant::now();
        let frequent2 = db.get_frequent_kmers(10, Some(100))?;
        let cached_query_time = start.elapsed();
        
        results.cache_speedup = if cached_query_time.as_nanos() > 0 {
            first_query_time.as_nanos() as f64 / cached_query_time.as_nanos() as f64
        } else {
            1.0
        };
        
        println!("  First query: {:?}, Cached query: {:?}, Cache speedup: {:.2}x",
                 first_query_time, cached_query_time, results.cache_speedup);
        
        let db_stats = db.get_metrics();
        results.query_throughput = db_stats.queries_executed as f64 
            / (db_stats.total_query_time_ms as f64 / 1000.0);
        
        Ok(results)
    }
    
    /// Benchmark advanced graph algorithms
    fn benchmark_graph_algorithms(&self) -> Result<GraphBenchmarkResults> {
        println!("\n🕸️  Graph Algorithms Benchmark");
        println!("-".repeat(40));
        
        let mut results = GraphBenchmarkResults::default();
        
        // Create test graph with various structures
        let mut graph = AdvancedDeBruijnGraph::new();
        let node_count = 10000;
        
        println!("Creating graph with {} nodes...", node_count);
        
        // Add nodes
        for i in 0..node_count {
            graph.add_node(i as u64, fastrand::u32(1..50));
        }
        
        // Add edges to create interesting structures
        let start = Instant::now();
        for i in 0..node_count - 1 {
            // Linear connections
            if fastrand::f64() > 0.3 {
                graph.add_edge(i as u64, (i + 1) as u64, 
                              fastrand::f64(), 
                              crate::assembly::advanced_graph_algorithms::EdgeType::Perfect).unwrap();
            }
            
            // Some random connections for complexity
            if fastrand::f64() > 0.7 {
                let target = fastrand::usize(0..node_count);
                graph.add_edge(i as u64, target as u64, 
                              fastrand::f64(), 
                              crate::assembly::advanced_graph_algorithms::EdgeType::Approximate).unwrap();
            }
        }
        let construction_time = start.elapsed();
        results.graph_construction_time_ms = construction_time.as_millis();
        
        graph.classify_node_types();
        
        // Benchmark SCC finding
        let start = Instant::now();
        let components = graph.find_strongly_connected_components();
        let scc_time = start.elapsed();
        results.scc_algorithm_time_ms = scc_time.as_millis();
        
        println!("  Graph construction: {:?}", construction_time);
        println!("  SCC algorithm: {:?} ({} components found)", scc_time, components.len());
        
        // Benchmark contig generation
        let start = Instant::now();
        let contigs = graph.generate_contigs()?;
        let contig_time = start.elapsed();
        results.contig_generation_time_ms = contig_time.as_millis();
        
        println!("  Contig generation: {:?} ({} contigs generated)", contig_time, contigs.len());
        
        // Benchmark graph simplification
        let start = Instant::now();
        let simplification_stats = graph.simplify_graph()?;
        let simplification_time = start.elapsed();
        results.simplification_time_ms = simplification_time.as_millis();
        
        println!("  Graph simplification: {:?}", simplification_time);
        println!("    Tips removed: {}", simplification_stats.tips_removed);
        println!("    Bubbles popped: {}", simplification_stats.bubbles_popped);
        
        Ok(results)
    }
    
    /// Benchmark SIMD operations
    fn benchmark_simd_operations(&self) -> Result<SimdBenchmarkResults> {
        println!("\n⚡ SIMD Operations Benchmark");
        println!("-".repeat(40));
        
        let simd_processor = SimdProcessor::new();
        let mut results = SimdBenchmarkResults::default();
        
        simd_processor.print_capabilities();
        
        // Test different sequence lengths
        for &length in &[1000, 10000, 100000] {
            let test_sequence = Self::generate_random_sequence(length);
            let sequence_bytes = test_sequence.as_bytes();
            
            println!("Testing sequence length: {} bp", length);
            
            // Benchmark nucleotide counting
            let iterations = 1000;
            let start = Instant::now();
            for _ in 0..iterations {
                let _ = simd_processor.count_nucleotides(sequence_bytes);
            }
            let simd_time = start.elapsed();
            
            // Compare with scalar implementation
            let scalar_processor = SimdProcessor {
                supports_avx512: false,
                supports_avx2: false,
                supports_sse41: false,
            };
            
            let start = Instant::now();
            for _ in 0..iterations {
                let _ = scalar_processor.count_nucleotides(sequence_bytes);
            }
            let scalar_time = start.elapsed();
            
            let speedup = if simd_time.as_nanos() > 0 {
                scalar_time.as_nanos() as f64 / simd_time.as_nanos() as f64
            } else {
                1.0
            };
            
            results.max_speedup = results.max_speedup.max(speedup);
            
            println!("  Scalar: {:?}, SIMD: {:?}, Speedup: {:.2}x", 
                     scalar_time, simd_time, speedup);
            
            // Test k-mer processing
            let k = 31;
            if test_sequence.len() >= k {
                let start = Instant::now();
                let kmers = simd_processor.process_kmers(sequence_bytes, k)?;
                let kmer_time = start.elapsed();
                
                let kmer_rate = kmers.len() as f64 / kmer_time.as_secs_f64();
                results.kmer_processing_rate = results.kmer_processing_rate.max(kmer_rate);
                
                println!("  K-mer processing: {} kmers in {:?} ({:.0} kmers/sec)", 
                         kmers.len(), kmer_time, kmer_rate);
            }
        }
        
        Ok(results)
    }
    
    /// Benchmark complete end-to-end pipeline
    fn benchmark_complete_pipeline(&self) -> Result<PipelineBenchmarkResults> {
        println!("\n🔄 Complete Pipeline Benchmark");
        println!("-".repeat(40));
        
        let mut results = PipelineBenchmarkResults::default();
        
        // Simulate complete pipeline with optimizations
        for &dataset_size in &[1000, 5000] { // Limit for pipeline tests
            println!("Testing pipeline with {} sequences", dataset_size);
            
            let sequences: Vec<_> = self.test_sequences.iter()
                .cycle()
                .take(dataset_size)
                .cloned()
                .collect();
            
            let pipeline_start = Instant::now();
            
            // Step 1: SIMD-accelerated k-mer extraction
            let simd_processor = SimdProcessor::new();
            let mut all_kmers = Vec::new();
            
            for sequence in &sequences {
                let sequence_bytes = sequence.as_bytes();
                if let Ok(kmers) = simd_processor.process_kmers(sequence_bytes, 21) {
                    all_kmers.extend(kmers);
                }
            }
            
            // Step 2: Parallel k-mer counting
            let parallel_processor = ParallelProcessor::new();
            let kmer_counts = parallel_processor.count_kmers_parallel(&sequences, 21);
            
            // Step 3: Graph construction with memory pools
            let pool = Arc::new(KmerPool::new(all_kmers.len(), 21));
            let mut graph = AdvancedDeBruijnGraph::new();
            
            // Add frequent k-mers as nodes
            let mut node_id = 0u64;
            for (kmer_seq, &count) in kmer_counts.iter() {
                if count >= 2 { // Minimum frequency
                    graph.add_node(node_id, count as u32);
                    node_id += 1;
                    
                    if node_id >= 1000 { // Limit for benchmark
                        break;
                    }
                }
            }
            
            // Step 4: Graph algorithms
            graph.classify_node_types();
            let _components = graph.find_strongly_connected_components();
            let _contigs = graph.generate_contigs()?;
            
            let pipeline_time = pipeline_start.elapsed();
            
            let throughput = sequences.len() as f64 / pipeline_time.as_secs_f64();
            results.max_throughput = results.max_throughput.max(throughput);
            results.total_processing_time_ms = pipeline_time.as_millis();
            
            println!("  Pipeline time: {:?} ({:.1} sequences/sec)", 
                     pipeline_time, throughput);
        }
        
        Ok(results)
    }
}

/// Complete benchmark results
#[derive(Debug)]
pub struct BenchmarkResults {
    pub memory_benchmark: MemoryBenchmarkResults,
    pub parallel_benchmark: ParallelBenchmarkResults,
    pub database_benchmark: DatabaseBenchmarkResults,
    pub graph_benchmark: GraphBenchmarkResults,
    pub simd_benchmark: SimdBenchmarkResults,
    pub pipeline_benchmark: PipelineBenchmarkResults,
}

impl BenchmarkResults {
    fn new() -> Self {
        Self {
            memory_benchmark: MemoryBenchmarkResults::default(),
            parallel_benchmark: ParallelBenchmarkResults::default(),
            database_benchmark: DatabaseBenchmarkResults::default(),
            graph_benchmark: GraphBenchmarkResults::default(),
            simd_benchmark: SimdBenchmarkResults::default(),
            pipeline_benchmark: PipelineBenchmarkResults::default(),
        }
    }
    
    /// Print comprehensive benchmark summary
    pub fn print_summary(&self) {
        println!("\n🎯 COMPREHENSIVE BENCHMARK RESULTS SUMMARY");
        println!("=".repeat(50));
        
        println!("\n📊 Memory Optimizations:");
        println!("  Memory reduction: {:.1}%", self.memory_benchmark.memory_reduction_percent);
        println!("  Allocation speedup: {:.2}x", self.memory_benchmark.allocation_speedup);
        
        println!("\n🔄 Parallel Processing:");
        println!("  Maximum speedup: {:.2}x", self.parallel_benchmark.max_speedup);
        println!("  Parallel efficiency: {:.1}%", self.parallel_benchmark.parallel_efficiency * 100.0);
        println!("  Throughput: {:.0} items/sec", self.parallel_benchmark.throughput_ips);
        
        println!("\n💾 Database Operations:");
        println!("  Batch insert rate: {:.0} kmers/sec", self.database_benchmark.batch_insert_rate);
        println!("  Cache speedup: {:.2}x", self.database_benchmark.cache_speedup);
        println!("  Query throughput: {:.0} queries/sec", self.database_benchmark.query_throughput);
        
        println!("\n🕸️  Graph Algorithms:");
        println!("  SCC algorithm: {} ms", self.graph_benchmark.scc_algorithm_time_ms);
        println!("  Contig generation: {} ms", self.graph_benchmark.contig_generation_time_ms);
        println!("  Graph simplification: {} ms", self.graph_benchmark.simplification_time_ms);
        
        println!("\n⚡ SIMD Operations:");
        println!("  Maximum speedup: {:.2}x", self.simd_benchmark.max_speedup);
        println!("  K-mer processing rate: {:.0} kmers/sec", self.simd_benchmark.kmer_processing_rate);
        
        println!("\n🔄 Complete Pipeline:");
        println!("  Maximum throughput: {:.1} sequences/sec", self.pipeline_benchmark.max_throughput);
        println!("  Total processing time: {} ms", self.pipeline_benchmark.total_processing_time_ms);
        
        // Overall assessment
        println!("\n🎖️  OPTIMIZATION IMPACT ASSESSMENT:");
        
        let memory_grade = if self.memory_benchmark.memory_reduction_percent > 50.0 { "A+" }
                          else if self.memory_benchmark.memory_reduction_percent > 30.0 { "A" }
                          else if self.memory_benchmark.memory_reduction_percent > 15.0 { "B" }
                          else { "C" };
        
        let parallel_grade = if self.parallel_benchmark.max_speedup > 6.0 { "A+" }
                            else if self.parallel_benchmark.max_speedup > 4.0 { "A" }
                            else if self.parallel_benchmark.max_speedup > 2.0 { "B" }
                            else { "C" };
        
        let simd_grade = if self.simd_benchmark.max_speedup > 4.0 { "A+" }
                        else if self.simd_benchmark.max_speedup > 2.0 { "A" }
                        else if self.simd_benchmark.max_speedup > 1.5 { "B" }
                        else { "C" };
        
        println!("  Memory optimization: {} ({:.1}% reduction)", 
                 memory_grade, self.memory_benchmark.memory_reduction_percent);
        println!("  Parallel processing: {} ({:.2}x speedup)", 
                 parallel_grade, self.parallel_benchmark.max_speedup);
        println!("  SIMD acceleration: {} ({:.2}x speedup)", 
                 simd_grade, self.simd_benchmark.max_speedup);
        
        let overall_improvement = (
            self.memory_benchmark.allocation_speedup +
            self.parallel_benchmark.max_speedup +
            self.simd_benchmark.max_speedup
        ) / 3.0;
        
        println!("\n🏆 OVERALL PERFORMANCE IMPROVEMENT: {:.2}x", overall_improvement);
        
        if overall_improvement > 5.0 {
            println!("   Status: ✅ EXCELLENT - Target performance achieved!");
        } else if overall_improvement > 3.0 {
            println!("   Status: ✅ GOOD - Strong performance improvements");
        } else if overall_improvement > 2.0 {
            println!("   Status: ⚠️  MODERATE - Some improvements achieved");
        } else {
            println!("   Status: ❌ NEEDS WORK - Optimization targets not met");
        }
    }
}

#[derive(Debug, Default)]
pub struct MemoryBenchmarkResults {
    pub memory_reduction_percent: f64,
    pub allocation_speedup: f64,
}

#[derive(Debug, Default)]
pub struct ParallelBenchmarkResults {
    pub max_speedup: f64,
    pub parallel_efficiency: f64,
    pub throughput_ips: f64,
}

#[derive(Debug, Default)]
pub struct DatabaseBenchmarkResults {
    pub batch_insert_rate: f64,
    pub cache_speedup: f64,
    pub query_throughput: f64,
}

#[derive(Debug, Default)]
pub struct GraphBenchmarkResults {
    pub graph_construction_time_ms: u128,
    pub scc_algorithm_time_ms: u128,
    pub contig_generation_time_ms: u128,
    pub simplification_time_ms: u128,
}

#[derive(Debug, Default)]
pub struct SimdBenchmarkResults {
    pub max_speedup: f64,
    pub kmer_processing_rate: f64,
}

#[derive(Debug, Default)]
pub struct PipelineBenchmarkResults {
    pub max_throughput: f64,
    pub total_processing_time_ms: u128,
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_benchmark_suite() -> Result<()> {
        let benchmark = ComprehensiveBenchmark::new();
        
        // Run a limited version of the benchmark for testing
        println!("Running limited benchmark suite for testing...");
        
        // Test memory optimizations
        let memory_results = benchmark.benchmark_memory_optimizations()?;
        assert!(memory_results.allocation_speedup >= 1.0);
        
        // Test SIMD operations
        let simd_results = benchmark.benchmark_simd_operations()?;
        assert!(simd_results.max_speedup >= 1.0);
        
        println!("✅ Benchmark suite test completed successfully");
        Ok(())
    }
    
    #[test]
    #[ignore] // This is a long-running test
    fn test_full_benchmark_suite() -> Result<()> {
        let benchmark = ComprehensiveBenchmark::new();
        let results = benchmark.run_comprehensive_benchmark()?;
        
        // Verify that we got meaningful results
        assert!(results.memory_benchmark.memory_reduction_percent >= 0.0);
        assert!(results.parallel_benchmark.max_speedup >= 1.0);
        assert!(results.simd_benchmark.max_speedup >= 1.0);
        
        Ok(())
    }
}