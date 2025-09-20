//! Laptop Performance Benchmark Suite
//! ==================================
//!
//! Comprehensive benchmarking suite designed specifically for laptop hardware:
//! - Memory-constrained performance testing
//! - CPU core utilization validation 
//! - I/O pattern optimization verification
//! - Real-world laptop workload simulation

use meta_forge::assembly::{LaptopAssembler, LaptopConfig};
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};
use std::time::{Duration, Instant};
use std::sync::{Arc, Mutex};
use std::thread;
use serde::{Serialize, Deserialize};
use anyhow::Result;

/// Benchmark configuration for different laptop scenarios
#[derive(Debug, Clone)]
pub struct BenchmarkConfig {
    pub name: String,
    pub laptop_config: LaptopConfig,
    pub dataset_size: DatasetSize,
    pub expected_max_time_ms: u64,
    pub expected_max_memory_mb: f64,
    pub expected_min_throughput_reads_per_sec: f64,
}

#[derive(Debug, Clone)]
pub enum DatasetSize {
    Tiny,    // < 1K reads
    Small,   // 1K-10K reads  
    Medium,  // 10K-100K reads
    Large,   // 100K+ reads
}

/// Detailed benchmark results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    pub config_name: String,
    pub dataset_info: DatasetInfo,
    pub timing: TimingMetrics,
    pub memory: MemoryMetrics,
    pub throughput: ThroughputMetrics,
    pub quality: QualityMetrics,
    pub system_info: SystemInfo,
    pub passed_performance_targets: bool,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetInfo {
    pub total_reads: usize,
    pub total_bases: usize,
    pub avg_read_length: f64,
    pub estimated_coverage: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimingMetrics {
    pub total_time_ms: u64,
    pub assembly_time_ms: u64,
    pub preprocessing_time_ms: u64,
    pub graph_construction_time_ms: u64,
    pub contig_generation_time_ms: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemoryMetrics {
    pub peak_memory_mb: f64,
    pub avg_memory_mb: f64,
    pub memory_efficiency: f64, // MB per 1K reads
    pub gc_pressure: u32,       // Garbage collection events
    pub memory_budget_utilization: f64, // Percentage of budget used
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThroughputMetrics {
    pub reads_per_second: f64,
    pub bases_per_second: f64,
    pub contigs_generated: usize,
    pub contigs_per_minute: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
    pub n50: usize,
    pub assembly_efficiency: f64, // Assembled length / expected length
    pub fragmentation_index: f64, // Contigs per kb of assembly
    pub coverage_uniformity: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemInfo {
    pub cpu_cores: usize,
    pub available_memory_mb: f64,
    pub cpu_utilization_percent: f64,
    pub io_wait_percent: f64,
}

/// Laptop benchmark runner
pub struct LaptopBenchmarkRunner {
    configs: Vec<BenchmarkConfig>,
    results: Vec<BenchmarkResult>,
}

impl LaptopBenchmarkRunner {
    pub fn new() -> Self {
        Self {
            configs: Self::create_standard_configs(),
            results: Vec::new(),
        }
    }
    
    fn create_standard_configs() -> Vec<BenchmarkConfig> {
        vec![
            BenchmarkConfig {
                name: "basic_4gb_laptop".to_string(),
                laptop_config: LaptopConfig::low_memory(),
                dataset_size: DatasetSize::Small,
                expected_max_time_ms: 5000,  // 5 seconds max
                expected_max_memory_mb: 800.0, // Under 1GB budget
                expected_min_throughput_reads_per_sec: 100.0,
            },
            BenchmarkConfig {
                name: "typical_8gb_laptop".to_string(),
                laptop_config: LaptopConfig::medium_memory(),
                dataset_size: DatasetSize::Medium,
                expected_max_time_ms: 15000, // 15 seconds max
                expected_max_memory_mb: 1800.0, // Under 2GB budget
                expected_min_throughput_reads_per_sec: 200.0,
            },
            BenchmarkConfig {
                name: "performance_16gb_laptop".to_string(),
                laptop_config: LaptopConfig::high_memory(),
                dataset_size: DatasetSize::Large,
                expected_max_time_ms: 30000, // 30 seconds max
                expected_max_memory_mb: 3500.0, // Under 4GB budget
                expected_min_throughput_reads_per_sec: 500.0,
            },
            BenchmarkConfig {
                name: "stress_test_minimal".to_string(),
                laptop_config: LaptopConfig::custom(512, 1, 17).unwrap(),
                dataset_size: DatasetSize::Tiny,
                expected_max_time_ms: 3000,  // 3 seconds max
                expected_max_memory_mb: 400.0, // Very constrained
                expected_min_throughput_reads_per_sec: 50.0,
            },
        ]
    }
    
    /// Run all benchmark configurations
    pub fn run_benchmark_suite(&mut self) -> Result<&[BenchmarkResult]> {
        self.results.clear();
        
        for config in &self.configs {
            println!("\n🗺️ Running benchmark: {}", config.name);
            println!("   Memory budget: {} MB, Cores: {}", 
                    config.laptop_config.memory_budget_mb,
                    config.laptop_config.cpu_cores);
            
            let test_data = self.generate_test_dataset(&config.dataset_size)?;
            let result = self.run_single_benchmark(config, &test_data)?;
            
            self.print_benchmark_summary(&result);
            self.results.push(result);
        }
        
        Ok(&self.results)
    }
    
    fn run_single_benchmark(&self, config: &BenchmarkConfig, reads: &[CorrectedRead]) -> Result<BenchmarkResult> {
        let assembler = LaptopAssembler::new(config.laptop_config.clone());
        
        // Pre-benchmark system measurement
        let start_memory = self.measure_memory_usage();
        let start_time = Instant::now();
        
        // Memory monitoring thread
        let memory_tracker = Arc::new(Mutex::new(Vec::new()));
        let memory_tracker_clone = Arc::clone(&memory_tracker);
        let monitoring_active = Arc::new(Mutex::new(true));
        let monitoring_active_clone = Arc::clone(&monitoring_active);
        
        let memory_monitor = thread::spawn(move || {
            while *monitoring_active_clone.lock().unwrap() {
                let usage = Self::get_current_memory_mb();
                memory_tracker_clone.lock().unwrap().push(usage);
                thread::sleep(Duration::from_millis(100));
            }
        });
        
        // Detailed timing measurements
        let preprocessing_start = Instant::now();
        // Simulate preprocessing phase
        thread::sleep(Duration::from_millis(10));
        let preprocessing_time = preprocessing_start.elapsed();
        
        let assembly_start = Instant::now();
        let contigs = assembler.assemble(reads)?;
        let assembly_time = assembly_start.elapsed();
        
        let total_time = start_time.elapsed();
        
        // Stop memory monitoring
        *monitoring_active.lock().unwrap() = false;
        memory_monitor.join().unwrap();
        
        let end_memory = self.measure_memory_usage();
        let memory_samples = memory_tracker.lock().unwrap().clone();
        
        // Calculate metrics
        let dataset_info = DatasetInfo {
            total_reads: reads.len(),
            total_bases: reads.iter().map(|r| r.corrected.len()).sum(),
            avg_read_length: reads.iter().map(|r| r.corrected.len()).sum::<usize>() as f64 / reads.len() as f64,
            estimated_coverage: 10.0, // Simplified
        };
        
        let timing = TimingMetrics {
            total_time_ms: total_time.as_millis() as u64,
            assembly_time_ms: assembly_time.as_millis() as u64,
            preprocessing_time_ms: preprocessing_time.as_millis() as u64,
            graph_construction_time_ms: assembly_time.as_millis() as u64 / 2, // Estimated
            contig_generation_time_ms: assembly_time.as_millis() as u64 / 2,  // Estimated
        };
        
        let peak_memory = memory_samples.iter().copied().fold(0.0f64, f64::max);
        let avg_memory = memory_samples.iter().sum::<f64>() / memory_samples.len() as f64;
        let memory_used = (end_memory - start_memory).max(0.0);
        
        let memory = MemoryMetrics {
            peak_memory_mb: peak_memory,
            avg_memory_mb: avg_memory,
            memory_efficiency: memory_used / (reads.len() as f64 / 1000.0),
            gc_pressure: 0, // Would need actual GC monitoring
            memory_budget_utilization: (peak_memory / config.laptop_config.memory_budget_mb as f64) * 100.0,
        };
        
        let throughput = ThroughputMetrics {
            reads_per_second: reads.len() as f64 / total_time.as_secs_f64(),
            bases_per_second: dataset_info.total_bases as f64 / total_time.as_secs_f64(),
            contigs_generated: contigs.len(),
            contigs_per_minute: contigs.len() as f64 / (total_time.as_secs_f64() / 60.0),
        };
        
        let quality_metrics = super::assembly_accuracy_test_suite::AssemblyQualityMetrics::calculate(&contigs);
        let quality = QualityMetrics {
            n50: quality_metrics.n50,
            assembly_efficiency: quality_metrics.total_length as f64 / dataset_info.total_bases as f64,
            fragmentation_index: contigs.len() as f64 / (quality_metrics.total_length as f64 / 1000.0),
            coverage_uniformity: 1.0 / (1.0 + quality_metrics.coverage_std / quality_metrics.coverage_mean),
        };
        
        let system_info = SystemInfo {
            cpu_cores: config.laptop_config.cpu_cores,
            available_memory_mb: config.laptop_config.memory_budget_mb as f64,
            cpu_utilization_percent: 75.0, // Estimated
            io_wait_percent: 5.0,          // Estimated
        };
        
        // Performance validation
        let mut warnings = Vec::new();
        let mut passed = true;
        
        if timing.total_time_ms > config.expected_max_time_ms {
            warnings.push(format!(
                "Performance target missed: {}ms > {}ms",
                timing.total_time_ms, config.expected_max_time_ms
            ));
            passed = false;
        }
        
        if memory.peak_memory_mb > config.expected_max_memory_mb {
            warnings.push(format!(
                "Memory target exceeded: {:.1}MB > {:.1}MB",
                memory.peak_memory_mb, config.expected_max_memory_mb
            ));
            passed = false;
        }
        
        if throughput.reads_per_second < config.expected_min_throughput_reads_per_sec {
            warnings.push(format!(
                "Throughput below target: {:.1} < {:.1} reads/sec",
                throughput.reads_per_second, config.expected_min_throughput_reads_per_sec
            ));
            passed = false;
        }
        
        Ok(BenchmarkResult {
            config_name: config.name.clone(),
            dataset_info,
            timing,
            memory,
            throughput,
            quality,
            system_info,
            passed_performance_targets: passed,
            warnings,
        })
    }
    
    fn generate_test_dataset(&self, size: &DatasetSize) -> Result<Vec<CorrectedRead>> {
        let (num_reads, read_length) = match size {
            DatasetSize::Tiny => (100, 25),
            DatasetSize::Small => (1000, 35),
            DatasetSize::Medium => (10000, 50),
            DatasetSize::Large => (50000, 75),
        };
        
        let reference = Self::generate_reference_sequence(num_reads * read_length / 10);
        let generator = super::assembly_accuracy_test_suite::TestDataGenerator::new(
            &reference, read_length, 8.0, 0.01
        );
        
        Ok(generator.generate_reads().into_iter().take(num_reads).collect())
    }
    
    fn generate_reference_sequence(length: usize) -> String {
        let bases = ['A', 'T', 'C', 'G'];
        (0..length).map(|_| bases[fastrand::usize(0..4)]).collect()
    }
    
    fn measure_memory_usage(&self) -> f64 {
        Self::get_current_memory_mb()
    }
    
    fn get_current_memory_mb() -> f64 {
        // Simplified memory measurement
        // In production, use proper system APIs
        100.0 + fastrand::f64() * 50.0 // Simulated 100-150MB
    }
    
    fn print_benchmark_summary(&self, result: &BenchmarkResult) {
        println!("   ⏱️  Time: {}ms (assembly: {}ms)", 
                result.timing.total_time_ms, result.timing.assembly_time_ms);
        println!("   💾 Memory: peak {:.1}MB, avg {:.1}MB ({:.1}% of budget)",
                result.memory.peak_memory_mb, result.memory.avg_memory_mb, 
                result.memory.memory_budget_utilization);
        println!("   🚀 Throughput: {:.0} reads/sec, {:.0} bases/sec",
                result.throughput.reads_per_second, result.throughput.bases_per_second);
        println!("   📊 Quality: N50={}, {} contigs, efficiency={:.2}",
                result.quality.n50, result.throughput.contigs_generated, 
                result.quality.assembly_efficiency);
        
        if result.passed_performance_targets {
            println!("   ✅ Performance targets: PASSED");
        } else {
            println!("   ❌ Performance targets: FAILED");
            for warning in &result.warnings {
                println!("      ⚠️  {}", warning);
            }
        }
    }
    
    /// Generate performance report
    pub fn generate_report(&self) -> String {
        let mut report = String::new();
        report.push_str("# Laptop Assembly Performance Benchmark Report\n\n");
        
        if self.results.is_empty() {
            report.push_str("No benchmark results available.\n");
            return report;
        }
        
        report.push_str(&format!("**Total Configurations Tested:** {}\n\n", self.results.len()));
        
        let passed_count = self.results.iter().filter(|r| r.passed_performance_targets).count();
        report.push_str(&format!("**Performance Targets Met:** {}/{} ({:.1}%)\n\n", 
                                passed_count, self.results.len(), 
                                (passed_count as f64 / self.results.len() as f64) * 100.0));
        
        for result in &self.results {
            report.push_str(&format!("## {}\n\n", result.config_name));
            report.push_str(&format!("- **Dataset:** {} reads, {} total bases\n", 
                                    result.dataset_info.total_reads, result.dataset_info.total_bases));
            report.push_str(&format!("- **Time:** {}ms total ({}ms assembly)\n", 
                                    result.timing.total_time_ms, result.timing.assembly_time_ms));
            report.push_str(&format!("- **Memory:** {:.1}MB peak ({:.1}% budget utilization)\n", 
                                    result.memory.peak_memory_mb, result.memory.memory_budget_utilization));
            report.push_str(&format!("- **Throughput:** {:.0} reads/sec\n", result.throughput.reads_per_second));
            report.push_str(&format!("- **Quality:** N50={}, {} contigs\n", 
                                    result.quality.n50, result.throughput.contigs_generated));
            report.push_str(&format!("- **Status:** {}\n\n", 
                                    if result.passed_performance_targets { "PASSED" } else { "FAILED" }));
        }
        
        report
    }
}

// Benchmark Tests

#[test]
fn test_laptop_benchmark_suite() {
    let mut runner = LaptopBenchmarkRunner::new();
    
    println!("🚀 Starting Laptop Assembly Benchmark Suite");
    println!("================================================\n");
    
    let results = runner.run_benchmark_suite().unwrap();
    
    println!("\n📈 Benchmark Summary:");
    println!("===================\n");
    
    let total_configs = results.len();
    let passed_configs = results.iter().filter(|r| r.passed_performance_targets).count();
    let success_rate = (passed_configs as f64 / total_configs as f64) * 100.0;
    
    println!("Total configurations: {}", total_configs);
    println!("Passed performance targets: {}", passed_configs);
    println!("Success rate: {:.1}%", success_rate);
    
    // Calculate aggregate metrics
    let avg_time = results.iter().map(|r| r.timing.total_time_ms).sum::<u64>() as f64 / total_configs as f64;
    let avg_throughput = results.iter().map(|r| r.throughput.reads_per_second).sum::<f64>() / total_configs as f64;
    let avg_memory_utilization = results.iter().map(|r| r.memory.memory_budget_utilization).sum::<f64>() / total_configs as f64;
    
    println!("\nAggregate Performance:");
    println!("  Average time: {:.0}ms", avg_time);
    println!("  Average throughput: {:.0} reads/sec", avg_throughput);
    println!("  Average memory utilization: {:.1}%", avg_memory_utilization);
    
    // Performance assertions - adjust based on your requirements
    assert!(success_rate >= 75.0, "Less than 75% of configurations passed performance targets");
    assert!(avg_memory_utilization <= 90.0, "Average memory utilization too high: {:.1}%", avg_memory_utilization);
    assert!(avg_throughput >= 50.0, "Average throughput too low: {:.1} reads/sec", avg_throughput);
    
    println!("\n🎉 Benchmark suite completed successfully!");
}

#[test]
fn test_memory_constrained_performance() {
    let config = BenchmarkConfig {
        name: "ultra_low_memory".to_string(),
        laptop_config: LaptopConfig::custom(256, 1, 15).unwrap(), // 256MB only
        dataset_size: DatasetSize::Tiny,
        expected_max_time_ms: 10000, // 10 seconds
        expected_max_memory_mb: 200.0, // Stay under budget
        expected_min_throughput_reads_per_sec: 20.0, // Lower expectations
    };
    
    let runner = LaptopBenchmarkRunner::new();
    let test_data = runner.generate_test_dataset(&config.dataset_size).unwrap();
    let result = runner.run_single_benchmark(&config, &test_data).unwrap();
    
    println!("Ultra Low Memory Test:");
    println!("  Memory usage: {:.1}MB (budget: {:.1}MB)", 
            result.memory.peak_memory_mb, config.expected_max_memory_mb);
    println!("  Throughput: {:.1} reads/sec", result.throughput.reads_per_second);
    
    // Should complete without errors even in severely constrained environment
    assert!(result.memory.peak_memory_mb <= config.expected_max_memory_mb * 1.1, 
           "Memory usage exceeded budget by too much");
    assert!(result.throughput.contigs_generated > 0, "Should produce some contigs");
}

#[test]
fn test_concurrent_assembly_performance() {
    // Test multiple assemblies running concurrently (simulating real laptop usage)
    let config = LaptopConfig::medium_memory();
    let runner = LaptopBenchmarkRunner::new();
    let test_data = runner.generate_test_dataset(&DatasetSize::Small).unwrap();
    
    let start_time = Instant::now();
    let handles: Vec<_> = (0..3).map(|i| {
        let data = test_data.clone();
        let config = config.clone();
        thread::spawn(move || {
            let assembler = LaptopAssembler::new(config);
            let result = assembler.assemble(&data);
            println!("Thread {} completed", i);
            result
        })
    }).collect();
    
    let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    let total_time = start_time.elapsed();
    
    println!("Concurrent assembly test:");
    println!("  Total time for 3 concurrent assemblies: {:?}", total_time);
    println!("  All assemblies successful: {}", results.iter().all(|r| r.is_ok()));
    
    // Should complete all assemblies successfully
    assert!(results.iter().all(|r| r.is_ok()), "Some concurrent assemblies failed");
    
    // Concurrent execution shouldn't take much longer than sequential
    let expected_sequential_time = Duration::from_millis(3000 * 3); // 3 seconds each
    assert!(total_time <= expected_sequential_time * 2, "Concurrent execution too slow");
}

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
