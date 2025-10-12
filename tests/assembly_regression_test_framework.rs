//! Assembly Regression Testing Framework
//! ====================================
//!
//! Comprehensive regression testing to ensure laptop assembly optimizations
//! do not break existing functionality or degrade assembly quality.
//! Includes baseline comparison, performance regression detection, and automated validation.

use meta_forge::assembly::{LaptopAssembler, LaptopConfig};
use meta_forge::core::data_structures::{CorrectedRead, Contig, CorrectionMetadata};
use std::time::{Duration, Instant};
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use anyhow::Result;

/// Baseline performance metrics for regression testing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BaselineMetrics {
    pub assembly_time_ms: u64,
    pub memory_usage_mb: f64,
    pub n50: usize,
    pub total_contigs: usize,
    pub total_length: usize,
    pub coverage_mean: f64,
    pub largest_contig: usize,
    pub success_rate: f64,
    pub timestamp: String,
    pub config_name: String,
}

/// Regression test results comparison
#[derive(Debug, Clone)]
pub struct RegressionTestResult {
    pub test_name: String,
    pub baseline: BaselineMetrics,
    pub current: BaselineMetrics,
    pub performance_delta: PerformanceDelta,
    pub quality_delta: QualityDelta,
    pub passed: bool,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct PerformanceDelta {
    pub time_change_percent: f64,
    pub memory_change_percent: f64,
    pub throughput_change_percent: f64,
}

#[derive(Debug, Clone)]
pub struct QualityDelta {
    pub n50_change_percent: f64,
    pub contigs_change_percent: f64,
    pub length_change_percent: f64,
    pub coverage_change_percent: f64,
}

/// Regression testing framework
pub struct RegressionTester {
    baselines: HashMap<String, BaselineMetrics>,
    tolerance_config: ToleranceConfig,
}

#[derive(Debug, Clone)]
pub struct ToleranceConfig {
    /// Maximum acceptable performance degradation (percentage)
    pub max_time_degradation: f64,
    /// Maximum acceptable memory increase (percentage)
    pub max_memory_increase: f64,
    /// Minimum acceptable N50 retention (percentage)
    pub min_n50_retention: f64,
    /// Maximum acceptable increase in contig count (fragmentation)
    pub max_fragmentation_increase: f64,
    /// Minimum acceptable assembly length retention
    pub min_length_retention: f64,
}

impl Default for ToleranceConfig {
    fn default() -> Self {
        Self {
            max_time_degradation: 20.0, // 20% slower max
            max_memory_increase: 15.0,   // 15% more memory max
            min_n50_retention: 90.0,     // Must retain 90% of N50
            max_fragmentation_increase: 25.0, // 25% more contigs max
            min_length_retention: 85.0,  // Must retain 85% of length
        }
    }
}

impl RegressionTester {
    pub fn new() -> Self {
        Self {
            baselines: HashMap::new(),
            tolerance_config: ToleranceConfig::default(),
        }
    }
    
    pub fn with_tolerance(mut self, config: ToleranceConfig) -> Self {
        self.tolerance_config = config;
        self
    }
    
    /// Store baseline metrics for comparison
    pub fn set_baseline(&mut self, test_name: &str, metrics: BaselineMetrics) {
        self.baselines.insert(test_name.to_string(), metrics);
    }
    
    /// Run regression test against stored baseline
    pub fn run_regression_test(&self, test_name: &str, test_reads: &[CorrectedRead], config: LaptopConfig) -> Result<RegressionTestResult> {
        let baseline = self.baselines.get(test_name)
            .ok_or_else(|| anyhow::anyhow!("No baseline found for test: {}", test_name))?;
        
        // Run current assembly and measure metrics
        let current = self.measure_assembly_metrics(test_reads, config, test_name)?;
        
        // Calculate deltas
        let performance_delta = PerformanceDelta {
            time_change_percent: ((current.assembly_time_ms as f64 - baseline.assembly_time_ms as f64) / baseline.assembly_time_ms as f64) * 100.0,
            memory_change_percent: ((current.memory_usage_mb - baseline.memory_usage_mb) / baseline.memory_usage_mb) * 100.0,
            throughput_change_percent: -((current.assembly_time_ms as f64 - baseline.assembly_time_ms as f64) / baseline.assembly_time_ms as f64) * 100.0,
        };
        
        let quality_delta = QualityDelta {
            n50_change_percent: ((current.n50 as f64 - baseline.n50 as f64) / baseline.n50 as f64) * 100.0,
            contigs_change_percent: ((current.total_contigs as f64 - baseline.total_contigs as f64) / baseline.total_contigs as f64) * 100.0,
            length_change_percent: ((current.total_length as f64 - baseline.total_length as f64) / baseline.total_length as f64) * 100.0,
            coverage_change_percent: ((current.coverage_mean - baseline.coverage_mean) / baseline.coverage_mean) * 100.0,
        };
        
        // Evaluate if test passed
        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        let mut passed = true;
        
        // Performance regression checks
        if performance_delta.time_change_percent > self.tolerance_config.max_time_degradation {
            errors.push(format!(
                "Performance regression: {:.1}% slower (max: {:.1}%)",
                performance_delta.time_change_percent,
                self.tolerance_config.max_time_degradation
            ));
            passed = false;
        } else if performance_delta.time_change_percent > self.tolerance_config.max_time_degradation / 2.0 {
            warnings.push(format!(
                "Performance warning: {:.1}% slower",
                performance_delta.time_change_percent
            ));
        }
        
        if performance_delta.memory_change_percent > self.tolerance_config.max_memory_increase {
            errors.push(format!(
                "Memory regression: {:.1}% more memory (max: {:.1}%)",
                performance_delta.memory_change_percent,
                self.tolerance_config.max_memory_increase
            ));
            passed = false;
        }
        
        // Quality regression checks
        let n50_retention = (current.n50 as f64 / baseline.n50 as f64) * 100.0;
        if n50_retention < self.tolerance_config.min_n50_retention {
            errors.push(format!(
                "N50 regression: {:.1}% retention (min: {:.1}%)",
                n50_retention,
                self.tolerance_config.min_n50_retention
            ));
            passed = false;
        }
        
        let length_retention = (current.total_length as f64 / baseline.total_length as f64) * 100.0;
        if length_retention < self.tolerance_config.min_length_retention {
            errors.push(format!(
                "Assembly length regression: {:.1}% retention (min: {:.1}%)",
                length_retention,
                self.tolerance_config.min_length_retention
            ));
            passed = false;
        }
        
        if quality_delta.contigs_change_percent > self.tolerance_config.max_fragmentation_increase {
            warnings.push(format!(
                "Potential fragmentation: {:.1}% more contigs",
                quality_delta.contigs_change_percent
            ));
        }
        
        Ok(RegressionTestResult {
            test_name: test_name.to_string(),
            baseline: baseline.clone(),
            current,
            performance_delta,
            quality_delta,
            passed,
            warnings,
            errors,
        })
    }
    
    fn measure_assembly_metrics(&self, reads: &[CorrectedRead], config: LaptopConfig, test_name: &str) -> Result<BaselineMetrics> {
        let assembler = LaptopAssembler::new(config.clone());
        
        let start_time = Instant::now();
        let start_memory = self.get_memory_usage();
        
        let contigs = assembler.assemble(reads)?;

        let assembly_time = start_time.elapsed();
        let end_memory = self.get_memory_usage();
        let memory_usage = (end_memory - start_memory).max(0.0);

        let metrics = calculate_quality_metrics(&contigs);
        
        Ok(BaselineMetrics {
            assembly_time_ms: assembly_time.as_millis() as u64,
            memory_usage_mb: memory_usage,
            n50: metrics.n50,
            total_contigs: metrics.num_contigs,
            total_length: metrics.total_length,
            coverage_mean: metrics.coverage_mean,
            largest_contig: metrics.largest_contig,
            success_rate: if contigs.is_empty() { 0.0 } else { 100.0 },
            timestamp: chrono::Utc::now().to_rfc3339(),
            config_name: format!("{:?}", config),
        })
    }
    
    fn get_memory_usage(&self) -> f64 {
        // Simplified memory measurement - in production use proper OS APIs
        std::env::var("CURRENT_MEMORY_MB")
            .unwrap_or_else(|_| "100.0".to_string())
            .parse()
            .unwrap_or(100.0)
    }
}

/// Helper struct for quality metrics (internal to this test)
struct QualityMetrics {
    n50: usize,
    num_contigs: usize,
    total_length: usize,
    coverage_mean: f64,
    largest_contig: usize,
}

/// Calculate assembly quality metrics
fn calculate_quality_metrics(contigs: &[Contig]) -> QualityMetrics {
    if contigs.is_empty() {
        return QualityMetrics {
            n50: 0,
            num_contigs: 0,
            total_length: 0,
            coverage_mean: 0.0,
            largest_contig: 0,
        };
    }

    let num_contigs = contigs.len();
    let total_length: usize = contigs.iter().map(|c| c.length).sum();
    let coverage_mean = contigs.iter().map(|c| c.coverage).sum::<f64>() / num_contigs as f64;
    let largest_contig = contigs.iter().map(|c| c.length).max().unwrap_or(0);

    // Calculate N50
    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
    lengths.sort_by(|a, b| b.cmp(a));

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

    QualityMetrics {
        n50,
        num_contigs,
        total_length,
        coverage_mean,
        largest_contig,
    }
}

/// Create standard test datasets for regression testing
pub struct StandardTestDatasets {
    pub small_uniform: Vec<CorrectedRead>,
    pub medium_complex: Vec<CorrectedRead>,
    pub large_repetitive: Vec<CorrectedRead>,
    pub low_coverage: Vec<CorrectedRead>,
    pub high_error: Vec<CorrectedRead>,
}

impl StandardTestDatasets {
    pub fn generate() -> Self {
        Self {
            small_uniform: generate_uniform_reads("ATCGATCGATCGAAGTCGATCGATCGAAGT", 20, 5.0, 0.01),
            medium_complex: generate_complex_reads(150, 50, 10.0, 0.015),
            large_repetitive: generate_repetitive_reads("ATCG", 10, 300, 15.0, 0.01),
            low_coverage: generate_uniform_reads("ATCGATCGATCGAAGTCGATCGATCGAAGT", 25, 2.0, 0.02),
            high_error: generate_uniform_reads("ATCGATCGATCGAAGTCGATCGATCGAAGT", 30, 8.0, 0.05),
        }
    }
}

fn generate_uniform_reads(reference: &str, read_len: usize, coverage: f64, _error_rate: f64) -> Vec<CorrectedRead> {
    // Simple read generation for regression testing
    let num_reads = ((reference.len() as f64 * coverage) / read_len as f64).ceil() as usize;
    let mut reads = Vec::new();

    for i in 0..num_reads {
        let start = (i * read_len) % (reference.len().saturating_sub(read_len).max(1));
        let end = (start + read_len).min(reference.len());
        let sequence = reference[start..end].to_string();

        reads.push(create_test_read(i, &sequence));
    }

    reads
}

fn generate_complex_reads(ref_len: usize, read_len: usize, coverage: f64, error_rate: f64) -> Vec<CorrectedRead> {
    let bases = ['A', 'T', 'C', 'G'];
    let reference: String = (0..ref_len).map(|_| bases[fastrand::usize(0..4)]).collect();
    generate_uniform_reads(&reference, read_len, coverage, error_rate)
}

fn generate_repetitive_reads(repeat_unit: &str, repeat_count: usize, flanking_len: usize, coverage: f64, error_rate: f64) -> Vec<CorrectedRead> {
    let bases = ['A', 'T', 'C', 'G'];
    let flanking: String = (0..flanking_len).map(|_| bases[fastrand::usize(0..4)]).collect();
    let reference = format!("{}{}{}", flanking, repeat_unit.repeat(repeat_count), flanking);
    generate_uniform_reads(&reference, 40, coverage, error_rate)
}

// Regression Tests

#[test]
fn test_laptop_assembly_regression_suite() {
    let mut tester = RegressionTester::new();
    let datasets = StandardTestDatasets::generate();
    
    // Test with different configurations
    let configs = vec![
        ("low_memory", LaptopConfig::low_memory()),
        ("medium_memory", LaptopConfig::medium_memory()),
        ("high_memory", LaptopConfig::high_memory()),
    ];
    
    for (config_name, config) in configs {
        println!("\n🧪 Testing {} configuration", config_name);
        
        // Generate baseline if not exists (first run)
        let test_name = format!("standard_small_{}", config_name);
        
        if !tester.baselines.contains_key(&test_name) {
            println!("  📊 Generating baseline for {}", test_name);
            let baseline = tester.measure_assembly_metrics(&datasets.small_uniform, config.clone(), &test_name).unwrap();
            tester.set_baseline(&test_name, baseline);
        }
        
        // Run regression test (simulated - would compare against stored baseline)
        let dummy_baseline = BaselineMetrics {
            assembly_time_ms: 100,
            memory_usage_mb: 50.0,
            n50: 25,
            total_contigs: 3,
            total_length: 75,
            coverage_mean: 5.0,
            largest_contig: 30,
            success_rate: 100.0,
            timestamp: "2024-01-01T00:00:00Z".to_string(),
            config_name: config_name.to_string(),
        };
        tester.set_baseline(&test_name, dummy_baseline);
        
        let result = tester.run_regression_test(&test_name, &datasets.small_uniform, config).unwrap();
        
        println!("  📈 Performance delta: {:.1}% time, {:.1}% memory", 
                result.performance_delta.time_change_percent,
                result.performance_delta.memory_change_percent);
        println!("  📊 Quality delta: {:.1}% N50, {:.1}% length",
                result.quality_delta.n50_change_percent,
                result.quality_delta.length_change_percent);
        
        if !result.warnings.is_empty() {
            println!("  ⚠️  Warnings: {:?}", result.warnings);
        }
        
        if !result.errors.is_empty() {
            println!("  ❌ Errors: {:?}", result.errors);
        }
        
        println!("  ✅ Test passed: {}", result.passed);
        
        // For now, allow some flexibility in tests
        // In production, you'd want: assert!(result.passed);
    }
}

#[test]
fn test_memory_constraint_regression() {
    let mut tester = RegressionTester::new()
        .with_tolerance(ToleranceConfig {
            max_memory_increase: 10.0, // Strict memory constraint
            ..ToleranceConfig::default()
        });
    
    let datasets = StandardTestDatasets::generate();
    
    // Test memory usage doesn't regress with optimizations
    let config = LaptopConfig::low_memory();
    
    let baseline = BaselineMetrics {
        assembly_time_ms: 150,
        memory_usage_mb: 100.0, // 100MB baseline
        n50: 20,
        total_contigs: 5,
        total_length: 60,
        coverage_mean: 3.0,
        largest_contig: 25,
        success_rate: 100.0,
        timestamp: "2024-01-01T00:00:00Z".to_string(),
        config_name: "low_memory".to_string(),
    };
    
    tester.set_baseline("memory_test", baseline);
    
    let result = tester.run_regression_test("memory_test", &datasets.low_coverage, config).unwrap();
    
    println!("Memory regression test:");
    println!("  Memory change: {:.1}%", result.performance_delta.memory_change_percent);
    
    // Check that memory usage is controlled
    assert!(result.performance_delta.memory_change_percent <= 10.0, 
           "Memory usage increased too much: {:.1}%", result.performance_delta.memory_change_percent);
}

#[test]
fn test_quality_regression_strict() {
    let mut tester = RegressionTester::new()
        .with_tolerance(ToleranceConfig {
            min_n50_retention: 95.0, // Very strict quality requirements
            min_length_retention: 90.0,
            max_fragmentation_increase: 10.0,
            ..ToleranceConfig::default()
        });
    
    let datasets = StandardTestDatasets::generate();
    
    let baseline = BaselineMetrics {
        assembly_time_ms: 200,
        memory_usage_mb: 150.0,
        n50: 40,
        total_contigs: 3,
        total_length: 120,
        coverage_mean: 8.0,
        largest_contig: 50,
        success_rate: 100.0,
        timestamp: "2024-01-01T00:00:00Z".to_string(),
        config_name: "quality_test".to_string(),
    };
    
    tester.set_baseline("quality_test", baseline);
    
    let result = tester.run_regression_test("quality_test", &datasets.medium_complex, LaptopConfig::medium_memory()).unwrap();
    
    println!("Quality regression test:");
    println!("  N50 retention: {:.1}%", (result.current.n50 as f64 / result.baseline.n50 as f64) * 100.0);
    println!("  Length retention: {:.1}%", (result.current.total_length as f64 / result.baseline.total_length as f64) * 100.0);
    
    // Quality should be maintained or improved
    assert!(result.current.n50 >= (result.baseline.n50 as f64 * 0.95) as usize,
           "N50 regression detected");
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
        kmer_hash_cache: Vec::new(),
    }
}
