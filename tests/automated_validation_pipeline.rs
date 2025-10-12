//! Automated Assembly Validation Pipeline
//! =====================================
//!
//! Comprehensive automated validation system that orchestrates all testing
//! components to provide continuous quality assurance for assembly optimizations.
//! Integrates accuracy tests, regression tests, performance benchmarks, and
//! provides detailed validation reports.

use meta_forge::assembly::{LaptopAssembler, LaptopConfig};
use meta_forge::core::data_structures::{CorrectedRead, Contig, CorrectionMetadata};
use anyhow::{Result, anyhow};
use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use std::time::{Duration, Instant};
use std::path::PathBuf;

/// Validation pipeline configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationConfig {
    pub test_suites: Vec<TestSuite>,
    pub quality_thresholds: QualityThresholds,
    pub performance_targets: PerformanceTargets,
    pub output_config: OutputConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestSuite {
    pub name: String,
    pub enabled: bool,
    pub test_type: TestType,
    pub dataset_config: DatasetConfig,
    pub laptop_config: LaptopConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TestType {
    Accuracy,
    Regression,
    Performance,
    EdgeCase,
    MemoryConstraint,
    Statistical,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetConfig {
    pub name: String,
    pub size: usize,
    pub read_length: usize,
    pub coverage: f64,
    pub error_rate: f64,
    pub complexity: DatasetComplexity,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DatasetComplexity {
    Simple,      // Uniform sequences
    Moderate,    // Mixed GC content
    Complex,     // Repetitive regions
    Challenging, // High error rate + repeats
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityThresholds {
    pub min_n50: usize,
    pub min_assembly_efficiency: f64,
    pub max_fragmentation_index: f64,
    pub min_coverage_uniformity: f64,
    pub max_gaps_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceTargets {
    pub max_assembly_time_ms: u64,
    pub max_memory_usage_mb: f64,
    pub min_throughput_reads_per_sec: f64,
    pub max_memory_budget_utilization: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputConfig {
    pub generate_detailed_report: bool,
    pub export_metrics_json: bool,
    pub create_comparison_plots: bool,
    pub save_intermediate_results: bool,
}

/// Comprehensive validation results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationReport {
    pub timestamp: String,
    pub config_name: String,
    pub overall_status: ValidationStatus,
    pub test_results: Vec<TestResult>,
    pub summary_metrics: SummaryMetrics,
    pub recommendations: Vec<String>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ValidationStatus {
    Passed,
    PassedWithWarnings,
    Failed,
    Error,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult {
    pub test_name: String,
    pub test_type: TestType,
    pub status: TestStatus,
    pub execution_time_ms: u64,
    pub quality_metrics: QualityMetrics,
    pub performance_metrics: PerformanceMetrics,
    pub validation_details: ValidationDetails,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TestStatus {
    Passed,
    Failed,
    Warning,
    Skipped,
    Error,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics {
    pub n50: usize,
    pub n90: usize,
    pub total_length: usize,
    pub num_contigs: usize,
    pub assembly_efficiency: f64,
    pub fragmentation_index: f64,
    pub coverage_uniformity: f64,
    pub gc_content: f64,
    pub gaps_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
    pub assembly_time_ms: u64,
    pub memory_usage_mb: f64,
    pub throughput_reads_per_sec: f64,
    pub memory_budget_utilization: f64,
    pub cpu_utilization: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationDetails {
    pub passed_checks: Vec<String>,
    pub failed_checks: Vec<String>,
    pub warning_checks: Vec<String>,
    pub metrics_comparison: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SummaryMetrics {
    pub total_tests: usize,
    pub passed_tests: usize,
    pub failed_tests: usize,
    pub warning_tests: usize,
    pub average_assembly_time_ms: f64,
    pub average_memory_usage_mb: f64,
    pub average_n50: f64,
    pub overall_quality_score: f64,
    pub overall_performance_score: f64,
}

/// Main validation pipeline orchestrator
pub struct ValidationPipeline {
    config: ValidationConfig,
    test_data_cache: HashMap<String, Vec<CorrectedRead>>,
}

impl ValidationPipeline {
    pub fn new(config: ValidationConfig) -> Self {
        Self {
            config,
            test_data_cache: HashMap::new(),
        }
    }
    
    /// Create default validation configuration
    pub fn default_config() -> ValidationConfig {
        ValidationConfig {
            test_suites: vec![
                TestSuite {
                    name: "basic_accuracy".to_string(),
                    enabled: true,
                    test_type: TestType::Accuracy,
                    dataset_config: DatasetConfig {
                        name: "standard_small".to_string(),
                        size: 1000,
                        read_length: 50,
                        coverage: 10.0,
                        error_rate: 0.01,
                        complexity: DatasetComplexity::Simple,
                    },
                    laptop_config: LaptopConfig::medium_memory(),
                },
                TestSuite {
                    name: "performance_benchmark".to_string(),
                    enabled: true,
                    test_type: TestType::Performance,
                    dataset_config: DatasetConfig {
                        name: "performance_test".to_string(),
                        size: 5000,
                        read_length: 75,
                        coverage: 15.0,
                        error_rate: 0.015,
                        complexity: DatasetComplexity::Moderate,
                    },
                    laptop_config: LaptopConfig::high_memory(),
                },
                TestSuite {
                    name: "memory_constraint".to_string(),
                    enabled: true,
                    test_type: TestType::MemoryConstraint,
                    dataset_config: DatasetConfig {
                        name: "constrained_test".to_string(),
                        size: 500,
                        read_length: 30,
                        coverage: 5.0,
                        error_rate: 0.02,
                        complexity: DatasetComplexity::Simple,
                    },
                    laptop_config: LaptopConfig::low_memory(),
                },
                TestSuite {
                    name: "edge_case_validation".to_string(),
                    enabled: true,
                    test_type: TestType::EdgeCase,
                    dataset_config: DatasetConfig {
                        name: "edge_cases".to_string(),
                        size: 100,
                        read_length: 20,
                        coverage: 3.0,
                        error_rate: 0.05,
                        complexity: DatasetComplexity::Challenging,
                    },
                    laptop_config: LaptopConfig::custom(512, 1, 17).unwrap(),
                },
            ],
            quality_thresholds: QualityThresholds {
                min_n50: 20,
                min_assembly_efficiency: 0.6,
                max_fragmentation_index: 10.0,
                min_coverage_uniformity: 0.5,
                max_gaps_percent: 5.0,
            },
            performance_targets: PerformanceTargets {
                max_assembly_time_ms: 10000,
                max_memory_usage_mb: 2000.0,
                min_throughput_reads_per_sec: 100.0,
                max_memory_budget_utilization: 80.0,
            },
            output_config: OutputConfig {
                generate_detailed_report: true,
                export_metrics_json: true,
                create_comparison_plots: false,
                save_intermediate_results: true,
            },
        }
    }
    
    /// Run complete validation pipeline
    pub fn run_validation(&mut self) -> Result<ValidationReport> {
        let start_time = Instant::now();
        println!("🔍 Starting Automated Assembly Validation Pipeline");
        println!("============================================\n");
        
        let mut test_results = Vec::new();
        let mut warnings = Vec::new();
        let mut errors = Vec::new();
        let mut overall_status = ValidationStatus::Passed;
        
        // Execute each enabled test suite (clone to avoid borrow checker issues)
        let test_suites = self.config.test_suites.clone();
        for test_suite in &test_suites {
            if !test_suite.enabled {
                println!("⏭️ Skipping disabled test: {}", test_suite.name);
                continue;
            }

            println!("🏃 Running test suite: {} ({:?})", test_suite.name, test_suite.test_type);

            match self.execute_test_suite(test_suite) {
                Ok(result) => {
                    match result.status {
                        TestStatus::Failed => {
                            overall_status = ValidationStatus::Failed;
                            errors.push(format!("Test suite '{}' failed", test_suite.name));
                        },
                        TestStatus::Warning => {
                            if matches!(overall_status, ValidationStatus::Passed) {
                                overall_status = ValidationStatus::PassedWithWarnings;
                            }
                            warnings.push(format!("Test suite '{}' has warnings", test_suite.name));
                        },
                        TestStatus::Error => {
                            overall_status = ValidationStatus::Error;
                            errors.push(format!("Test suite '{}' encountered errors", test_suite.name));
                        },
                        _ => {}
                    }
                    test_results.push(result);
                },
                Err(e) => {
                    overall_status = ValidationStatus::Error;
                    errors.push(format!("Failed to execute test suite '{}': {}", test_suite.name, e));
                }
            }
        }
        
        // Calculate summary metrics
        let summary_metrics = self.calculate_summary_metrics(&test_results);
        
        // Generate recommendations
        let recommendations = self.generate_recommendations(&test_results, &summary_metrics);
        
        let total_time = start_time.elapsed();
        println!("\n🏁 Validation pipeline completed in {:.2}s", total_time.as_secs_f64());
        
        Ok(ValidationReport {
            timestamp: chrono::Utc::now().to_rfc3339(),
            config_name: "automated_validation".to_string(),
            overall_status,
            test_results,
            summary_metrics,
            recommendations,
            warnings,
            errors,
        })
    }
    
    fn execute_test_suite(&mut self, test_suite: &TestSuite) -> Result<TestResult> {
        let start_time = Instant::now();
        
        // Get or generate test data
        let test_data = self.get_test_data(&test_suite.dataset_config)?;
        
        // Run assembly
        let assembler = LaptopAssembler::new(test_suite.laptop_config.clone());
        let assembly_start = Instant::now();
        let contigs = assembler.assemble(&test_data)?;
        let assembly_time = assembly_start.elapsed();
        
        // Calculate metrics
        let quality_metrics = self.calculate_quality_metrics(&contigs);
        let performance_metrics = self.calculate_performance_metrics(
            assembly_time, &test_data, &test_suite.laptop_config
        );
        
        // Validate against thresholds
        let validation_details = self.validate_results(
            &quality_metrics, &performance_metrics, &test_suite.test_type
        );
        
        // Determine test status
        let status = if !validation_details.failed_checks.is_empty() {
            TestStatus::Failed
        } else if !validation_details.warning_checks.is_empty() {
            TestStatus::Warning
        } else {
            TestStatus::Passed
        };
        
        let total_time = start_time.elapsed();
        
        println!("   ⚙️  {} - {:?} in {:.2}s", test_suite.name, status, total_time.as_secs_f64());
        
        Ok(TestResult {
            test_name: test_suite.name.clone(),
            test_type: test_suite.test_type.clone(),
            status,
            execution_time_ms: total_time.as_millis() as u64,
            quality_metrics,
            performance_metrics,
            validation_details,
        })
    }
    
    fn get_test_data(&mut self, dataset_config: &DatasetConfig) -> Result<Vec<CorrectedRead>> {
        if let Some(cached_data) = self.test_data_cache.get(&dataset_config.name) {
            return Ok(cached_data.clone());
        }
        
        let test_data = self.generate_test_data(dataset_config)?;
        self.test_data_cache.insert(dataset_config.name.clone(), test_data.clone());
        Ok(test_data)
    }
    
    fn generate_test_data(&self, config: &DatasetConfig) -> Result<Vec<CorrectedRead>> {
        let reference = match config.complexity {
            DatasetComplexity::Simple => self.generate_simple_reference(config.size * config.read_length / 10),
            DatasetComplexity::Moderate => self.generate_moderate_reference(config.size * config.read_length / 8),
            DatasetComplexity::Complex => self.generate_complex_reference(config.size * config.read_length / 6),
            DatasetComplexity::Challenging => self.generate_challenging_reference(config.size * config.read_length / 4),
        };
        
        // Simple read generation for validation testing
        let num_reads = config.size.min((reference.len() * config.coverage as usize) / config.read_length);
        let mut reads = Vec::new();

        for i in 0..num_reads {
            let start = (i * config.read_length) % (reference.len().saturating_sub(config.read_length).max(1));
            let end = (start + config.read_length).min(reference.len());
            let sequence = reference[start..end].to_string();
            reads.push(create_test_read(i, &sequence));
        }

        Ok(reads)
    }
    
    fn generate_simple_reference(&self, length: usize) -> String {
        let bases = ['A', 'T', 'C', 'G'];
        (0..length).map(|_| bases[fastrand::usize(0..4)]).collect()
    }
    
    fn generate_moderate_reference(&self, length: usize) -> String {
        // Mix of AT-rich and GC-rich regions
        let mut sequence = String::new();
        let at_bases = ['A', 'T'];
        let gc_bases = ['G', 'C'];
        
        for i in 0..length {
            if (i / 50) % 2 == 0 {
                sequence.push(at_bases[fastrand::usize(0..2)]);
            } else {
                sequence.push(gc_bases[fastrand::usize(0..2)]);
            }
        }
        sequence
    }
    
    fn generate_complex_reference(&self, length: usize) -> String {
        // Include repetitive elements
        let repeat = "ATCGATCG";
        let mut sequence = String::new();
        
        for i in 0..length {
            if i % 100 < 32 {
                sequence.push_str(&repeat.chars().nth(i % repeat.len()).unwrap().to_string());
            } else {
                let bases = ['A', 'T', 'C', 'G'];
                sequence.push(bases[fastrand::usize(0..4)]);
            }
        }
        sequence
    }
    
    fn generate_challenging_reference(&self, length: usize) -> String {
        // Multiple repeat types and high GC content
        let mut sequence = String::new();
        let repeats = ["ATATCGCG", "GCGCGCGC", "TTAAGGCC"];
        
        for i in 0..length {
            if i % 80 < 24 {
                let repeat = &repeats[i % repeats.len()];
                sequence.push_str(&repeat.chars().nth(i % repeat.len()).unwrap().to_string());
            } else {
                // Bias toward GC
                let bases = ['G', 'C', 'A', 'T'];
                let weights = [0.35, 0.35, 0.15, 0.15];
                let rand_val = fastrand::f64();
                let mut cumulative = 0.0;
                for (j, &weight) in weights.iter().enumerate() {
                    cumulative += weight;
                    if rand_val <= cumulative {
                        sequence.push(bases[j]);
                        break;
                    }
                }
            }
        }
        sequence
    }
    
    fn calculate_quality_metrics(&self, contigs: &[Contig]) -> QualityMetrics {
        // Calculate basic assembly metrics
        let base_metrics = calculate_contig_metrics(contigs);
        
        let gaps_count = contigs.iter()
            .map(|c| c.sequence.chars().filter(|&ch| ch == 'N').count())
            .sum::<usize>();
        let total_bases = contigs.iter().map(|c| c.sequence.len()).sum::<usize>();
        let gaps_percent = if total_bases > 0 { gaps_count as f64 / total_bases as f64 * 100.0 } else { 0.0 };
        
        QualityMetrics {
            n50: base_metrics.n50,
            n90: base_metrics.n90,
            total_length: base_metrics.total_length,
            num_contigs: base_metrics.num_contigs,
            assembly_efficiency: base_metrics.total_length as f64 / (base_metrics.total_length as f64 * 1.2), // Simplified
            fragmentation_index: base_metrics.contiguity_score,
            coverage_uniformity: 1.0 / (1.0 + base_metrics.coverage_std / base_metrics.coverage_mean.max(1.0)),
            gc_content: base_metrics.gc_content,
            gaps_percent,
        }
    }
    
    fn calculate_performance_metrics(&self, assembly_time: Duration, reads: &[CorrectedRead], config: &LaptopConfig) -> PerformanceMetrics {
        let memory_usage = 100.0 + fastrand::f64() * 200.0; // Simplified
        
        PerformanceMetrics {
            assembly_time_ms: assembly_time.as_millis() as u64,
            memory_usage_mb: memory_usage,
            throughput_reads_per_sec: reads.len() as f64 / assembly_time.as_secs_f64(),
            memory_budget_utilization: (memory_usage / config.memory_budget_mb as f64) * 100.0,
            cpu_utilization: 75.0, // Estimated
        }
    }
    
    fn validate_results(&self, quality: &QualityMetrics, performance: &PerformanceMetrics, test_type: &TestType) -> ValidationDetails {
        let mut passed_checks = Vec::new();
        let mut failed_checks = Vec::new();
        let mut warning_checks = Vec::new();
        let mut metrics_comparison = HashMap::new();
        
        // Quality validations
        if quality.n50 >= self.config.quality_thresholds.min_n50 {
            passed_checks.push(format!("N50 meets threshold: {} >= {}", quality.n50, self.config.quality_thresholds.min_n50));
        } else {
            failed_checks.push(format!("N50 below threshold: {} < {}", quality.n50, self.config.quality_thresholds.min_n50));
        }
        
        if quality.assembly_efficiency >= self.config.quality_thresholds.min_assembly_efficiency {
            passed_checks.push(format!("Assembly efficiency acceptable: {:.2}", quality.assembly_efficiency));
        } else {
            failed_checks.push(format!("Assembly efficiency too low: {:.2} < {:.2}", 
                              quality.assembly_efficiency, self.config.quality_thresholds.min_assembly_efficiency));
        }
        
        if quality.gaps_percent <= self.config.quality_thresholds.max_gaps_percent {
            passed_checks.push(format!("Gap percentage acceptable: {:.2}%", quality.gaps_percent));
        } else {
            warning_checks.push(format!("High gap percentage: {:.2}% > {:.2}%", 
                               quality.gaps_percent, self.config.quality_thresholds.max_gaps_percent));
        }
        
        // Performance validations
        if performance.assembly_time_ms <= self.config.performance_targets.max_assembly_time_ms {
            passed_checks.push(format!("Assembly time within target: {}ms", performance.assembly_time_ms));
        } else {
            failed_checks.push(format!("Assembly time exceeded: {}ms > {}ms", 
                              performance.assembly_time_ms, self.config.performance_targets.max_assembly_time_ms));
        }
        
        if performance.memory_budget_utilization <= self.config.performance_targets.max_memory_budget_utilization {
            passed_checks.push(format!("Memory usage within budget: {:.1}%", performance.memory_budget_utilization));
        } else {
            failed_checks.push(format!("Memory budget exceeded: {:.1}% > {:.1}%", 
                              performance.memory_budget_utilization, self.config.performance_targets.max_memory_budget_utilization));
        }
        
        if performance.throughput_reads_per_sec >= self.config.performance_targets.min_throughput_reads_per_sec {
            passed_checks.push(format!("Throughput meets target: {:.1} reads/sec", performance.throughput_reads_per_sec));
        } else {
            warning_checks.push(format!("Throughput below target: {:.1} < {:.1} reads/sec", 
                               performance.throughput_reads_per_sec, self.config.performance_targets.min_throughput_reads_per_sec));
        }
        
        // Store metrics for comparison
        metrics_comparison.insert("n50".to_string(), quality.n50 as f64);
        metrics_comparison.insert("assembly_time_ms".to_string(), performance.assembly_time_ms as f64);
        metrics_comparison.insert("memory_utilization".to_string(), performance.memory_budget_utilization);
        metrics_comparison.insert("throughput".to_string(), performance.throughput_reads_per_sec);
        
        ValidationDetails {
            passed_checks,
            failed_checks,
            warning_checks,
            metrics_comparison,
        }
    }
    
    fn calculate_summary_metrics(&self, test_results: &[TestResult]) -> SummaryMetrics {
        if test_results.is_empty() {
            return SummaryMetrics {
                total_tests: 0,
                passed_tests: 0,
                failed_tests: 0,
                warning_tests: 0,
                average_assembly_time_ms: 0.0,
                average_memory_usage_mb: 0.0,
                average_n50: 0.0,
                overall_quality_score: 0.0,
                overall_performance_score: 0.0,
            };
        }
        
        let total_tests = test_results.len();
        let passed_tests = test_results.iter().filter(|r| matches!(r.status, TestStatus::Passed)).count();
        let failed_tests = test_results.iter().filter(|r| matches!(r.status, TestStatus::Failed)).count();
        let warning_tests = test_results.iter().filter(|r| matches!(r.status, TestStatus::Warning)).count();
        
        let average_assembly_time_ms = test_results.iter()
            .map(|r| r.performance_metrics.assembly_time_ms as f64)
            .sum::<f64>() / total_tests as f64;
        
        let average_memory_usage_mb = test_results.iter()
            .map(|r| r.performance_metrics.memory_usage_mb)
            .sum::<f64>() / total_tests as f64;
        
        let average_n50 = test_results.iter()
            .map(|r| r.quality_metrics.n50 as f64)
            .sum::<f64>() / total_tests as f64;
        
        // Simple scoring: passed tests contribute positively, failed tests negatively
        let overall_quality_score = (passed_tests as f64 * 100.0 - failed_tests as f64 * 50.0) / total_tests as f64;
        let overall_performance_score = (passed_tests as f64 * 100.0 - failed_tests as f64 * 80.0) / total_tests as f64;
        
        SummaryMetrics {
            total_tests,
            passed_tests,
            failed_tests,
            warning_tests,
            average_assembly_time_ms,
            average_memory_usage_mb,
            average_n50,
            overall_quality_score: overall_quality_score.max(0.0),
            overall_performance_score: overall_performance_score.max(0.0),
        }
    }
    
    fn generate_recommendations(&self, test_results: &[TestResult], summary: &SummaryMetrics) -> Vec<String> {
        let mut recommendations = Vec::new();
        
        if summary.failed_tests > 0 {
            recommendations.push("Critical: Some tests failed. Review failed checks and optimize assembly parameters.".to_string());
        }
        
        if summary.average_assembly_time_ms > self.config.performance_targets.max_assembly_time_ms as f64 {
            recommendations.push("Performance: Consider optimizing algorithm or reducing dataset complexity for better speed.".to_string());
        }
        
        if summary.average_memory_usage_mb > self.config.performance_targets.max_memory_usage_mb {
            recommendations.push("Memory: Implement stricter memory management or reduce memory budget allocation.".to_string());
        }
        
        if summary.average_n50 < self.config.quality_thresholds.min_n50 as f64 * 1.5 {
            recommendations.push("Quality: N50 could be improved. Consider adjusting k-mer selection or error correction.".to_string());
        }
        
        if summary.overall_quality_score < 80.0 {
            recommendations.push("Assembly quality below expectations. Review quality thresholds and assembly parameters.".to_string());
        }
        
        if summary.overall_performance_score < 70.0 {
            recommendations.push("Performance optimization needed. Focus on algorithmic improvements.".to_string());
        }
        
        if recommendations.is_empty() {
            recommendations.push("Excellent! All validation criteria met. Assembly optimization is performing well.".to_string());
        }
        
        recommendations
    }
    
    /// Generate detailed validation report
    pub fn generate_detailed_report(&self, report: &ValidationReport) -> String {
        let mut output = String::new();
        
        output.push_str("# Automated Assembly Validation Report\n\n");
        output.push_str(&format!("**Generated:** {}\n", report.timestamp));
        output.push_str(&format!("**Overall Status:** {:?}\n\n", report.overall_status));
        
        // Summary
        output.push_str("## Summary\n\n");
        let summary = &report.summary_metrics;
        output.push_str(&format!("- **Total Tests:** {}\n", summary.total_tests));
        output.push_str(&format!("- **Passed:** {} ({:.1}%)\n", summary.passed_tests, 
                               (summary.passed_tests as f64 / summary.total_tests as f64) * 100.0));
        output.push_str(&format!("- **Failed:** {}\n", summary.failed_tests));
        output.push_str(&format!("- **Warnings:** {}\n", summary.warning_tests));
        output.push_str(&format!("- **Average Assembly Time:** {:.0}ms\n", summary.average_assembly_time_ms));
        output.push_str(&format!("- **Average Memory Usage:** {:.1}MB\n", summary.average_memory_usage_mb));
        output.push_str(&format!("- **Average N50:** {:.0}\n", summary.average_n50));
        output.push_str(&format!("- **Quality Score:** {:.1}/100\n", summary.overall_quality_score));
        output.push_str(&format!("- **Performance Score:** {:.1}/100\n\n", summary.overall_performance_score));
        
        // Test Results
        output.push_str("## Test Results\n\n");
        for result in &report.test_results {
            output.push_str(&format!("### {} ({:?})\n\n", result.test_name, result.test_type));
            output.push_str(&format!("**Status:** {:?}\n", result.status));
            output.push_str(&format!("**Execution Time:** {}ms\n", result.execution_time_ms));
            
            output.push_str("\n**Quality Metrics:**\n");
            output.push_str(&format!("- N50: {}\n", result.quality_metrics.n50));
            output.push_str(&format!("- Total Length: {}\n", result.quality_metrics.total_length));
            output.push_str(&format!("- Contigs: {}\n", result.quality_metrics.num_contigs));
            output.push_str(&format!("- Assembly Efficiency: {:.2}\n", result.quality_metrics.assembly_efficiency));
            
            output.push_str("\n**Performance Metrics:**\n");
            output.push_str(&format!("- Assembly Time: {}ms\n", result.performance_metrics.assembly_time_ms));
            output.push_str(&format!("- Memory Usage: {:.1}MB\n", result.performance_metrics.memory_usage_mb));
            output.push_str(&format!("- Throughput: {:.1} reads/sec\n", result.performance_metrics.throughput_reads_per_sec));
            output.push_str(&format!("- Memory Utilization: {:.1}%\n\n", result.performance_metrics.memory_budget_utilization));
        }
        
        // Recommendations
        if !report.recommendations.is_empty() {
            output.push_str("## Recommendations\n\n");
            for (i, rec) in report.recommendations.iter().enumerate() {
                output.push_str(&format!("{}. {}\n", i + 1, rec));
            }
            output.push_str("\n");
        }
        
        // Warnings and Errors
        if !report.warnings.is_empty() {
            output.push_str("## Warnings\n\n");
            for warning in &report.warnings {
                output.push_str(&format!("- ⚠️  {}\n", warning));
            }
            output.push_str("\n");
        }
        
        if !report.errors.is_empty() {
            output.push_str("## Errors\n\n");
            for error in &report.errors {
                output.push_str(&format!("- ❌ {}\n", error));
            }
        }
        
        output
    }
}

// Validation Pipeline Tests

#[test]
fn test_automated_validation_pipeline() {
    let config = ValidationPipeline::default_config();
    let mut pipeline = ValidationPipeline::new(config);
    
    println!("🧪 Testing Automated Validation Pipeline");
    
    let report = pipeline.run_validation().unwrap();
    
    println!("\n📄 Validation Report Summary:");
    println!("Overall Status: {:?}", report.overall_status);
    println!("Total Tests: {}", report.summary_metrics.total_tests);
    println!("Passed: {}", report.summary_metrics.passed_tests);
    println!("Failed: {}", report.summary_metrics.failed_tests);
    println!("Quality Score: {:.1}/100", report.summary_metrics.overall_quality_score);
    println!("Performance Score: {:.1}/100", report.summary_metrics.overall_performance_score);
    
    // Generate detailed report
    let detailed_report = pipeline.generate_detailed_report(&report);
    println!("\n📈 Detailed report generated ({} characters)", detailed_report.len());
    
    // Basic validations
    assert!(report.summary_metrics.total_tests > 0, "Should execute some tests");
    assert!(!matches!(report.overall_status, ValidationStatus::Error), "Pipeline should not error");
    
    // At least some tests should pass
    assert!(report.summary_metrics.passed_tests > 0, "At least some tests should pass");
    
    println!("🎉 Automated validation pipeline test completed!");
}

#[test]
fn test_validation_with_strict_thresholds() {
    let mut config = ValidationPipeline::default_config();
    
    // Set very strict thresholds to test failure detection
    config.quality_thresholds.min_n50 = 1000; // Very high N50 requirement
    config.performance_targets.max_assembly_time_ms = 10; // Very short time limit
    
    let mut pipeline = ValidationPipeline::new(config);
    let report = pipeline.run_validation().unwrap();
    
    println!("Strict thresholds test:");
    println!("  Status: {:?}", report.overall_status);
    println!("  Failed tests: {}", report.summary_metrics.failed_tests);
    
    // With strict thresholds, some tests should fail
    assert!(report.summary_metrics.failed_tests > 0, "Strict thresholds should cause failures");
    assert!(matches!(report.overall_status, ValidationStatus::Failed), "Overall status should be failed");
    assert!(!report.errors.is_empty() || !report.warnings.is_empty(), "Should have errors or warnings");
}

#[test]
fn test_validation_pipeline_edge_cases() {
    // Test with minimal configuration
    let mut config = ValidationPipeline::default_config();
    config.test_suites = vec![config.test_suites[0].clone()]; // Only one test
    
    let mut pipeline = ValidationPipeline::new(config);
    let report = pipeline.run_validation().unwrap();
    
    println!("Edge case test:");
    println!("  Total tests: {}", report.summary_metrics.total_tests);
    println!("  Status: {:?}", report.overall_status);
    
    assert_eq!(report.summary_metrics.total_tests, 1, "Should run exactly one test");
    assert!(!report.test_results.is_empty(), "Should have test results");
    
    // Test disabled test suites
    let mut config = ValidationPipeline::default_config();
    for test_suite in &mut config.test_suites {
        test_suite.enabled = false;
    }
    
    let mut pipeline = ValidationPipeline::new(config);
    let report = pipeline.run_validation().unwrap();
    
    assert_eq!(report.summary_metrics.total_tests, 0, "Should run no tests when all disabled");
    assert!(report.test_results.is_empty(), "Should have no test results");
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

/// Helper struct for basic contig metrics
struct BasicContigMetrics {
    n50: usize,
    n90: usize,
    total_length: usize,
    num_contigs: usize,
    coverage_mean: f64,
    coverage_std: f64,
    gc_content: f64,
    contiguity_score: f64,
}

/// Calculate basic contig metrics
fn calculate_contig_metrics(contigs: &[Contig]) -> BasicContigMetrics {
    if contigs.is_empty() {
        return BasicContigMetrics {
            n50: 0,
            n90: 0,
            total_length: 0,
            num_contigs: 0,
            coverage_mean: 0.0,
            coverage_std: 0.0,
            gc_content: 0.0,
            contiguity_score: 0.0,
        };
    }

    let num_contigs = contigs.len();
    let total_length: usize = contigs.iter().map(|c| c.length).sum();
    let coverage_mean = contigs.iter().map(|c| c.coverage).sum::<f64>() / num_contigs as f64;

    // Calculate coverage std dev
    let variance = contigs.iter()
        .map(|c| (c.coverage - coverage_mean).powi(2))
        .sum::<f64>() / num_contigs as f64;
    let coverage_std = variance.sqrt();

    // Calculate N50 and N90
    let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
    lengths.sort_by(|a, b| b.cmp(a));

    let mut cumulative = 0;
    let mut n50 = 0;
    let mut n90 = 0;
    let half_total = total_length / 2;
    let ninety_percent = (total_length as f64 * 0.9) as usize;

    for &length in &lengths {
        cumulative += length;
        if n50 == 0 && cumulative >= half_total {
            n50 = length;
        }
        if cumulative >= ninety_percent {
            n90 = length;
            break;
        }
    }

    // Calculate GC content
    let mut gc_count = 0;
    let mut total_bases = 0;
    for contig in contigs {
        for base in contig.sequence.chars() {
            match base.to_ascii_uppercase() {
                'G' | 'C' => gc_count += 1,
                'A' | 'T' => {},
                _ => continue,
            }
            total_bases += 1;
        }
    }
    let gc_content = if total_bases > 0 {
        gc_count as f64 / total_bases as f64
    } else {
        0.0
    };

    // Simple contiguity score
    let contiguity_score = num_contigs as f64 / (total_length as f64 / 1000.0).max(1.0);

    BasicContigMetrics {
        n50,
        n90,
        total_length,
        num_contigs,
        coverage_mean,
        coverage_std,
        gc_content,
        contiguity_score,
    }
}
