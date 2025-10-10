//! Performance Bottleneck Analysis Framework
//! ==========================================
//!
//! Comprehensive performance analysis and bottleneck detection for metagenomic assembly.
//! Provides real-time monitoring, resource profiling, and optimization recommendations.

use anyhow::Result;
use serde::{Serialize, Deserialize};
use std::collections::{HashMap, VecDeque};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

/// Performance analysis configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnalysisConfig {
    /// Sampling interval for metrics collection (ms)
    pub sampling_interval_ms: u64,
    /// Maximum number of samples to retain
    pub max_samples: usize,
    /// Memory usage threshold for warnings (MB)
    pub memory_warning_threshold_mb: usize,
    /// CPU usage threshold for warnings (%)
    pub cpu_warning_threshold_percent: f64,
    /// I/O wait threshold for warnings (ms)
    pub io_wait_threshold_ms: u64,
    /// Enable detailed tracing
    pub enable_detailed_tracing: bool,
}

impl Default for AnalysisConfig {
    fn default() -> Self {
        Self {
            sampling_interval_ms: 100,
            max_samples: 1000,
            memory_warning_threshold_mb: 1024,
            cpu_warning_threshold_percent: 80.0,
            io_wait_threshold_ms: 100,
            enable_detailed_tracing: false,
        }
    }
}

/// Real-time performance metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerformanceMetrics {
    pub timestamp: u64,
    pub memory_usage_mb: f64,
    pub cpu_usage_percent: f64,
    pub io_wait_ms: u64,
    pub active_threads: usize,
    pub heap_allocations: u64,
    pub task_queue_depth: usize,
    pub pipeline_stage: String,
    pub custom_metrics: HashMap<String, f64>,
}

/// Bottleneck detection result
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BottleneckAnalysis {
    pub bottleneck_type: BottleneckType,
    pub severity: Severity,
    pub impact_score: f64,
    pub description: String,
    pub root_cause: String,
    pub optimization_suggestions: Vec<OptimizationSuggestion>,
    pub estimated_improvement: f64,
    pub implementation_difficulty: Difficulty,
    pub affected_stages: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BottleneckType {
    Memory,
    CPU,
    IO,
    Coordination,
    Algorithm,
    DataStructure,
    NetworkIO,
    DiskIO,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Severity {
    Low,
    Medium,
    High,
    Critical,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Difficulty {
    Easy,
    Medium,
    Hard,
    Complex,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationSuggestion {
    pub title: String,
    pub description: String,
    pub expected_improvement: f64,
    pub implementation_effort: Difficulty,
    pub code_changes_required: Vec<String>,
    pub configuration_changes: HashMap<String, String>,
}

/// Performance profiler for real-time monitoring
pub struct PerformanceProfiler {
    config: AnalysisConfig,
    metrics_history: Arc<Mutex<VecDeque<PerformanceMetrics>>>,
    start_time: Instant,
    memory_baseline: AtomicU64,
    cpu_baseline: AtomicU64,
    task_counters: Arc<Mutex<HashMap<String, TaskCounter>>>,
    bottlenecks: Arc<Mutex<Vec<BottleneckAnalysis>>>,
}

#[derive(Debug)]
struct TaskCounter {
    total_executions: AtomicU64,
    total_duration_ms: AtomicU64,
    active_count: AtomicUsize,
    max_duration_ms: AtomicU64,
    min_duration_ms: AtomicU64,
}

impl PerformanceProfiler {
    pub fn new(config: AnalysisConfig) -> Self {
        Self {
            config,
            metrics_history: Arc::new(Mutex::new(VecDeque::new())),
            start_time: Instant::now(),
            memory_baseline: AtomicU64::new(0),
            cpu_baseline: AtomicU64::new(0),
            task_counters: Arc::new(Mutex::new(HashMap::new())),
            bottlenecks: Arc::new(Mutex::new(Vec::new())),
        }
    }

    /// Start monitoring performance metrics
    pub fn start_monitoring(&self) -> Result<()> {
        // Initialize baseline metrics
        let initial_metrics = self.collect_system_metrics("initialization")?;
        self.memory_baseline.store(initial_metrics.memory_usage_mb as u64, Ordering::Relaxed);
        self.cpu_baseline.store(initial_metrics.cpu_usage_percent as u64, Ordering::Relaxed);

        // Store initial metrics
        {
            let mut history = self.metrics_history.lock().unwrap();
            history.push_back(initial_metrics);
        }

        println!("🔍 Performance monitoring started");
        println!("   📊 Baseline memory: {:.1} MB", self.memory_baseline.load(Ordering::Relaxed));
        println!("   ⚡ Sampling interval: {} ms", self.config.sampling_interval_ms);

        Ok(())
    }

    /// Collect current system metrics
    pub fn collect_system_metrics(&self, stage: &str) -> Result<PerformanceMetrics> {
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)?
            .as_millis() as u64;

        // Memory usage estimation (simplified for cross-platform compatibility)
        let memory_usage_mb = self.estimate_memory_usage()?;

        // CPU usage estimation
        let cpu_usage_percent = self.estimate_cpu_usage()?;

        // Thread count
        let active_threads = self.count_active_threads();

        // Task queue depth
        let task_queue_depth = {
            let counters = self.task_counters.lock().unwrap();
            counters.values()
                .map(|c| c.active_count.load(Ordering::Relaxed))
                .sum()
        };

        let metrics = PerformanceMetrics {
            timestamp,
            memory_usage_mb,
            cpu_usage_percent,
            io_wait_ms: 0, // Simplified for now
            active_threads,
            heap_allocations: 0, // Would need specialized profiling
            task_queue_depth,
            pipeline_stage: stage.to_string(),
            custom_metrics: HashMap::new(),
        };

        // Store metrics with size limit
        {
            let mut history = self.metrics_history.lock().unwrap();
            history.push_back(metrics.clone());

            if history.len() > self.config.max_samples {
                history.pop_front();
            }
        }

        Ok(metrics)
    }

    /// Profile a specific task execution
    pub fn profile_task<F, R>(&self, task_name: &str, task: F) -> Result<(R, Duration)>
    where
        F: FnOnce() -> Result<R>,
    {
        let start_time = Instant::now();

        // Update task counter
        {
            let mut counters = self.task_counters.lock().unwrap();
            let counter = counters.entry(task_name.to_string())
                .or_insert_with(|| TaskCounter {
                    total_executions: AtomicU64::new(0),
                    total_duration_ms: AtomicU64::new(0),
                    active_count: AtomicUsize::new(0),
                    max_duration_ms: AtomicU64::new(0),
                    min_duration_ms: AtomicU64::new(u64::MAX),
                });

            counter.active_count.fetch_add(1, Ordering::Relaxed);
        }

        // Execute task
        let result = task();
        let duration = start_time.elapsed();
        let duration_ms = duration.as_millis() as u64;

        // Update counters after execution
        {
            let mut counters = self.task_counters.lock().unwrap();
            if let Some(counter) = counters.get_mut(task_name) {
                counter.active_count.fetch_sub(1, Ordering::Relaxed);
                counter.total_executions.fetch_add(1, Ordering::Relaxed);
                counter.total_duration_ms.fetch_add(duration_ms, Ordering::Relaxed);

                // Update min/max
                counter.max_duration_ms.fetch_max(duration_ms, Ordering::Relaxed);
                counter.min_duration_ms.fetch_min(duration_ms, Ordering::Relaxed);
            }
        }

        Ok((result?, duration))
    }

    /// Analyze performance data for bottlenecks
    pub fn analyze_bottlenecks(&self) -> Result<Vec<BottleneckAnalysis>> {
        let history = self.metrics_history.lock().unwrap();
        let counters = self.task_counters.lock().unwrap();

        if history.is_empty() {
            return Ok(Vec::new());
        }

        let mut bottlenecks = Vec::new();

        // Memory analysis
        if let Some(memory_bottleneck) = self.analyze_memory_bottleneck(&history)? {
            bottlenecks.push(memory_bottleneck);
        }

        // CPU analysis
        if let Some(cpu_bottleneck) = self.analyze_cpu_bottleneck(&history)? {
            bottlenecks.push(cpu_bottleneck);
        }

        // Task execution analysis
        bottlenecks.extend(self.analyze_task_bottlenecks(&counters)?);

        // Coordination overhead analysis
        if let Some(coordination_bottleneck) = self.analyze_coordination_overhead(&history)? {
            bottlenecks.push(coordination_bottleneck);
        }

        // Sort by impact score
        bottlenecks.sort_by(|a, b| b.impact_score.partial_cmp(&a.impact_score).unwrap());

        // Store results
        {
            let mut stored_bottlenecks = self.bottlenecks.lock().unwrap();
            *stored_bottlenecks = bottlenecks.clone();
        }

        Ok(bottlenecks)
    }

    /// Generate optimization report
    pub fn generate_optimization_report(&self) -> Result<OptimizationReport> {
        let bottlenecks = self.analyze_bottlenecks()?;
        let metrics_summary = self.generate_metrics_summary()?;
        let task_performance = self.generate_task_performance_summary()?;

        let report = OptimizationReport {
            timestamp: SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs(),
            overall_performance_score: self.calculate_performance_score(&metrics_summary),
            critical_bottlenecks: bottlenecks.into_iter()
                .filter(|b| matches!(b.severity, Severity::Critical | Severity::High))
                .collect(),
            metrics_summary,
            task_performance,
            optimization_roadmap: self.generate_optimization_roadmap()?,
            recommendations: self.generate_recommendations()?,
        };

        Ok(report)
    }

    // Private analysis methods

    fn analyze_memory_bottleneck(&self, history: &VecDeque<PerformanceMetrics>) -> Result<Option<BottleneckAnalysis>> {
        let recent_metrics: Vec<_> = history.iter().rev().take(10).collect();
        let avg_memory = recent_metrics.iter()
            .map(|m| m.memory_usage_mb)
            .sum::<f64>() / recent_metrics.len() as f64;

        let baseline = self.memory_baseline.load(Ordering::Relaxed) as f64;
        let memory_growth = ((avg_memory - baseline) / baseline * 100.0).max(0.0);

        if avg_memory > self.config.memory_warning_threshold_mb as f64 || memory_growth > 200.0 {
            let severity = if memory_growth > 500.0 {
                Severity::Critical
            } else if memory_growth > 300.0 {
                Severity::High
            } else {
                Severity::Medium
            };

            let suggestions = vec![
                OptimizationSuggestion {
                    title: "Implement streaming processing".to_string(),
                    description: "Process data in chunks to reduce peak memory usage".to_string(),
                    expected_improvement: 40.0,
                    implementation_effort: Difficulty::Medium,
                    code_changes_required: vec![
                        "Modify k-mer counting to use bounded counters".to_string(),
                        "Implement chunk-based graph construction".to_string(),
                    ],
                    configuration_changes: HashMap::from([
                        ("chunk_size".to_string(), "1000".to_string()),
                        ("memory_budget_mb".to_string(), "2048".to_string()),
                    ]),
                },
                OptimizationSuggestion {
                    title: "Use memory-efficient data structures".to_string(),
                    description: "Replace hash maps with compact alternatives where possible".to_string(),
                    expected_improvement: 25.0,
                    implementation_effort: Difficulty::Hard,
                    code_changes_required: vec![
                        "Replace AHashMap with Vec + binary search for small collections".to_string(),
                        "Use bit-packed k-mer representations".to_string(),
                    ],
                    configuration_changes: HashMap::new(),
                },
            ];

            return Ok(Some(BottleneckAnalysis {
                bottleneck_type: BottleneckType::Memory,
                severity,
                impact_score: memory_growth.min(100.0),
                description: format!("Memory usage increased by {:.1}% above baseline", memory_growth),
                root_cause: "Large k-mer tables and graph structures consuming excessive memory".to_string(),
                optimization_suggestions: suggestions,
                estimated_improvement: 45.0,
                implementation_difficulty: Difficulty::Medium,
                affected_stages: vec!["k-mer counting".to_string(), "graph construction".to_string()],
            }));
        }

        Ok(None)
    }

    fn analyze_cpu_bottleneck(&self, history: &VecDeque<PerformanceMetrics>) -> Result<Option<BottleneckAnalysis>> {
        let recent_metrics: Vec<_> = history.iter().rev().take(20).collect();
        let avg_cpu = recent_metrics.iter()
            .map(|m| m.cpu_usage_percent)
            .sum::<f64>() / recent_metrics.len() as f64;

        if avg_cpu > self.config.cpu_warning_threshold_percent {
            let severity = if avg_cpu > 95.0 {
                Severity::Critical
            } else if avg_cpu > 85.0 {
                Severity::High
            } else {
                Severity::Medium
            };

            let suggestions = vec![
                OptimizationSuggestion {
                    title: "Increase parallelization".to_string(),
                    description: "Better utilize multiple CPU cores for k-mer processing".to_string(),
                    expected_improvement: 30.0,
                    implementation_effort: Difficulty::Medium,
                    code_changes_required: vec![
                        "Use rayon parallel iterators for k-mer extraction".to_string(),
                        "Implement parallel graph construction".to_string(),
                    ],
                    configuration_changes: HashMap::from([
                        ("cpu_cores".to_string(), "auto".to_string()),
                        ("parallel_processing".to_string(), "true".to_string()),
                    ]),
                },
                OptimizationSuggestion {
                    title: "Optimize hot paths".to_string(),
                    description: "Profile and optimize frequently called functions".to_string(),
                    expected_improvement: 20.0,
                    implementation_effort: Difficulty::Hard,
                    code_changes_required: vec![
                        "Cache hash calculations".to_string(),
                        "Use SIMD for nucleotide processing".to_string(),
                    ],
                    configuration_changes: HashMap::new(),
                },
            ];

            return Ok(Some(BottleneckAnalysis {
                bottleneck_type: BottleneckType::CPU,
                severity,
                impact_score: (avg_cpu - self.config.cpu_warning_threshold_percent).max(0.0),
                description: format!("CPU usage at {:.1}%, exceeding threshold", avg_cpu),
                root_cause: "Sequential processing and inefficient algorithms causing CPU saturation".to_string(),
                optimization_suggestions: suggestions,
                estimated_improvement: 35.0,
                implementation_difficulty: Difficulty::Medium,
                affected_stages: vec!["k-mer processing".to_string(), "graph construction".to_string()],
            }));
        }

        Ok(None)
    }

    fn analyze_task_bottlenecks(&self, counters: &HashMap<String, TaskCounter>) -> Result<Vec<BottleneckAnalysis>> {
        let mut bottlenecks = Vec::new();

        for (task_name, counter) in counters {
            let total_executions = counter.total_executions.load(Ordering::Relaxed);
            let total_duration_ms = counter.total_duration_ms.load(Ordering::Relaxed);
            let max_duration_ms = counter.max_duration_ms.load(Ordering::Relaxed);
            let min_duration_ms = counter.min_duration_ms.load(Ordering::Relaxed);

            if total_executions == 0 {
                continue;
            }

            let avg_duration_ms = total_duration_ms as f64 / total_executions as f64;
            let variance = if min_duration_ms != u64::MAX {
                (max_duration_ms as f64 - min_duration_ms as f64) / avg_duration_ms
            } else {
                0.0
            };

            // Detect tasks with high execution time or high variance
            if avg_duration_ms > 1000.0 || variance > 2.0 {
                let severity = if avg_duration_ms > 10000.0 || variance > 5.0 {
                    Severity::High
                } else {
                    Severity::Medium
                };

                let suggestions = if task_name.contains("graph") {
                    vec![
                        OptimizationSuggestion {
                            title: "Optimize graph construction".to_string(),
                            description: "Use more efficient graph data structures".to_string(),
                            expected_improvement: 50.0,
                            implementation_effort: Difficulty::Medium,
                            code_changes_required: vec![
                                "Replace adjacency lists with edge arrays".to_string(),
                                "Batch edge insertions".to_string(),
                            ],
                            configuration_changes: HashMap::new(),
                        },
                    ]
                } else if task_name.contains("kmer") {
                    vec![
                        OptimizationSuggestion {
                            title: "Optimize k-mer processing".to_string(),
                            description: "Use rolling hash for k-mer generation".to_string(),
                            expected_improvement: 40.0,
                            implementation_effort: Difficulty::Easy,
                            code_changes_required: vec![
                                "Implement rolling hash for consecutive k-mers".to_string(),
                            ],
                            configuration_changes: HashMap::new(),
                        },
                    ]
                } else {
                    Vec::new()
                };

                bottlenecks.push(BottleneckAnalysis {
                    bottleneck_type: BottleneckType::Algorithm,
                    severity,
                    impact_score: (avg_duration_ms / 100.0).min(100.0),
                    description: format!("Task '{}' has high execution time ({:.1}ms avg) and variance ({:.1}x)",
                                        task_name, avg_duration_ms, variance),
                    root_cause: format!("Inefficient algorithm implementation in {}", task_name),
                    optimization_suggestions: suggestions,
                    estimated_improvement: 40.0,
                    implementation_difficulty: Difficulty::Medium,
                    affected_stages: vec![task_name.clone()],
                });
            }
        }

        Ok(bottlenecks)
    }

    fn analyze_coordination_overhead(&self, history: &VecDeque<PerformanceMetrics>) -> Result<Option<BottleneckAnalysis>> {
        let recent_metrics: Vec<_> = history.iter().rev().take(10).collect();
        let avg_queue_depth = recent_metrics.iter()
            .map(|m| m.task_queue_depth as f64)
            .sum::<f64>() / recent_metrics.len() as f64;

        if avg_queue_depth > 100.0 {
            let suggestions = vec![
                OptimizationSuggestion {
                    title: "Reduce coordination overhead".to_string(),
                    description: "Minimize inter-thread communication and task queuing".to_string(),
                    expected_improvement: 25.0,
                    implementation_effort: Difficulty::Medium,
                    code_changes_required: vec![
                        "Batch multiple operations per thread".to_string(),
                        "Use lock-free data structures".to_string(),
                    ],
                    configuration_changes: HashMap::from([
                        ("batch_size".to_string(), "100".to_string()),
                    ]),
                },
            ];

            return Ok(Some(BottleneckAnalysis {
                bottleneck_type: BottleneckType::Coordination,
                severity: Severity::Medium,
                impact_score: (avg_queue_depth / 10.0).min(100.0),
                description: format!("High task queue depth ({:.0} tasks)", avg_queue_depth),
                root_cause: "Excessive coordination overhead between threads".to_string(),
                optimization_suggestions: suggestions,
                estimated_improvement: 25.0,
                implementation_difficulty: Difficulty::Medium,
                affected_stages: vec!["parallel processing".to_string()],
            }));
        }

        Ok(None)
    }

    // Helper methods for metrics collection

    fn estimate_memory_usage(&self) -> Result<f64> {
        // Simplified memory estimation - in a real implementation,
        // this would use platform-specific APIs
        match std::env::var("MOCK_MEMORY_MB") {
            Ok(mem_str) => Ok(mem_str.parse().unwrap_or(512.0)),
            Err(_) => {
                // Fallback estimation based on runtime
                let elapsed_secs = self.start_time.elapsed().as_secs() as f64;
                let base_memory = 100.0;
                let growth_rate = 2.0; // MB per second (simulated)
                Ok(base_memory + elapsed_secs * growth_rate)
            }
        }
    }

    fn estimate_cpu_usage(&self) -> Result<f64> {
        // Simplified CPU estimation
        match std::env::var("MOCK_CPU_PERCENT") {
            Ok(cpu_str) => Ok(cpu_str.parse().unwrap_or(50.0)),
            Err(_) => {
                // Simulate CPU usage based on active tasks
                let counters = self.task_counters.lock().unwrap();
                let active_tasks: usize = counters.values()
                    .map(|c| c.active_count.load(Ordering::Relaxed))
                    .sum();
                Ok((active_tasks as f64 * 15.0).min(100.0))
            }
        }
    }

    fn count_active_threads(&self) -> usize {
        rayon::current_num_threads().max(1)
    }

    fn generate_metrics_summary(&self) -> Result<MetricsSummary> {
        let history = self.metrics_history.lock().unwrap();

        if history.is_empty() {
            return Ok(MetricsSummary::default());
        }

        let memory_values: Vec<f64> = history.iter().map(|m| m.memory_usage_mb).collect();
        let cpu_values: Vec<f64> = history.iter().map(|m| m.cpu_usage_percent).collect();

        Ok(MetricsSummary {
            memory_usage: StatsSummary {
                min: memory_values.iter().copied().fold(f64::INFINITY, f64::min),
                max: memory_values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                avg: memory_values.iter().sum::<f64>() / memory_values.len() as f64,
                median: Self::calculate_median(&memory_values),
            },
            cpu_usage: StatsSummary {
                min: cpu_values.iter().copied().fold(f64::INFINITY, f64::min),
                max: cpu_values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                avg: cpu_values.iter().sum::<f64>() / cpu_values.len() as f64,
                median: Self::calculate_median(&cpu_values),
            },
            total_samples: history.len(),
            monitoring_duration_secs: self.start_time.elapsed().as_secs(),
        })
    }

    fn generate_task_performance_summary(&self) -> Result<Vec<TaskPerformanceSummary>> {
        let counters = self.task_counters.lock().unwrap();
        let mut summaries = Vec::new();

        for (task_name, counter) in counters.iter() {
            let total_executions = counter.total_executions.load(Ordering::Relaxed);
            let total_duration_ms = counter.total_duration_ms.load(Ordering::Relaxed);

            if total_executions > 0 {
                summaries.push(TaskPerformanceSummary {
                    task_name: task_name.clone(),
                    total_executions,
                    total_duration_ms,
                    avg_duration_ms: total_duration_ms as f64 / total_executions as f64,
                    max_duration_ms: counter.max_duration_ms.load(Ordering::Relaxed),
                    min_duration_ms: if counter.min_duration_ms.load(Ordering::Relaxed) == u64::MAX {
                        0
                    } else {
                        counter.min_duration_ms.load(Ordering::Relaxed)
                    },
                });
            }
        }

        summaries.sort_by(|a, b| b.total_duration_ms.cmp(&a.total_duration_ms));
        Ok(summaries)
    }

    fn calculate_performance_score(&self, summary: &MetricsSummary) -> f64 {
        // Simple scoring based on resource efficiency
        let memory_score = (100.0 - (summary.memory_usage.avg / 10.0)).max(0.0);
        let cpu_score = (100.0 - summary.cpu_usage.avg).max(0.0);

        (memory_score + cpu_score) / 2.0
    }

    fn generate_optimization_roadmap(&self) -> Result<Vec<OptimizationPhase>> {
        Ok(vec![
            OptimizationPhase {
                phase_name: "Quick Wins".to_string(),
                estimated_duration_days: 2,
                expected_improvement: 20.0,
                optimizations: vec![
                    "Enable parallel processing".to_string(),
                    "Adjust chunk sizes".to_string(),
                    "Optimize k-mer size selection".to_string(),
                ],
            },
            OptimizationPhase {
                phase_name: "Memory Optimization".to_string(),
                estimated_duration_days: 7,
                expected_improvement: 40.0,
                optimizations: vec![
                    "Implement streaming processing".to_string(),
                    "Use memory-efficient data structures".to_string(),
                    "Add memory cleanup routines".to_string(),
                ],
            },
            OptimizationPhase {
                phase_name: "Algorithm Improvements".to_string(),
                estimated_duration_days: 14,
                expected_improvement: 60.0,
                optimizations: vec![
                    "Optimize graph construction algorithms".to_string(),
                    "Implement rolling hash for k-mers".to_string(),
                    "Add SIMD optimizations".to_string(),
                ],
            },
        ])
    }

    fn generate_recommendations(&self) -> Result<Vec<String>> {
        Ok(vec![
            "Monitor memory usage during peak processing phases".to_string(),
            "Use profiling tools to identify hot code paths".to_string(),
            "Test with different k-mer sizes to find optimal balance".to_string(),
            "Consider SSD storage for better I/O performance".to_string(),
            "Validate optimizations with realistic datasets".to_string(),
        ])
    }

    fn calculate_median(values: &[f64]) -> f64 {
        if values.is_empty() {
            return 0.0;
        }

        let mut sorted = values.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            (sorted[mid - 1] + sorted[mid]) / 2.0
        } else {
            sorted[mid]
        }
    }
}

// Supporting data structures for reports

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationReport {
    pub timestamp: u64,
    pub overall_performance_score: f64,
    pub critical_bottlenecks: Vec<BottleneckAnalysis>,
    pub metrics_summary: MetricsSummary,
    pub task_performance: Vec<TaskPerformanceSummary>,
    pub optimization_roadmap: Vec<OptimizationPhase>,
    pub recommendations: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MetricsSummary {
    pub memory_usage: StatsSummary,
    pub cpu_usage: StatsSummary,
    pub total_samples: usize,
    pub monitoring_duration_secs: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct StatsSummary {
    pub min: f64,
    pub max: f64,
    pub avg: f64,
    pub median: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaskPerformanceSummary {
    pub task_name: String,
    pub total_executions: u64,
    pub total_duration_ms: u64,
    pub avg_duration_ms: f64,
    pub max_duration_ms: u64,
    pub min_duration_ms: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationPhase {
    pub phase_name: String,
    pub estimated_duration_days: u32,
    pub expected_improvement: f64,
    pub optimizations: Vec<String>,
}

/// Helper macros for performance profiling

#[macro_export]
macro_rules! profile_task {
    ($profiler:expr, $task_name:expr, $task:expr) => {
        $profiler.profile_task($task_name, || Ok($task))
    };
}

#[macro_export]
macro_rules! collect_metrics {
    ($profiler:expr, $stage:expr) => {
        $profiler.collect_system_metrics($stage)
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    use std::time::Duration as StdDuration;

    #[test]
    fn test_performance_profiler_creation() {
        let config = AnalysisConfig::default();
        let profiler = PerformanceProfiler::new(config);

        assert!(profiler.start_monitoring().is_ok());
    }

    #[test]
    fn test_metrics_collection() {
        let config = AnalysisConfig::default();
        let profiler = PerformanceProfiler::new(config);

        profiler.start_monitoring().unwrap();

        let metrics = profiler.collect_system_metrics("test_stage").unwrap();
        assert!(!metrics.pipeline_stage.is_empty());
        assert!(metrics.timestamp > 0);
    }

    #[test]
    fn test_task_profiling() {
        let config = AnalysisConfig::default();
        let profiler = PerformanceProfiler::new(config);

        let (result, duration) = profiler.profile_task("test_task", || {
            thread::sleep(StdDuration::from_millis(10));
            Ok(42)
        }).unwrap();

        assert_eq!(result, 42);
        assert!(duration.as_millis() >= 10);
    }

    #[test]
    fn test_bottleneck_analysis() {
        std::env::set_var("MOCK_MEMORY_MB", "2048"); // High memory usage
        std::env::set_var("MOCK_CPU_PERCENT", "90");  // High CPU usage

        let config = AnalysisConfig {
            memory_warning_threshold_mb: 1024,
            cpu_warning_threshold_percent: 80.0,
            ..Default::default()
        };
        let profiler = PerformanceProfiler::new(config);

        profiler.start_monitoring().unwrap();

        // Collect some metrics
        for _ in 0..5 {
            profiler.collect_system_metrics("test").unwrap();
        }

        let bottlenecks = profiler.analyze_bottlenecks().unwrap();
        assert!(!bottlenecks.is_empty());

        // Should detect memory and CPU bottlenecks
        let has_memory_bottleneck = bottlenecks.iter()
            .any(|b| matches!(b.bottleneck_type, BottleneckType::Memory));
        let has_cpu_bottleneck = bottlenecks.iter()
            .any(|b| matches!(b.bottleneck_type, BottleneckType::CPU));

        assert!(has_memory_bottleneck || has_cpu_bottleneck);

        // Clean up
        std::env::remove_var("MOCK_MEMORY_MB");
        std::env::remove_var("MOCK_CPU_PERCENT");
    }

    #[test]
    fn test_optimization_report() {
        let config = AnalysisConfig::default();
        let profiler = PerformanceProfiler::new(config);

        profiler.start_monitoring().unwrap();
        profiler.collect_system_metrics("test").unwrap();

        let report = profiler.generate_optimization_report().unwrap();
        assert!(report.overall_performance_score >= 0.0);
        assert!(report.overall_performance_score <= 100.0);
        assert!(!report.optimization_roadmap.is_empty());
        assert!(!report.recommendations.is_empty());
    }
}