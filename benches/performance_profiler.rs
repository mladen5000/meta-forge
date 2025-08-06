use anyhow::Result;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::{Duration, Instant};

/// Comprehensive performance profiler for genomics pipeline
pub struct GenomicsProfiler {
    /// Performance metrics for different pipeline stages
    stage_metrics: HashMap<String, StageMetrics>,
    /// Overall pipeline timing
    pipeline_start: Option<Instant>,
    /// Memory usage tracking
    memory_tracker: MemoryTracker,
    /// Current active stage
    current_stage: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageMetrics {
    pub stage_name: String,
    pub total_duration: Duration,
    pub call_count: usize,
    pub avg_duration: Duration,
    pub min_duration: Duration,
    pub max_duration: Duration,
    pub memory_peak: usize,
    pub throughput_items_per_sec: f64,
    pub items_processed: usize,
}

pub struct MemoryTracker {
    peak_memory: usize,
    current_estimate: usize,
    samples: Vec<(Instant, usize)>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PerformanceReport {
    pub total_pipeline_time: Duration,
    pub stage_breakdown: Vec<StageMetrics>,
    pub bottlenecks: Vec<String>,
    pub memory_profile: MemoryProfile,
    pub recommendations: Vec<String>,
    pub efficiency_score: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct MemoryProfile {
    pub peak_memory_mb: f64,
    pub avg_memory_mb: f64,
    pub memory_efficiency: f64,
    pub gc_pressure_estimate: f64,
}

impl Default for GenomicsProfiler {
    fn default() -> Self {
        Self::new()
    }
}

impl GenomicsProfiler {
    pub fn new() -> Self {
        Self {
            stage_metrics: HashMap::new(),
            pipeline_start: None,
            memory_tracker: MemoryTracker::new(),
            current_stage: None,
        }
    }

    /// Start timing the entire pipeline
    pub fn start_pipeline(&mut self) {
        self.pipeline_start = Some(Instant::now());
        self.memory_tracker.sample_memory();
        println!("🚀 Pipeline profiling started");
    }

    /// Start timing a specific stage
    pub fn start_stage(&mut self, stage_name: &str) -> StageTimer {
        self.current_stage = Some(stage_name.to_string());
        self.memory_tracker.sample_memory();

        println!("⏱️  Starting stage: {stage_name}");
        StageTimer::new(stage_name.to_string(), Instant::now())
    }

    /// End timing a stage and record metrics
    pub fn end_stage(&mut self, timer: StageTimer, items_processed: usize) {
        let duration = timer.end_time.elapsed();
        let stage_name = timer.stage_name.clone();

        // Update or create stage metrics
        let entry = self
            .stage_metrics
            .entry(stage_name.clone())
            .or_insert(StageMetrics {
                stage_name: stage_name.clone(),
                total_duration: Duration::from_nanos(0),
                call_count: 0,
                avg_duration: Duration::from_nanos(0),
                min_duration: Duration::from_secs(u64::MAX),
                max_duration: Duration::from_nanos(0),
                memory_peak: 0,
                throughput_items_per_sec: 0.0,
                items_processed: 0,
            });

        // Update metrics
        entry.total_duration += duration;
        entry.call_count += 1;
        entry.avg_duration = entry.total_duration / entry.call_count as u32;
        entry.min_duration = entry.min_duration.min(duration);
        entry.max_duration = entry.max_duration.max(duration);
        entry.items_processed += items_processed;
        entry.memory_peak = entry.memory_peak.max(self.memory_tracker.current_estimate);

        // Calculate throughput
        if duration.as_secs_f64() > 0.0 {
            entry.throughput_items_per_sec = items_processed as f64 / duration.as_secs_f64();
        }

        self.current_stage = None;
        self.memory_tracker.sample_memory();

        println!(
            "✅ Completed stage: {} ({:.2}s, {} items, {:.0} items/sec)",
            stage_name,
            duration.as_secs_f64(),
            items_processed,
            entry.throughput_items_per_sec
        );
    }

    /// Profile a specific function call
    pub fn profile_function<F, R>(&mut self, stage_name: &str, items: usize, func: F) -> Result<R>
    where
        F: FnOnce() -> Result<R>,
    {
        let timer = self.start_stage(stage_name);
        let result = func()?;
        self.end_stage(timer, items);
        Ok(result)
    }

    /// End pipeline timing and generate report
    pub fn end_pipeline(&mut self) -> PerformanceReport {
        let total_time = self
            .pipeline_start
            .map(|start| start.elapsed())
            .unwrap_or(Duration::from_nanos(0));

        let memory_profile = self.memory_tracker.get_profile();
        let mut stage_breakdown: Vec<_> = self.stage_metrics.values().cloned().collect();
        stage_breakdown.sort_by(|a, b| b.total_duration.cmp(&a.total_duration));

        let bottlenecks = self.identify_bottlenecks(&stage_breakdown);
        let recommendations = self.generate_recommendations(&stage_breakdown, &memory_profile);
        let efficiency_score = self.calculate_efficiency_score(&stage_breakdown, total_time);

        println!("🏁 Pipeline completed in {:.2}s", total_time.as_secs_f64());

        PerformanceReport {
            total_pipeline_time: total_time,
            stage_breakdown,
            bottlenecks,
            memory_profile,
            recommendations,
            efficiency_score,
        }
    }

    fn identify_bottlenecks(&self, stages: &[StageMetrics]) -> Vec<String> {
        let mut bottlenecks = Vec::new();
        let total_time: Duration = stages.iter().map(|s| s.total_duration).sum();

        for stage in stages {
            let percentage = stage.total_duration.as_secs_f64() / total_time.as_secs_f64() * 100.0;
            if percentage > 25.0 {
                bottlenecks.push(format!(
                    "{} ({}% of total time)",
                    stage.stage_name, percentage as u32
                ));
            }
        }

        // Check for low throughput stages
        for stage in stages {
            if stage.throughput_items_per_sec < 100.0 && stage.items_processed > 1000 {
                bottlenecks.push(format!(
                    "{} (low throughput: {:.0} items/sec)",
                    stage.stage_name, stage.throughput_items_per_sec
                ));
            }
        }

        bottlenecks
    }

    fn generate_recommendations(
        &self,
        stages: &[StageMetrics],
        memory: &MemoryProfile,
    ) -> Vec<String> {
        let mut recommendations = Vec::new();

        // Memory recommendations
        if memory.peak_memory_mb > 8000.0 {
            recommendations.push("Consider enabling compression for intermediate data".to_string());
            recommendations
                .push("Implement streaming processing to reduce memory usage".to_string());
        }

        if memory.memory_efficiency < 0.7 {
            recommendations.push(
                "Memory usage patterns suggest fragmentation - consider custom allocators"
                    .to_string(),
            );
        }

        // Stage-specific recommendations
        for stage in stages {
            if stage.throughput_items_per_sec < 1000.0 && stage.stage_name.contains("kmer") {
                recommendations.push(format!(
                    "Optimize k-mer processing in {} - consider parallel hash tables",
                    stage.stage_name
                ));
            }

            if stage.max_duration.as_secs_f64() / stage.min_duration.as_secs_f64() > 10.0 {
                recommendations.push(format!(
                    "High variance in {} timing - check for batch size optimization",
                    stage.stage_name
                ));
            }

            if stage.stage_name.contains("assembly") && stage.throughput_items_per_sec < 100.0 {
                recommendations.push("Assembly stage is slow - consider adaptive k-mer sizing or graph simplification".to_string());
            }

            if stage.stage_name.contains("database") && stage.avg_duration.as_millis() > 100 {
                recommendations.push(
                    "Database queries are slow - check indexes and consider caching".to_string(),
                );
            }
        }

        // Parallelization recommendations
        let sequential_stages: Vec<_> = stages
            .iter()
            .filter(|s| s.throughput_items_per_sec < 500.0 && s.items_processed > 100)
            .collect();

        if !sequential_stages.is_empty() {
            recommendations
                .push("Consider parallelizing stages with low throughput using rayon".to_string());
        }

        recommendations
    }

    fn calculate_efficiency_score(&self, stages: &[StageMetrics], total_time: Duration) -> f64 {
        if stages.is_empty() || total_time.as_secs_f64() == 0.0 {
            return 0.0;
        }

        // Calculate based on throughput and memory efficiency
        let avg_throughput = stages
            .iter()
            .map(|s| s.throughput_items_per_sec)
            .sum::<f64>()
            / stages.len() as f64;

        let memory_efficiency = self.memory_tracker.get_profile().memory_efficiency;

        // Normalize and combine metrics (0-1 scale)
        let throughput_score = (avg_throughput.ln() / 10.0).min(1.0).max(0.0);
        let memory_score = memory_efficiency;
        let time_score = (1.0 / (total_time.as_secs_f64() / 60.0)).min(1.0); // Penalty for long runs

        (throughput_score * 0.4 + memory_score * 0.3 + time_score * 0.3)
            .min(1.0)
            .max(0.0)
    }

    /// Print a detailed performance report
    pub fn print_report(&self, report: &PerformanceReport) {
        println!("\n📊 ===== GENOMICS PIPELINE PERFORMANCE REPORT =====");
        println!(
            "Total Pipeline Time: {:.2} seconds",
            report.total_pipeline_time.as_secs_f64()
        );
        println!("Efficiency Score: {:.2}/1.0", report.efficiency_score);
        println!(
            "Peak Memory Usage: {:.1} MB",
            report.memory_profile.peak_memory_mb
        );

        println!("\n⏱️  Stage Breakdown:");
        for stage in &report.stage_breakdown {
            println!(
                "  {} - {:.2}s ({} calls, {:.0} items/sec)",
                stage.stage_name,
                stage.total_duration.as_secs_f64(),
                stage.call_count,
                stage.throughput_items_per_sec
            );
        }

        if !report.bottlenecks.is_empty() {
            println!("\n🚨 Bottlenecks Identified:");
            for bottleneck in &report.bottlenecks {
                println!("  - {bottleneck}");
            }
        }

        if !report.recommendations.is_empty() {
            println!("\n💡 Optimization Recommendations:");
            for (i, rec) in report.recommendations.iter().enumerate() {
                println!("  {}. {}", i + 1, rec);
            }
        }

        println!("\n📈 Memory Profile:");
        println!("  Peak: {:.1} MB", report.memory_profile.peak_memory_mb);
        println!("  Average: {:.1} MB", report.memory_profile.avg_memory_mb);
        println!(
            "  Efficiency: {:.2}",
            report.memory_profile.memory_efficiency
        );
    }

    /// Export detailed metrics to JSON
    pub fn export_metrics(&self, report: &PerformanceReport, path: &str) -> Result<()> {
        let json = serde_json::to_string_pretty(report)?;
        std::fs::write(path, json)?;
        println!("📄 Performance metrics exported to {path}");
        Ok(())
    }
}

pub struct StageTimer {
    stage_name: String,
    end_time: Instant,
}

impl StageTimer {
    fn new(stage_name: String, start_time: Instant) -> Self {
        Self {
            stage_name,
            end_time: start_time,
        }
    }
}

impl MemoryTracker {
    fn new() -> Self {
        Self {
            peak_memory: 0,
            current_estimate: 0,
            samples: Vec::new(),
        }
    }

    fn sample_memory(&mut self) {
        // Simple memory estimation (in real implementation, use system APIs)
        let estimated = Self::get_memory_usage();
        self.current_estimate = estimated;
        self.peak_memory = self.peak_memory.max(estimated);
        self.samples.push((Instant::now(), estimated));

        // Keep only recent samples
        if self.samples.len() > 1000 {
            self.samples.drain(0..500);
        }
    }

    fn get_memory_usage() -> usize {
        // Placeholder - use actual system memory APIs in production
        #[cfg(target_os = "linux")]
        {
            if let Ok(status) = std::fs::read_to_string("/proc/self/status") {
                for line in status.lines() {
                    if line.starts_with("VmRSS:") {
                        if let Some(kb_str) = line.split_whitespace().nth(1) {
                            if let Ok(kb) = kb_str.parse::<usize>() {
                                return kb * 1024; // Convert to bytes
                            }
                        }
                    }
                }
            }
        }

        // Fallback estimation based on allocations
        std::env::var("MEMORY_ESTIMATE")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(100 * 1024 * 1024) // 100 MB default
    }

    fn get_profile(&self) -> MemoryProfile {
        let peak_mb = self.peak_memory as f64 / (1024.0 * 1024.0);
        let avg_mb = if !self.samples.is_empty() {
            self.samples.iter().map(|(_, mem)| *mem as f64).sum::<f64>()
                / self.samples.len() as f64
                / (1024.0 * 1024.0)
        } else {
            0.0
        };
        let memory_efficiency = if peak_mb > 0.0 { avg_mb / peak_mb } else { 0.0 };
        let gc_pressure_estimate = if self.samples.len() > 1 {
            let total_duration = self.samples.last().unwrap().0 - self.samples.first().unwrap().0;
            (self.samples.len() as f64 / total_duration.as_secs_f64()).max(0.0)
        } else {
            0.0
        };

        MemoryProfile {
            peak_memory_mb: peak_mb,
            avg_memory_mb: avg_mb,
            memory_efficiency,
            gc_pressure_estimate,
        }
    }
}
