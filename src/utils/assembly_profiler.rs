//! Assembly Performance Profiler
//!
//! Fine-grained profiling for assembly pipeline phases:
//! - K-mer counting
//! - Graph construction
//! - Contig generation
//! - Coverage filtering

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::{Duration, Instant};

/// Phase timing information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhaseMetrics {
    pub phase_name: String,
    pub duration_ms: u64,
    pub memory_mb: f64,
    pub throughput: Option<f64>, // items/second
    pub metadata: HashMap<String, String>,
}

/// Complete assembly profiling report
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyProfile {
    pub total_duration_ms: u64,
    pub peak_memory_mb: f64,
    pub phases: Vec<PhaseMetrics>,
    pub summary: ProfileSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProfileSummary {
    pub reads_processed: usize,
    pub kmers_counted: usize,
    pub nodes_created: usize,
    pub edges_created: usize,
    pub contigs_generated: usize,
    pub kmers_per_second: f64,
    pub reads_per_second: f64,
}

/// Real-time profiler for assembly operations
pub struct AssemblyProfiler {
    start_time: Instant,
    phases: Vec<PhaseMetrics>,
    current_phase: Option<(String, Instant)>,
    peak_memory_mb: f64,
}

impl AssemblyProfiler {
    pub fn new() -> Self {
        Self {
            start_time: Instant::now(),
            phases: Vec::new(),
            current_phase: None,
            peak_memory_mb: 0.0,
        }
    }

    /// Start timing a new phase
    pub fn start_phase(&mut self, name: impl Into<String>) {
        if let Some((prev_name, prev_start)) = self.current_phase.take() {
            // Auto-complete previous phase
            self.end_phase_internal(prev_name, prev_start, HashMap::new());
        }

        self.current_phase = Some((name.into(), Instant::now()));
    }

    /// End current phase with metadata
    pub fn end_phase(&mut self, metadata: HashMap<String, String>) {
        if let Some((name, start)) = self.current_phase.take() {
            self.end_phase_internal(name, start, metadata);
        }
    }

    fn end_phase_internal(
        &mut self,
        name: String,
        start: Instant,
        metadata: HashMap<String, String>,
    ) {
        let duration = start.elapsed();
        let memory = Self::get_current_memory_mb();

        if memory > self.peak_memory_mb {
            self.peak_memory_mb = memory;
        }

        let throughput = metadata
            .get("items_processed")
            .and_then(|s| s.parse::<f64>().ok())
            .map(|items| items / duration.as_secs_f64());

        let metrics = PhaseMetrics {
            phase_name: name,
            duration_ms: duration.as_millis() as u64,
            memory_mb: memory,
            throughput,
            metadata,
        };

        self.phases.push(metrics);
    }

    /// Generate final profiling report
    pub fn report(mut self, summary: ProfileSummary) -> AssemblyProfile {
        // Complete any running phase
        if let Some((name, start)) = self.current_phase.take() {
            self.end_phase_internal(name, start, HashMap::new());
        }

        AssemblyProfile {
            total_duration_ms: self.start_time.elapsed().as_millis() as u64,
            peak_memory_mb: self.peak_memory_mb,
            phases: self.phases,
            summary,
        }
    }

    /// Get current process memory usage in MB
    fn get_current_memory_mb() -> f64 {
        #[cfg(target_os = "linux")]
        {
            use std::fs;
            if let Ok(status) = fs::read_to_string("/proc/self/status") {
                for line in status.lines() {
                    if line.starts_with("VmRSS:") {
                        if let Some(kb) = line.split_whitespace().nth(1) {
                            if let Ok(kb_val) = kb.parse::<f64>() {
                                return kb_val / 1024.0; // Convert KB to MB
                            }
                        }
                    }
                }
            }
        }

        #[cfg(target_os = "macos")]
        {
            use std::process::Command;
            if let Ok(output) = Command::new("ps")
                .args(["-o", "rss=", "-p", &std::process::id().to_string()])
                .output()
            {
                if let Ok(rss_str) = String::from_utf8(output.stdout) {
                    if let Ok(kb) = rss_str.trim().parse::<f64>() {
                        return kb / 1024.0; // Convert KB to MB
                    }
                }
            }
        }

        0.0 // Fallback if we can't determine memory
    }
}

impl Default for AssemblyProfiler {
    fn default() -> Self {
        Self::new()
    }
}

/// Helper macro for automatic phase profiling
#[macro_export]
macro_rules! profile_phase {
    ($profiler:expr, $name:expr, $metadata:expr, $body:block) => {{
        $profiler.start_phase($name);
        let result = $body;
        $profiler.end_phase($metadata);
        result
    }};
}

/// Pretty-print profiling report
impl AssemblyProfile {
    pub fn print_report(&self) {
        use colored::Colorize;

        println!(
            "\n{}",
            "═══════════════════════════════════════════════".bright_cyan()
        );
        println!(
            "{}",
            "    ASSEMBLY PERFORMANCE PROFILE".bright_white().bold()
        );
        println!(
            "{}",
            "═══════════════════════════════════════════════".bright_cyan()
        );

        println!("\n{}", "📊 Summary:".bright_yellow().bold());
        println!(
            "  Total time:       {} ms",
            self.total_duration_ms.to_string().bright_white()
        );
        println!("  Peak memory:      {:.1} MB", self.peak_memory_mb);
        println!(
            "  Reads processed:  {}",
            self.summary.reads_processed.to_string().bright_white()
        );
        println!(
            "  K-mers counted:   {}",
            self.summary.kmers_counted.to_string().bright_white()
        );
        println!(
            "  Contigs created:  {}",
            self.summary.contigs_generated.to_string().bright_white()
        );
        println!(
            "  Throughput:       {:.0} reads/s",
            self.summary.reads_per_second
        );
        println!(
            "                    {:.0} k-mers/s",
            self.summary.kmers_per_second
        );

        println!("\n{}", "⏱️  Phase Breakdown:".bright_yellow().bold());
        println!(
            "  {:<30} {:>12} {:>12} {:>15}",
            "Phase", "Time (ms)", "Memory (MB)", "Throughput"
        );
        println!("  {}", "─".repeat(72).bright_black());

        for phase in &self.phases {
            let throughput_str = phase
                .throughput
                .map(|t| format!("{:.0} items/s", t))
                .unwrap_or_else(|| "-".to_string());

            println!(
                "  {:<30} {:>12} {:>12.1} {:>15}",
                phase.phase_name.bright_cyan(),
                phase.duration_ms,
                phase.memory_mb,
                throughput_str
            );

            // Print metadata if present
            if !phase.metadata.is_empty() {
                for (key, value) in &phase.metadata {
                    println!(
                        "    {} {}: {}",
                        "→".bright_blue(),
                        key.bright_black(),
                        value.bright_black()
                    );
                }
            }
        }

        println!(
            "\n{}",
            "═══════════════════════════════════════════════".bright_cyan()
        );

        // Performance warnings
        self.print_warnings();
    }

    fn print_warnings(&self) {
        use colored::Colorize;

        let mut warnings = Vec::new();

        // Check for slow k-mer counting
        if let Some(kmer_phase) = self.phases.iter().find(|p| p.phase_name.contains("K-mer")) {
            if let Some(throughput) = kmer_phase.throughput {
                if throughput < 10000.0 {
                    warnings.push(format!(
                        "⚠️  K-mer counting is slow ({:.0} k-mers/s). Expected >10,000/s",
                        throughput
                    ));
                }
            }
        }

        // Check for memory issues
        if self.peak_memory_mb > 8192.0 {
            warnings.push(format!(
                "⚠️  High memory usage ({:.1} MB). Consider reducing --memory-budget",
                self.peak_memory_mb
            ));
        }

        // Check phase imbalance
        if let (Some(longest), Some(shortest)) = (
            self.phases.iter().max_by_key(|p| p.duration_ms),
            self.phases.iter().min_by_key(|p| p.duration_ms),
        ) {
            let ratio = longest.duration_ms as f64 / shortest.duration_ms.max(1) as f64;
            if ratio > 50.0 {
                warnings.push(format!(
                    "⚠️  Severe phase imbalance: {} is {}x slower than {}",
                    longest.phase_name, ratio as u64, shortest.phase_name
                ));
            }
        }

        if !warnings.is_empty() {
            println!("\n{}", "⚠️  Performance Warnings:".bright_yellow().bold());
            for warning in warnings {
                println!("  {}", warning.yellow());
            }
        }
    }

    /// Export report as JSON
    pub fn to_json(&self) -> serde_json::Result<String> {
        serde_json::to_string_pretty(self)
    }

    /// Save report to file
    pub fn save_to_file(&self, path: impl AsRef<std::path::Path>) -> anyhow::Result<()> {
        use std::fs;
        let json = self.to_json()?;
        fs::write(path, json)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_profiler_basic() {
        let mut profiler = AssemblyProfiler::new();

        profiler.start_phase("test_phase");
        std::thread::sleep(Duration::from_millis(100));
        profiler.end_phase(HashMap::new());

        let summary = ProfileSummary {
            reads_processed: 1000,
            kmers_counted: 50000,
            nodes_created: 10000,
            edges_created: 15000,
            contigs_generated: 100,
            kmers_per_second: 500000.0,
            reads_per_second: 10000.0,
        };

        let report = profiler.report(summary);
        assert_eq!(report.phases.len(), 1);
        assert!(report.phases[0].duration_ms >= 100);
    }

    #[test]
    fn test_multiple_phases() {
        let mut profiler = AssemblyProfiler::new();

        for i in 0..3 {
            profiler.start_phase(format!("phase_{}", i));
            std::thread::sleep(Duration::from_millis(50));
            profiler.end_phase(HashMap::new());
        }

        let summary = ProfileSummary {
            reads_processed: 1000,
            kmers_counted: 50000,
            nodes_created: 10000,
            edges_created: 15000,
            contigs_generated: 100,
            kmers_per_second: 500000.0,
            reads_per_second: 10000.0,
        };

        let report = profiler.report(summary);
        assert_eq!(report.phases.len(), 3);
    }
}
