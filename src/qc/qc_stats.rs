//! QC Statistics and Reporting
//!
//! Tracks and reports QC metrics for preprocessing

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Comprehensive QC statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QCStats {
    pub reads_input: usize,
    pub reads_passed: usize,
    pub reads_failed: usize,

    // Quality filtering stats
    pub reads_failed_quality: usize,
    pub reads_failed_length: usize,
    pub reads_failed_adapter: usize,

    // Trimming stats
    pub bases_trimmed_quality: usize,
    pub bases_trimmed_adapter: usize,
    pub total_bases_before: usize,
    pub total_bases_after: usize,

    // Adapter stats
    pub adapters_detected: usize,
    pub adapter_types: HashMap<String, usize>,

    // Quality stats
    pub mean_quality_before: f64,
    pub mean_quality_after: f64,
    pub mean_length_before: f64,
    pub mean_length_after: f64,

    // Q20/Q30 stats
    pub q20_percentage_before: f64,
    pub q30_percentage_before: f64,
    pub q20_percentage_after: f64,
    pub q30_percentage_after: f64,
}

impl Default for QCStats {
    fn default() -> Self {
        Self {
            reads_input: 0,
            reads_passed: 0,
            reads_failed: 0,
            reads_failed_quality: 0,
            reads_failed_length: 0,
            reads_failed_adapter: 0,
            bases_trimmed_quality: 0,
            bases_trimmed_adapter: 0,
            total_bases_before: 0,
            total_bases_after: 0,
            adapters_detected: 0,
            adapter_types: HashMap::new(),
            mean_quality_before: 0.0,
            mean_quality_after: 0.0,
            mean_length_before: 0.0,
            mean_length_after: 0.0,
            q20_percentage_before: 0.0,
            q30_percentage_before: 0.0,
            q20_percentage_after: 0.0,
            q30_percentage_after: 0.0,
        }
    }
}

impl QCStats {
    pub fn new() -> Self {
        Self::default()
    }

    /// Record input read
    pub fn record_input(&mut self, sequence_length: usize, quality: &[u8]) {
        self.reads_input += 1;
        self.total_bases_before += sequence_length;

        // Update quality stats
        let qual_stats = Self::calculate_quality_stats(quality);
        self.mean_quality_before += qual_stats.mean;
        self.q20_percentage_before += qual_stats.q20_pct;
        self.q30_percentage_before += qual_stats.q30_pct;
    }

    /// Record passed read
    pub fn record_passed(&mut self, sequence_length: usize, quality: &[u8]) {
        self.reads_passed += 1;
        self.total_bases_after += sequence_length;

        // Update quality stats
        let qual_stats = Self::calculate_quality_stats(quality);
        self.mean_quality_after += qual_stats.mean;
        self.q20_percentage_after += qual_stats.q20_pct;
        self.q30_percentage_after += qual_stats.q30_pct;
    }

    /// Record failed read with reason
    pub fn record_failed(&mut self, reason: FailureReason) {
        self.reads_failed += 1;

        match reason {
            FailureReason::Quality => self.reads_failed_quality += 1,
            FailureReason::Length => self.reads_failed_length += 1,
            FailureReason::Adapter => self.reads_failed_adapter += 1,
        }
    }

    /// Record quality trimming
    pub fn record_quality_trimming(&mut self, bases_trimmed: usize) {
        self.bases_trimmed_quality += bases_trimmed;
    }

    /// Record adapter trimming
    pub fn record_adapter_trimming(&mut self, adapter_name: String, bases_trimmed: usize) {
        self.adapters_detected += 1;
        self.bases_trimmed_adapter += bases_trimmed;
        *self.adapter_types.entry(adapter_name).or_insert(0) += 1;
    }

    /// Calculate final statistics (averages)
    pub fn finalize(&mut self) {
        if self.reads_input > 0 {
            self.mean_quality_before /= self.reads_input as f64;
            self.mean_length_before = self.total_bases_before as f64 / self.reads_input as f64;
            self.q20_percentage_before /= self.reads_input as f64;
            self.q30_percentage_before /= self.reads_input as f64;
        }

        if self.reads_passed > 0 {
            self.mean_quality_after /= self.reads_passed as f64;
            self.mean_length_after = self.total_bases_after as f64 / self.reads_passed as f64;
            self.q20_percentage_after /= self.reads_passed as f64;
            self.q30_percentage_after /= self.reads_passed as f64;
        }
    }

    /// Calculate quality statistics for a single read
    fn calculate_quality_stats(quality: &[u8]) -> QualStats {
        if quality.is_empty() {
            return QualStats::default();
        }

        let qualities: Vec<u8> = quality.iter().map(|&q| q.saturating_sub(33)).collect();

        let mean = qualities.iter().map(|&q| q as f64).sum::<f64>() / qualities.len() as f64;
        let q20_count = qualities.iter().filter(|&&q| q >= 20).count();
        let q30_count = qualities.iter().filter(|&&q| q >= 30).count();

        QualStats {
            mean,
            q20_pct: (q20_count as f64 / qualities.len() as f64) * 100.0,
            q30_pct: (q30_count as f64 / qualities.len() as f64) * 100.0,
        }
    }

    /// Generate human-readable report
    pub fn generate_report(&self) -> QCReport {
        QCReport::new(self)
    }
}

#[derive(Debug, Clone, Default)]
struct QualStats {
    mean: f64,
    q20_pct: f64,
    q30_pct: f64,
}

/// Failure reasons
#[derive(Debug, Clone, Copy)]
pub enum FailureReason {
    Quality,
    Length,
    Adapter,
}

/// Human-readable QC report
#[derive(Debug, Clone)]
pub struct QCReport {
    pub summary: String,
    pub detailed: String,
    pub json: String,
}

impl QCReport {
    pub fn new(stats: &QCStats) -> Self {
        let summary = Self::generate_summary(stats);
        let detailed = Self::generate_detailed(stats);
        let json = serde_json::to_string_pretty(stats).unwrap_or_default();

        Self {
            summary,
            detailed,
            json,
        }
    }

    fn generate_summary(stats: &QCStats) -> String {
        use colored::Colorize;

        let pass_rate = if stats.reads_input > 0 {
            (stats.reads_passed as f64 / stats.reads_input as f64) * 100.0
        } else {
            0.0
        };

        format!(
            r#"
{}
    QC SUMMARY
{}

{}
  Input reads:        {}
  Passed reads:       {} ({:.1}%)
  Failed reads:       {} ({:.1}%)

{}
  Failed - Quality:   {}
  Failed - Length:    {}
  Failed - Adapter:   {}

{}
  Adapters detected:  {}
  Bases trimmed (Q):  {}
  Bases trimmed (A):  {}

{}
  Before: {:.1} bp avg, Q{:.1} avg
  After:  {:.1} bp avg, Q{:.1} avg
  Q20/Q30: {:.1}%/{:.1}% → {:.1}%/{:.1}%
"#,
            "═".repeat(50).bright_cyan(),
            "═".repeat(50).bright_cyan(),
            "📊 Read Statistics:".bright_yellow(),
            stats.reads_input.to_string().bright_white(),
            stats.reads_passed.to_string().bright_green(),
            pass_rate,
            stats.reads_failed.to_string().bright_red(),
            100.0 - pass_rate,
            "❌ Failure Breakdown:".bright_yellow(),
            stats.reads_failed_quality.to_string().bright_white(),
            stats.reads_failed_length.to_string().bright_white(),
            stats.reads_failed_adapter.to_string().bright_white(),
            "✂️  Trimming Statistics:".bright_yellow(),
            stats.adapters_detected.to_string().bright_white(),
            stats.bases_trimmed_quality.to_string().bright_white(),
            stats.bases_trimmed_adapter.to_string().bright_white(),
            "📈 Quality Improvement:".bright_yellow(),
            stats.mean_length_before,
            stats.mean_quality_before,
            stats.mean_length_after,
            stats.mean_quality_after,
            stats.q20_percentage_before,
            stats.q30_percentage_before,
            stats.q20_percentage_after,
            stats.q30_percentage_after,
        )
    }

    fn generate_detailed(stats: &QCStats) -> String {
        let mut report = String::new();

        report.push_str(&format!("Detailed QC Report\n"));
        report.push_str(&format!("==================\n\n"));

        report.push_str(&format!("Input Statistics:\n"));
        report.push_str(&format!("  Total reads: {}\n", stats.reads_input));
        report.push_str(&format!("  Total bases: {}\n", stats.total_bases_before));
        report.push_str(&format!(
            "  Mean length: {:.1} bp\n",
            stats.mean_length_before
        ));
        report.push_str(&format!(
            "  Mean quality: Q{:.1}\n\n",
            stats.mean_quality_before
        ));

        report.push_str(&format!("Filtering Results:\n"));
        report.push_str(&format!("  Passed: {} reads\n", stats.reads_passed));
        report.push_str(&format!("  Failed: {} reads\n", stats.reads_failed));
        report.push_str(&format!(
            "    - Quality: {}\n",
            stats.reads_failed_quality
        ));
        report.push_str(&format!("    - Length: {}\n", stats.reads_failed_length));
        report.push_str(&format!(
            "    - Adapter: {}\n\n",
            stats.reads_failed_adapter
        ));

        report.push_str(&format!("Adapter Detection:\n"));
        report.push_str(&format!("  Total adapters found: {}\n", stats.adapters_detected));
        for (adapter, count) in &stats.adapter_types {
            report.push_str(&format!("    {}: {} occurrences\n", adapter, count));
        }
        report.push_str("\n");

        report.push_str(&format!("Trimming Statistics:\n"));
        report.push_str(&format!(
            "  Quality-based trimming: {} bp\n",
            stats.bases_trimmed_quality
        ));
        report.push_str(&format!(
            "  Adapter-based trimming: {} bp\n",
            stats.bases_trimmed_adapter
        ));
        report.push_str(&format!(
            "  Total bases retained: {}\n\n",
            stats.total_bases_after
        ));

        report.push_str(&format!("Output Statistics:\n"));
        report.push_str(&format!("  Total reads: {}\n", stats.reads_passed));
        report.push_str(&format!("  Total bases: {}\n", stats.total_bases_after));
        report.push_str(&format!("  Mean length: {:.1} bp\n", stats.mean_length_after));
        report.push_str(&format!("  Mean quality: Q{:.1}\n", stats.mean_quality_after));
        report.push_str(&format!(
            "  Q20 percentage: {:.1}%\n",
            stats.q20_percentage_after
        ));
        report.push_str(&format!(
            "  Q30 percentage: {:.1}%\n",
            stats.q30_percentage_after
        ));

        report
    }

    /// Print summary to terminal
    pub fn print_summary(&self) {
        println!("{}", self.summary);
    }

    /// Print detailed report to terminal
    pub fn print_detailed(&self) {
        println!("{}", self.detailed);
    }

    /// Save report to file
    pub fn save_to_file(&self, path: &std::path::Path) -> anyhow::Result<()> {
        use std::fs;
        fs::write(path, &self.detailed)?;
        Ok(())
    }

    /// Save JSON to file
    pub fn save_json(&self, path: &std::path::Path) -> anyhow::Result<()> {
        use std::fs;
        fs::write(path, &self.json)?;
        Ok(())
    }
}
