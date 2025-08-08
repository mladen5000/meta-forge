//! GenomicDataValidator
//!
//! A self-contained quality-control toolkit for FASTQ / raw-string
//! validation.  Designed for integration in CLI or pipeline code.

use std::collections::HashMap;

use anyhow::{Result, bail};
use bio::io::fastq;
use serde::{Deserialize, Serialize};

/* ------------------------------------------------------------------------- */
/*                               CORE STRUCTS                                */
/* ------------------------------------------------------------------------- */

/// High-level validator holding tunable thresholds plus run statistics.
pub struct GenomicDataValidator {
    thresholds: ValidationThresholds,
    stats: ValidationStats,
}

/// Thresholds controlling pass/fail decisions.
///
/// The defaults are lenient enough for short-read Illumina data.
#[derive(Clone, Serialize, Deserialize)]
pub struct ValidationThresholds {
    pub min_sequence_length: usize,
    pub max_sequence_length: usize,
    pub min_quality_score: u8,
    pub max_n_content: f64, // proportion 0–1
    pub min_gc_content: f64,
    pub max_gc_content: f64,
    pub max_homopolymer_len: usize,
}

/// Cumulative run statistics.
#[derive(Default, Serialize, Deserialize)]
pub struct ValidationStats {
    pub sequences_validated: usize,
    pub sequences_passed: usize,
    pub sequences_failed: usize,
    pub common_failures: HashMap<String, usize>,
}

/// Per-sequence QC result.
#[derive(Debug, Clone)]
pub struct ValidationResult {
    pub passed: bool,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
    pub metrics: SequenceMetrics,
}

/// Per-sequence computed metrics.
#[derive(Debug, Clone)]
pub struct SequenceMetrics {
    pub length: usize,
    pub gc_content: f64,
    pub n_content: f64,
    pub mean_quality: f64,
    pub max_homopolymer: usize,
    pub complexity: f64,
}

/* ------------------------------------------------------------------------- */
/*                               IMPLEMENTATION                              */
/* ------------------------------------------------------------------------- */

impl Default for GenomicDataValidator {
    fn default() -> Self {
        Self::new()
    }
}

impl GenomicDataValidator {
    /// Create a new validator with default thresholds.
    pub fn new() -> Self {
        Self {
            thresholds: ValidationThresholds::default(),
            stats: ValidationStats::default(),
        }
    }

    /// Main entry point: validate a single sequence + optional qualities.
    pub fn validate_sequence(
        &mut self,
        sequence: &str,
        quality: Option<&[u8]>,
    ) -> ValidationResult {
        self.stats.sequences_validated += 1;

        let metrics = self.calculate_metrics(sequence, quality);
        let mut res = ValidationResult {
            passed: true,
            warnings: Vec::new(),
            errors: Vec::new(),
            metrics: metrics.clone(),
        };

        // ---- Hard failures -------------------------------------------------
        if metrics.length < self.thresholds.min_sequence_length {
            res.fail(format!(
                "Sequence too short: {} bp (min {})",
                metrics.length, self.thresholds.min_sequence_length
            ));
        }
        if metrics.n_content > self.thresholds.max_n_content {
            res.fail(format!(
                "N content {:.2}% > {:.2}%",
                metrics.n_content * 100.0,
                self.thresholds.max_n_content * 100.0
            ));
        }
        for (i, ch) in sequence.chars().enumerate() {
            if !matches!(ch.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'N') {
                res.fail(format!("Invalid base '{ch}' at position {i}"));
            }
        }

        // ---- Soft warnings -------------------------------------------------
        if metrics.length > self.thresholds.max_sequence_length {
            res.warn(format!(
                "Unusually long: {} bp (recommended max {})",
                metrics.length, self.thresholds.max_sequence_length
            ));
        }
        if metrics.gc_content < self.thresholds.min_gc_content {
            res.warn(format!(
                "GC {:.2}% < {:.2}%",
                metrics.gc_content * 100.0,
                self.thresholds.min_gc_content * 100.0
            ));
        }
        if metrics.gc_content > self.thresholds.max_gc_content {
            res.warn(format!(
                "GC {:.2}% > {:.2}%",
                metrics.gc_content * 100.0,
                self.thresholds.max_gc_content * 100.0
            ));
        }
        if metrics.mean_quality < self.thresholds.min_quality_score as f64 {
            res.warn(format!(
                "Mean Q {:.1} < {}",
                metrics.mean_quality, self.thresholds.min_quality_score
            ));
        }
        if metrics.max_homopolymer > self.thresholds.max_homopolymer_len {
            res.warn(format!(
                "Long homopolymer: {} bp (max {})",
                metrics.max_homopolymer, self.thresholds.max_homopolymer_len
            ));
        }
        if metrics.complexity < 0.5 {
            res.warn(format!("Low complexity: {:.3}", metrics.complexity));
        }

        // ---- Stats aggregation --------------------------------------------
        if res.passed {
            self.stats.sequences_passed += 1;
        } else {
            self.stats.sequences_failed += 1;
            for e in &res.errors {
                *self.stats.common_failures.entry(e.clone()).or_insert(0) += 1;
            }
        }
        res
    }

    /// Fast sanity check for k-mer sizes.
    pub fn validate_kmer_size(&self, seq_len: usize, k: usize) -> Result<()> {
        if k > seq_len {
            bail!("k={} exceeds sequence length {}", k, seq_len);
        }
        if !(15..=255).contains(&k) {
            bail!("k={} outside supported 15–255 range", k);
        }
        if k % 2 == 0 {
            eprintln!("⚠️  Even k={k} can hinder assembly");
        }
        Ok(())
    }

    /// Iterate through a FASTQ and return a roll-up summary.
    pub fn validate_fastq_file(&mut self, path: &str) -> Result<ValidationSummary> {
        let mut summary = ValidationSummary::default();
        let reader = fastq::Reader::from_file(path)?;

        for (idx, rec) in reader.records().enumerate() {
            let rec = rec?;
            let seq = std::str::from_utf8(rec.seq())?;
            let res = self.validate_sequence(seq, Some(rec.qual()));
            summary.add(res);

            if idx % 10_000 == 9 {
                eprintln!("  validated {} reads …", idx + 1);
            }
        }
        Ok(summary)
    }

    /* --------------------- internal helper functions --------------------- */

    fn calculate_metrics(&self, seq: &str, qual: Option<&[u8]>) -> SequenceMetrics {
        let len = seq.len();
        if len == 0 {
            return SequenceMetrics {
                length: 0,
                gc_content: 0.0,
                n_content: 0.0,
                mean_quality: 0.0,
                max_homopolymer: 0,
                complexity: 0.0,
            };
        }

        let mut gc = 0;
        let mut n = 0;
        let mut max_hpoly = 0;
        let mut cur = 0;
        let mut prev = '\0';

        for ch in seq.chars() {
            let upper = ch.to_ascii_uppercase();
            match upper {
                'G' | 'C' => gc += 1,
                'N' => n += 1,
                _ => {}
            };

            if upper == prev {
                cur += 1;
            } else {
                max_hpoly = max_hpoly.max(cur);
                cur = 1;
                prev = upper;
            }
        }
        max_hpoly = max_hpoly.max(cur);

        let mean_q = qual
            .map(|q| q.iter().map(|&b| b as f64).sum::<f64>() / q.len() as f64)
            .unwrap_or(0.0);

        SequenceMetrics {
            length: len,
            gc_content: gc as f64 / len as f64,
            n_content: n as f64 / len as f64,
            mean_quality: mean_q,
            max_homopolymer: max_hpoly,
            complexity: Self::shannon_entropy(seq),
        }
    }

    /// Shannon entropy normalised to [0, 1].
    fn shannon_entropy(seq: &str) -> f64 {
        let mut counts = [0u32; 4]; // A,C,G,T
        let mut total = 0u32;
        for ch in seq.chars() {
            match ch.to_ascii_uppercase() {
                'A' => counts[0] += 1,
                'C' => counts[1] += 1,
                'G' => counts[2] += 1,
                'T' => counts[3] += 1,
                _ => {}
            }
            total += 1;
        }
        if total == 0 {
            return 0.0;
        }
        let entropy = counts
            .iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / total as f64;
                -p * p.log2()
            })
            .sum::<f64>();
        entropy / 2.0 // log2(4)=2
    }

    /// Public accessor for cumulative stats.
    pub fn stats(&self) -> &ValidationStats {
        &self.stats
    }

    /// Pretty console report.
    pub fn print_report(&self) {
        println!(
            "\nValidated: {}, Passed: {}, Failed: {}  ({:.2} % pass)",
            self.stats.sequences_validated,
            self.stats.sequences_passed,
            self.stats.sequences_failed,
            100.0 * self.stats.sequences_passed as f64
                / self.stats.sequences_validated.max(1) as f64
        );
        if !self.stats.common_failures.is_empty() {
            println!("Top failure causes:");
            let mut v: Vec<_> = self.stats.common_failures.iter().collect();
            v.sort_by(|a, b| b.1.cmp(a.1));
            for (err, n) in v.into_iter().take(5) {
                println!("  {n} – {err}");
            }
        }
    }
}

/* ------------------------------------------------------------------------- */
/*                             SUPPORT STRUCTS                               */
/* ------------------------------------------------------------------------- */

impl ValidationResult {
    fn warn(&mut self, m: String) {
        self.warnings.push(m);
    }
    fn fail(&mut self, m: String) {
        self.errors.push(m);
        self.passed = false;
    }
}

impl Default for ValidationThresholds {
    fn default() -> Self {
        Self {
            min_sequence_length: 50,
            max_sequence_length: 1_000_000,
            min_quality_score: 20,
            max_n_content: 0.10,
            min_gc_content: 0.20,
            max_gc_content: 0.80,
            max_homopolymer_len: 20,
        }
    }
}

/// Run-level roll-up.
#[derive(Default)]
pub struct ValidationSummary {
    pub total: usize,
    pub passed: usize,
    pub warnings: usize,
    pub errors: usize,
    pub avg_len: f64,
    pub avg_gc: f64,
    pub avg_quality: f64,
}

impl ValidationSummary {
    fn add(&mut self, r: ValidationResult) {
        self.total += 1;
        self.passed += usize::from(r.passed);
        self.warnings += r.warnings.len();
        self.errors += r.errors.len();

        let n = self.total as f64;
        self.avg_len = (self.avg_len * (n - 1.0) + r.metrics.length as f64) / n;
        self.avg_gc = (self.avg_gc * (n - 1.0) + r.metrics.gc_content) / n;
        self.avg_quality = (self.avg_quality * (n - 1.0) + r.metrics.mean_quality) / n;
    }
}

/* ------------------------------------------------------------------------- */
/*                                    TESTS                                  */
/* ------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_pass() {
        let mut v = GenomicDataValidator::new();
        let res = v.validate_sequence("ATCGATCGATCG", None);
        assert!(res.passed);
        assert!(res.errors.is_empty());
    }

    #[test]
    fn catches_bad_base_and_length() {
        let mut v = GenomicDataValidator::new();
        let r = v.validate_sequence("AXT", None);
        assert!(!r.passed);
        assert!(r.errors.iter().any(|e| e.contains("Invalid base")));
        assert!(r.errors.iter().any(|e| e.contains("too short")));
    }

    #[test]
    fn entropy_boundaries() {
        let v = GenomicDataValidator::new();
        let complex = v.calculate_metrics("ACGTACGT", None).complexity;
        let low = v.calculate_metrics("AAAAAAAA", None).complexity;
        assert!(complex > low);
        assert!((0.0..=1.0).contains(&complex));
    }

    use proptest::proptest;
    
    proptest! {
        #[test]
        fn fuzz_no_panic(seq in "[ACGTN]{0,200}") {
            let mut v = GenomicDataValidator::new();
            let _ = v.validate_sequence(&seq, None);
        }
    }
}
