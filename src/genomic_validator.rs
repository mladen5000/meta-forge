use std::collections::HashMap;
use anyhow::{Result, anyhow};
use serde::{Serialize, Deserialize};

/// Comprehensive validation tool for genomic data processing
pub struct GenomicDataValidator {
    /// Quality thresholds for different data types
    thresholds: ValidationThresholds,
    /// Validation statistics
    stats: ValidationStats,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct ValidationThresholds {
    pub min_sequence_length: usize,
    pub max_sequence_length: usize,
    pub min_quality_score: u8,
    pub max_n_content: f64,
    pub min_gc_content: f64,
    pub max_gc_content: f64,
    pub max_homopolymer_length: usize,
}

#[derive(Default, Serialize, Deserialize)]
pub struct ValidationStats {
    pub sequences_validated: usize,
    pub sequences_passed: usize,
    pub sequences_failed: usize,
    pub common_failures: HashMap<String, usize>,
}

#[derive(Debug, Clone)]
pub struct ValidationResult {
    pub passed: bool,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
    pub metrics: SequenceMetrics,
}

#[derive(Debug, Clone)]
pub struct SequenceMetrics {
    pub length: usize,
    pub gc_content: f64,
    pub n_content: f64,
    pub mean_quality: f64,
    pub max_homopolymer: usize,
    pub complexity_score: f64,
}

impl GenomicDataValidator {
    pub fn new() -> Self {
        Self {
            thresholds: ValidationThresholds::default(),
            stats: ValidationStats::default(),
        }
    }

    /// Validate a DNA sequence with quality scores
    pub fn validate_sequence(&mut self, sequence: &str, quality: Option<&[u8]>) -> ValidationResult {
        self.stats.sequences_validated += 1;
        
        let mut result = ValidationResult {
            passed: true,
            warnings: Vec::new(),
            errors: Vec::new(),
            metrics: self.calculate_metrics(sequence, quality),
        };

        // Length validation
        if sequence.len() < self.thresholds.min_sequence_length {
            result.errors.push(format!("Sequence too short: {} bp (min: {})", 
                sequence.len(), self.thresholds.min_sequence_length));
            result.passed = false;
        }

        if sequence.len() > self.thresholds.max_sequence_length {
            result.warnings.push(format!("Sequence very long: {} bp (max recommended: {})", 
                sequence.len(), self.thresholds.max_sequence_length));
        }

        // GC content validation
        if result.metrics.gc_content < self.thresholds.min_gc_content {
            result.warnings.push(format!("Low GC content: {:.2}% (min: {:.2}%)", 
                result.metrics.gc_content * 100.0, self.thresholds.min_gc_content * 100.0));
        }

        if result.metrics.gc_content > self.thresholds.max_gc_content {
            result.warnings.push(format!("High GC content: {:.2}% (max: {:.2}%)", 
                result.metrics.gc_content * 100.0, self.thresholds.max_gc_content * 100.0));
        }

        // N content validation
        if result.metrics.n_content > self.thresholds.max_n_content {
            result.errors.push(format!("Too many ambiguous bases: {:.2}% (max: {:.2}%)", 
                result.metrics.n_content * 100.0, self.thresholds.max_n_content * 100.0));
            result.passed = false;
        }

        // Quality validation
        if let Some(_) = quality {
            if result.metrics.mean_quality < self.thresholds.min_quality_score as f64 {
                result.warnings.push(format!("Low mean quality: {:.1} (min: {})", 
                    result.metrics.mean_quality, self.thresholds.min_quality_score));
            }
        }

        // Homopolymer validation
        if result.metrics.max_homopolymer > self.thresholds.max_homopolymer_length {
            result.warnings.push(format!("Long homopolymer run detected: {} bp", 
                result.metrics.max_homopolymer));
        }

        // Complexity validation
        if result.metrics.complexity_score < 0.5 {
            result.warnings.push(format!("Low sequence complexity: {:.3}", 
                result.metrics.complexity_score));
        }

        // Character validation
        for (i, ch) in sequence.chars().enumerate() {
            match ch.to_ascii_uppercase() {
                'A' | 'C' | 'G' | 'T' | 'N' => {},
                _ => {
                    result.errors.push(format!("Invalid character '{}' at position {}", ch, i));
                    result.passed = false;
                }
            }
        }

        // Update statistics
        if result.passed {
            self.stats.sequences_passed += 1;
        } else {
            self.stats.sequences_failed += 1;
            for error in &result.errors {
                *self.stats.common_failures.entry(error.clone()).or_insert(0) += 1;
            }
        }

        result
    }

    /// Validate k-mer size compatibility
    pub fn validate_kmer_size(&self, sequence_length: usize, k: usize) -> Result<()> {
        if k > sequence_length {
            return Err(anyhow!("K-mer size {} exceeds sequence length {}", k, sequence_length));
        }
        
        if k < 15 {
            return Err(anyhow!("K-mer size {} too small (minimum recommended: 15)", k));
        }
        
        if k > 255 {
            return Err(anyhow!("K-mer size {} too large (maximum supported: 255)", k));
        }
        
        if k % 2 == 0 {
            println!("Warning: Even k-mer size {} may cause assembly issues", k);
        }
        
        Ok(())
    }

    /// Validate file format and content
    pub fn validate_fastq_file(&mut self, file_path: &str) -> Result<ValidationSummary> {
        use bio::io::fastq;
        
        let mut summary = ValidationSummary::default();
        let reader = fastq::Reader::from_file(file_path)?;
        
        for (read_id, record_result) in reader.records().enumerate() {
            let record = record_result?;
            let sequence = std::str::from_utf8(record.seq())?;
            let quality = record.qual();
            
            let result = self.validate_sequence(sequence, Some(quality));
            summary.add_result(result);
            
            if read_id % 10000 == 0 && read_id > 0 {
                println!("Validated {} sequences...", read_id);
            }
        }
        
        Ok(summary)
    }

    fn calculate_metrics(&self, sequence: &str, quality: Option<&[u8]>) -> SequenceMetrics {
        let length = sequence.len();
        let mut gc_count = 0;
        let mut n_count = 0;
        let mut max_homopolymer = 0;
        let mut current_homopolymer = 1;
        let mut prev_char = '\0';

        // Count nucleotides and homopolymers
        for ch in sequence.chars() {
            match ch.to_ascii_uppercase() {
                'G' | 'C' => gc_count += 1,
                'N' => n_count += 1,
                _ => {}
            }

            if ch == prev_char {
                current_homopolymer += 1;
            } else {
                max_homopolymer = max_homopolymer.max(current_homopolymer);
                current_homopolymer = 1;
                prev_char = ch;
            }
        }
        max_homopolymer = max_homopolymer.max(current_homopolymer);

        // Calculate quality metrics
        let mean_quality = if let Some(qual) = quality {
            qual.iter().map(|&q| q as f64).sum::<f64>() / qual.len() as f64
        } else {
            0.0
        };

        // Calculate complexity using Shannon entropy
        let complexity_score = self.calculate_shannon_entropy(sequence);

        SequenceMetrics {
            length,
            gc_content: gc_count as f64 / length as f64,
            n_content: n_count as f64 / length as f64,
            mean_quality,
            max_homopolymer,
            complexity_score,
        }
    }

    fn calculate_shannon_entropy(&self, sequence: &str) -> f64 {
        let mut counts = [0u32; 4]; // A, C, G, T
        let mut total = 0u32;

        for ch in sequence.chars() {
            match ch.to_ascii_uppercase() {
                'A' => { counts[0] += 1; total += 1; },
                'C' => { counts[1] += 1; total += 1; },
                'G' => { counts[2] += 1; total += 1; },
                'T' => { counts[3] += 1; total += 1; },
                _ => {} // Skip N's and other characters
            }
        }

        if total == 0 { return 0.0; }

        let entropy = counts.iter()
            .filter(|&&c| c > 0)
            .map(|&c| {
                let p = c as f64 / total as f64;
                -p * p.log2()
            })
            .sum::<f64>();

        entropy / 2.0 // Normalize to [0, 1]
    }

    pub fn get_stats(&self) -> &ValidationStats {
        &self.stats
    }

    pub fn print_report(&self) {
        println!("\n=== Genomic Data Validation Report ===");
        println!("Total sequences validated: {}", self.stats.sequences_validated);
        println!("Passed validation: {}", self.stats.sequences_passed);
        println!("Failed validation: {}", self.stats.sequences_failed);
        println!("Success rate: {:.2}%", 
            (self.stats.sequences_passed as f64 / self.stats.sequences_validated as f64) * 100.0);
        
        if !self.stats.common_failures.is_empty() {
            println!("\nCommon failure types:");
            let mut failures: Vec<_> = self.stats.common_failures.iter().collect();
            failures.sort_by(|a, b| b.1.cmp(a.1));
            
            for (error, count) in failures.iter().take(5) {
                println!("  {}: {} occurrences", error, count);
            }
        }
    }
}

impl Default for ValidationThresholds {
    fn default() -> Self {
        Self {
            min_sequence_length: 50,
            max_sequence_length: 1_000_000,
            min_quality_score: 20,
            max_n_content: 0.1, // 10%
            min_gc_content: 0.2, // 20%
            max_gc_content: 0.8, // 80%
            max_homopolymer_length: 20,
        }
    }
}

#[derive(Default)]
pub struct ValidationSummary {
    pub total_sequences: usize,
    pub passed_sequences: usize,
    pub total_warnings: usize,
    pub total_errors: usize,
    pub avg_length: f64,
    pub avg_gc_content: f64,
    pub avg_quality: f64,
}

impl ValidationSummary {
    fn add_result(&mut self, result: ValidationResult) {
        self.total_sequences += 1;
        if result.passed {
            self.passed_sequences += 1;
        }
        self.total_warnings += result.warnings.len();
        self.total_errors += result.errors.len();
        
        // Update running averages
        let n = self.total_sequences as f64;
        self.avg_length = (self.avg_length * (n - 1.0) + result.metrics.length as f64) / n;
        self.avg_gc_content = (self.avg_gc_content * (n - 1.0) + result.metrics.gc_content) / n;
        self.avg_quality = (self.avg_quality * (n - 1.0) + result.metrics.mean_quality) / n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_validation() {
        let mut validator = GenomicDataValidator::new();
        
        // Valid sequence
        let result = validator.validate_sequence("ATCGATCGATCGATCG", None);
        assert!(result.passed);
        assert!(result.errors.is_empty());
        
        // Invalid character
        let result = validator.validate_sequence("ATCGATCXATCGATCG", None);
        assert!(!result.passed);
        assert!(!result.errors.is_empty());
        
        // Too short
        let result = validator.validate_sequence("ATCG", None);
        assert!(!result.passed);
    }

    #[test]
    fn test_metrics_calculation() {
        let validator = GenomicDataValidator::new();
        let metrics = validator.calculate_metrics("AAATTTCCCGGG", None);
        
        assert_eq!(metrics.length, 12);
        assert!((metrics.gc_content - 0.5).abs() < 0.01);
        assert!(metrics.complexity_score > 0.8); // Should be high complexity
    }
}