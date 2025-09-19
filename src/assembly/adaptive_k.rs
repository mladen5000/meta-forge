//! Adaptive K-mer Selection for Laptop Assembly
//! ============================================
//!
//! Simple adaptive k-mer size selection based on read characteristics.
//! Optimized for typical laptop memory constraints and processing power.

use anyhow::Result;
use crate::core::data_structures::CorrectedRead;

/// Configuration for adaptive k-mer selection
#[derive(Debug, Clone)]
pub struct AdaptiveKConfig {
    /// Minimum k-mer size to consider
    pub min_k: usize,
    /// Maximum k-mer size to consider
    pub max_k: usize,
    /// Sample size for read analysis
    pub sample_size: usize,
    /// Memory budget in MB for k-mer selection
    pub memory_budget_mb: usize,
}

impl Default for AdaptiveKConfig {
    fn default() -> Self {
        Self {
            min_k: 15,
            max_k: 31,
            sample_size: 1000,
            memory_budget_mb: 512,
        }
    }
}

/// Adaptive k-mer selector
pub struct AdaptiveKSelector {
    config: AdaptiveKConfig,
}

impl AdaptiveKSelector {
    /// Create new adaptive k-mer selector
    pub fn new(config: AdaptiveKConfig) -> Self {
        Self { config }
    }

    /// Select optimal k-mer size based on read characteristics
    pub fn select_optimal_k(&self, reads: &[CorrectedRead]) -> Result<usize> {
        if reads.is_empty() {
            return Ok(21); // Reasonable default
        }

        // Sample reads for analysis
        let sample_size = self.config.sample_size.min(reads.len());
        let sample: Vec<_> = reads.iter().take(sample_size).collect();

        // Analyze read characteristics
        let stats = self.analyze_read_characteristics(&sample);

        // Select k based on characteristics and memory constraints
        let optimal_k = self.calculate_optimal_k(&stats);

        // Ensure k is within configured bounds
        let bounded_k = optimal_k
            .max(self.config.min_k)
            .min(self.config.max_k);

        println!("📏 Selected k-mer size: {} (avg_len={:.1}, gc={:.1}%)",
                bounded_k, stats.avg_length, stats.gc_content);

        Ok(bounded_k)
    }

    /// Analyze characteristics of read sample
    fn analyze_read_characteristics(&self, sample: &[&CorrectedRead]) -> ReadCharacteristics {
        let total_length: usize = sample.iter().map(|r| r.corrected.len()).sum();
        let avg_length = total_length as f64 / sample.len() as f64;

        // Calculate GC content
        let mut gc_count = 0;
        let mut total_bases = 0;

        for read in sample {
            for c in read.corrected.chars() {
                match c.to_ascii_uppercase() {
                    'G' | 'C' => {
                        gc_count += 1;
                        total_bases += 1;
                    }
                    'A' | 'T' => {
                        total_bases += 1;
                    }
                    _ => {} // Skip ambiguous bases
                }
            }
        }

        let gc_content = if total_bases > 0 {
            (gc_count as f64 / total_bases as f64) * 100.0
        } else {
            50.0 // Default assumption
        };

        // Estimate complexity by checking repeat patterns
        let complexity = self.estimate_sequence_complexity(sample);

        ReadCharacteristics {
            avg_length,
            gc_content,
            complexity,
            total_reads: sample.len(),
        }
    }

    /// Estimate sequence complexity (simple approach)
    fn estimate_sequence_complexity(&self, sample: &[&CorrectedRead]) -> f64 {
        // Simple approach: look at 4-mer diversity in a subsample
        let mut kmer_4_counts = std::collections::HashMap::new();
        let mut total_4mers = 0;

        for read in sample.iter().take(100) { // Limit for performance
            let seq = &read.corrected;
            if seq.len() >= 4 {
                for i in 0..=seq.len() - 4 {
                    let kmer = &seq[i..i + 4];
                    *kmer_4_counts.entry(kmer.to_string()).or_insert(0) += 1;
                    total_4mers += 1;
                }
            }
        }

        if total_4mers == 0 {
            return 1.0; // Default complexity
        }

        // Calculate Shannon entropy
        let mut entropy = 0.0;
        for count in kmer_4_counts.values() {
            let p = *count as f64 / total_4mers as f64;
            entropy -= p * p.log2();
        }

        // Normalize to 0-1 range (4-mer entropy ranges from 0 to ~8)
        (entropy / 8.0).min(1.0).max(0.0)
    }

    /// Calculate optimal k-mer size based on characteristics
    fn calculate_optimal_k(&self, stats: &ReadCharacteristics) -> usize {
        // Base k-mer size on average read length
        let length_based_k: usize = if stats.avg_length < 50.0 {
            15 // Short reads
        } else if stats.avg_length < 100.0 {
            21 // Medium reads
        } else if stats.avg_length < 200.0 {
            31 // Long reads
        } else {
            41 // Very long reads
        };

        // Adjust for GC content (extreme GC can cause problems)
        let gc_adjusted_k = if stats.gc_content < 30.0 || stats.gc_content > 70.0 {
            // Extreme GC content - use smaller k for better sensitivity
            length_based_k.saturating_sub(4)
        } else {
            length_based_k
        };

        // Adjust for sequence complexity
        let complexity_adjusted_k = if stats.complexity < 0.5 {
            // Low complexity sequences need smaller k to avoid gaps
            gc_adjusted_k.saturating_sub(6)
        } else if stats.complexity > 0.8 {
            // High complexity can use larger k for specificity
            gc_adjusted_k.saturating_add(4)
        } else {
            gc_adjusted_k
        };

        // Memory constraint check
        let memory_limited_k = self.apply_memory_constraints(complexity_adjusted_k);

        memory_limited_k
    }

    /// Apply memory constraints to k-mer selection
    fn apply_memory_constraints(&self, candidate_k: usize) -> usize {
        // Rough estimate: 4^k possible k-mers, each taking ~12 bytes
        // Use safe arithmetic to prevent overflow
        let base_4_power = 4_u64.saturating_pow(candidate_k as u32);
        let estimated_memory_mb = (base_4_power.saturating_mul(12)) as f64 / (1024.0 * 1024.0);

        // If estimated memory exceeds budget, reduce k
        if estimated_memory_mb > self.config.memory_budget_mb as f64 {
            // Find largest k that fits in memory budget
            for k in (self.config.min_k..=candidate_k).rev() {
                let base_4_power = 4_u64.saturating_pow(k as u32);
                let mem_mb = (base_4_power.saturating_mul(12)) as f64 / (1024.0 * 1024.0);
                if mem_mb <= self.config.memory_budget_mb as f64 {
                    return k;
                }
            }
            self.config.min_k // Fallback to minimum
        } else {
            candidate_k
        }
    }
}

/// Characteristics of read sample
#[derive(Debug)]
struct ReadCharacteristics {
    avg_length: f64,
    gc_content: f64,
    complexity: f64,
    total_reads: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::data_structures::CorrectionMetadata;

    fn create_test_read(sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id: 0,
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

    #[test]
    fn test_adaptive_k_selection() {
        let config = AdaptiveKConfig::default();
        let selector = AdaptiveKSelector::new(config);

        // Test with short reads
        let short_reads = vec![
            create_test_read("ATCGATCGATCG"),
            create_test_read("GCTAGCTAGCTA"),
            create_test_read("TTAACCGGAATT"),
        ];

        let k = selector.select_optimal_k(&short_reads).unwrap();
        assert!(k >= 15 && k <= 31);

        // Test with long reads
        let long_reads = vec![
            create_test_read(&"ATCGATCGATCG".repeat(10)),
            create_test_read(&"GCTAGCTAGCTA".repeat(10)),
        ];

        let k_long = selector.select_optimal_k(&long_reads).unwrap();
        assert!(k_long >= k); // Longer reads should allow larger k
    }

    #[test]
    fn test_gc_content_adjustment() {
        let config = AdaptiveKConfig::default();
        let selector = AdaptiveKSelector::new(config);

        // High GC content reads
        let high_gc_reads = vec![
            create_test_read(&"GCGCGCGCGCGC".repeat(5)),
            create_test_read(&"CCGGCCGGCCGG".repeat(5)),
        ];

        let k_high_gc = selector.select_optimal_k(&high_gc_reads).unwrap();

        // Normal GC content reads
        let normal_reads = vec![
            create_test_read(&"ATCGATCGATCG".repeat(5)),
            create_test_read(&"TACGTACGTACG".repeat(5)),
        ];

        let k_normal = selector.select_optimal_k(&normal_reads).unwrap();

        // High GC should tend to use smaller k
        assert!(k_high_gc <= k_normal + 2);
    }

    #[test]
    fn test_memory_constraints() {
        let config = AdaptiveKConfig {
            min_k: 15,
            max_k: 25, // Reduced max to test constraint
            memory_budget_mb: 1, // Very low budget
            ..Default::default()
        };
        let selector = AdaptiveKSelector::new(config);

        let reads = vec![
            create_test_read(&"ATCGATCGATCG".repeat(10)),
        ];

        let k = selector.select_optimal_k(&reads).unwrap();
        assert!(k <= 25); // Should respect max_k constraint
    }
}