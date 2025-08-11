use anyhow::Result;
use rayon::prelude::*;
use std::collections::HashMap;

#[cfg(test)]
use fastrand;

/// Hybrid streaming sketch combining HyperLogLog (F₀) with L₀ sampling
#[derive(Clone)]
pub struct HybridAbundanceEstimator {
    /// HyperLogLog for cardinality estimation
    hyperloglog: HyperLogLog,
    /// L₀ sampler for representative k-mers
    l0_sampler: L0Sampler,
    /// Combined sketch parameters
    config: SketchConfig,
    /// Statistics tracking
    stats: SketchStats,
}

/// HyperLogLog implementation for F₀ (unique count) estimation
#[derive(Clone)]
pub struct HyperLogLog {
    /// Number of buckets (2^precision)
    num_buckets: usize,
    /// Precision parameter (typically 10-16)
    precision: u8,
    /// Bucket registers storing max leading zeros
    buckets: Vec<u8>,
    /// Alpha constant for bias correction
    alpha: f64,
}

/// L₀ sampler for retaining representative elements
#[derive(Clone)]
pub struct L0Sampler {
    /// Current sample set with priorities
    samples: HashMap<u64, f64>,
    /// Maximum sample size
    max_samples: usize,
    /// Threshold for inclusion (updated dynamically)
    threshold: f64,
    /// Random number generator seed
    rng_state: u64,
}

#[derive(Clone)]
pub struct SketchConfig {
    /// HyperLogLog precision (10-16)
    hll_precision: u8,
    /// Maximum L₀ samples to retain
    max_l0_samples: usize,
    /// K-mer size for sketching
    k: usize,
    /// Enable threaded processing
    parallel_processing: bool,
    /// Memory limit in MB
    memory_limit_mb: usize,
}

#[derive(Default, Clone)]
pub struct SketchStats {
    /// Total k-mers processed
    total_kmers: u64,
    /// Estimated unique k-mers (F₀)
    estimated_unique: u64,
    /// Actual L₀ samples retained
    l0_samples_count: usize,
    /// Processing time statistics
    processing_time_ms: u64,
    /// Memory usage estimate
    memory_usage_bytes: usize,
}

/// Result of sketch merging operation
pub struct MergeResult {
    pub merged_sketch: HybridAbundanceEstimator,
    pub merge_stats: MergeStats,
}

#[derive(Default)]
pub struct MergeStats {
    pub sketches_merged: usize,
    pub total_samples_before: usize,
    pub total_samples_after: usize,
    pub estimated_unique_before: Vec<u64>,
    pub estimated_unique_after: u64,
}

impl HybridAbundanceEstimator {
    pub fn new(config: SketchConfig) -> Self {
        let hyperloglog = HyperLogLog::new(config.hll_precision);
        let l0_sampler = L0Sampler::new(config.max_l0_samples);

        Self {
            hyperloglog,
            l0_sampler,
            config,
            stats: SketchStats::default(),
        }
    }

    /// Process a single k-mer and update both sketches
    pub fn add_kmer(&mut self, kmer_hash: u64) {
        // Update HyperLogLog
        self.hyperloglog.add(kmer_hash);

        // Update L₀ sampler
        self.l0_sampler.add(kmer_hash);

        // Update statistics
        self.stats.total_kmers += 1;
        self.update_memory_estimate();
    }

    /// Process a sequence and extract k-mers
    pub fn add_sequence(&mut self, sequence: &str) -> Result<()> {
        let seq_bytes = sequence.as_bytes();

        if seq_bytes.len() < self.config.k {
            return Ok(());
        }

        if self.config.parallel_processing && sequence.len() > 10000 {
            self.add_sequence_parallel(sequence)
        } else {
            self.add_sequence_sequential(sequence)
        }
    }

    fn add_sequence_sequential(&mut self, sequence: &str) -> Result<()> {
        let seq_bytes = sequence.as_bytes();

        for window in seq_bytes.windows(self.config.k) {
            let kmer_hash = self.hash_kmer(window);
            self.add_kmer(kmer_hash);
        }

        Ok(())
    }

    fn add_sequence_parallel(&mut self, sequence: &str) -> Result<()> {
        let seq_bytes = sequence.as_bytes();
        let chunk_size = 1000; // Process in chunks for parallelization

        // Collect k-mer hashes in parallel
        let hashes: Vec<u64> = seq_bytes
            .windows(self.config.k)
            .collect::<Vec<_>>()
            .par_chunks(chunk_size)
            .map(|chunk| {
                chunk
                    .iter()
                    .map(|window| self.hash_kmer(window))
                    .collect::<Vec<_>>()
            })
            .flatten()
            .collect();

        // Add all hashes to sketches
        for hash in hashes {
            self.add_kmer(hash);
        }

        Ok(())
    }

    /// Get current F₀ estimate (unique k-mer count)
    pub fn estimate_unique_count(&mut self) -> u64 {
        let estimate = self.hyperloglog.estimate();
        self.stats.estimated_unique = estimate;
        estimate
    }

    /// Get current L₀ samples with abundance estimates
    pub fn get_representative_samples(&self) -> HashMap<u64, f64> {
        self.l0_sampler.get_samples_with_weights()
    }

    /// Get abundance estimate for a specific k-mer
    pub fn estimate_abundance(&self, kmer_hash: u64) -> Option<f64> {
        self.l0_sampler.get_sample_weight(kmer_hash)
    }

    /// Merge multiple sketches from parallel threads
    pub fn merge_sketches(sketches: Vec<Self>) -> Result<MergeResult> {
        if sketches.is_empty() {
            return Err(anyhow::anyhow!("Cannot merge empty sketch list"));
        }

        let mut merge_stats = MergeStats {
            sketches_merged: sketches.len(),
            total_samples_before: sketches.iter().map(|s| s.l0_sampler.samples.len()).sum(),
            estimated_unique_before: sketches.iter().map(|s| s.stats.estimated_unique).collect(),
            ..Default::default()
        };

        // Use first sketch as base
        let mut merged = sketches[0].clone();

        // Merge remaining sketches
        for sketch in sketches.into_iter().skip(1) {
            merged.merge_with(sketch)?;
        }

        // Update final statistics
        merge_stats.total_samples_after = merged.l0_sampler.samples.len();
        merge_stats.estimated_unique_after = merged.estimate_unique_count();

        Ok(MergeResult {
            merged_sketch: merged,
            merge_stats,
        })
    }

    fn merge_with(&mut self, other: Self) -> Result<()> {
        // Merge HyperLogLog sketches
        self.hyperloglog.merge(&other.hyperloglog)?;

        // Merge L₀ samplers
        self.l0_sampler.merge(&other.l0_sampler);

        // Combine statistics
        self.stats.total_kmers += other.stats.total_kmers;
        self.stats.processing_time_ms += other.stats.processing_time_ms;

        self.update_memory_estimate();
        Ok(())
    }

    /// Export sketch for serialization/storage
    pub fn export_sketch(&self) -> SketchExport {
        SketchExport {
            hll_buckets: self.hyperloglog.buckets.clone(),
            hll_precision: self.hyperloglog.precision,
            l0_samples: self.l0_sampler.samples.clone(),
            l0_threshold: self.l0_sampler.threshold,
            config: self.config.clone(),
            stats: self.stats.clone(),
        }
    }

    /// Import sketch from serialized data
    pub fn import_sketch(export: SketchExport) -> Self {
        let hyperloglog = HyperLogLog::from_buckets(export.hll_buckets, export.hll_precision);
        let l0_sampler = L0Sampler::from_samples(
            export.l0_samples,
            export.l0_threshold,
            export.config.max_l0_samples,
        );

        Self {
            hyperloglog,
            l0_sampler,
            config: export.config,
            stats: export.stats,
        }
    }

    fn hash_kmer(&self, kmer: &[u8]) -> u64 {
        // Use xxHash or similar high-quality hash function
        // For demo, using FNV-1a
        const FNV_OFFSET: u64 = 0xcbf29ce484222325;
        const FNV_PRIME: u64 = 0x100000001b3;
        kmer.iter().fold(FNV_OFFSET, |acc, &b| {
            (acc ^ b as u64).wrapping_mul(FNV_PRIME)
        })
    }

    fn update_memory_estimate(&mut self) {
        self.stats.memory_usage_bytes = self.hyperloglog.memory_usage()
            + self.l0_sampler.memory_usage()
            + std::mem::size_of::<SketchStats>();

        self.stats.l0_samples_count = self.l0_sampler.samples.len();
    }

    pub fn get_stats(&self) -> &SketchStats {
        &self.stats
    }

    /// Parameter tuning recommendations based on observed data
    pub fn recommend_parameters(&self) -> ParameterRecommendations {
        let current_error = self.estimate_error_rate();
        let memory_efficiency =
            self.stats.memory_usage_bytes as f64 / self.stats.total_kmers as f64;

        ParameterRecommendations {
            current_error_rate: current_error,
            memory_per_kmer: memory_efficiency,
            recommended_hll_precision: self.recommend_hll_precision(current_error),
            recommended_l0_samples: self.recommend_l0_size(),
            memory_usage_ok: self.stats.memory_usage_bytes
                < (self.config.memory_limit_mb * 1024 * 1024),
        }
    }

    fn estimate_error_rate(&self) -> f64 {
        // Estimate error based on HyperLogLog theory
        1.04 / (self.hyperloglog.num_buckets as f64).sqrt()
    }

    fn recommend_hll_precision(&self, current_error: f64) -> u8 {
        if current_error > 0.05 {
            (self.config.hll_precision + 1).min(16)
        } else if current_error < 0.01 {
            (self.config.hll_precision - 1).max(10)
        } else {
            self.config.hll_precision
        }
    }

    fn recommend_l0_size(&self) -> usize {
        let unique_estimate = self.stats.estimated_unique as f64;
        let current_samples = self.stats.l0_samples_count as f64;

        if current_samples / unique_estimate < 0.001 {
            // Too few samples for good abundance estimation
            (self.config.max_l0_samples * 2).min(100_000)
        } else if current_samples / unique_estimate > 0.1 {
            // Too many samples, wasting memory
            (self.config.max_l0_samples / 2).max(1000)
        } else {
            self.config.max_l0_samples
        }
    }
}

impl HyperLogLog {
    fn new(precision: u8) -> Self {
        let num_buckets = 1 << precision; // 2^precision
        let alpha = Self::calculate_alpha(num_buckets);

        Self {
            num_buckets,
            precision,
            buckets: vec![0; num_buckets],
            alpha,
        }
    }

    fn from_buckets(buckets: Vec<u8>, precision: u8) -> Self {
        let num_buckets = buckets.len();
        let alpha = Self::calculate_alpha(num_buckets);

        Self {
            num_buckets,
            precision,
            buckets,
            alpha,
        }
    }

    fn add(&mut self, hash: u64) {
        // Use first p bits for bucket selection
        let bucket_mask = (1u64 << self.precision) - 1;
        let bucket = (hash & bucket_mask) as usize;

        // Use remaining bits for leading zero count
        let remaining_bits = hash >> self.precision;
        let leading_zeros = if remaining_bits == 0 {
            64 - self.precision
        } else {
            remaining_bits.leading_zeros() as u8 + 1
        };

        // Update bucket with maximum leading zeros seen
        self.buckets[bucket] = self.buckets[bucket].max(leading_zeros);
    }

    fn estimate(&self) -> u64 {
        // Calculate raw estimate
        let sum: f64 = self
            .buckets
            .iter()
            .map(|&bucket_val| 2.0_f64.powi(-(bucket_val as i32)))
            .sum();

        let raw_estimate = self.alpha * (self.num_buckets as f64).powi(2) / sum;

        // Apply bias correction
        if raw_estimate <= 2.5 * self.num_buckets as f64 {
            // Small range correction
            let zeros = self.buckets.iter().filter(|&&x| x == 0).count();
            if zeros != 0 {
                (self.num_buckets as f64 * (self.num_buckets as f64 / zeros as f64).ln()) as u64
            } else {
                raw_estimate as u64
            }
        } else if raw_estimate <= (1.0 / 30.0) * (1u64 << 32) as f64 {
            // Intermediate range - no correction
            raw_estimate as u64
        } else {
            // Large range correction
            let corrected =
                -((1u64 << 32) as f64) * (1.0 - raw_estimate / (1u64 << 32) as f64).ln();
            corrected as u64
        }
    }

    fn merge(&mut self, other: &Self) -> Result<()> {
        if self.precision != other.precision {
            return Err(anyhow::anyhow!(
                "Cannot merge HyperLogLog with different precisions"
            ));
        }

        for (i, &other_val) in other.buckets.iter().enumerate() {
            self.buckets[i] = self.buckets[i].max(other_val);
        }

        Ok(())
    }

    fn calculate_alpha(m: usize) -> f64 {
        match m {
            16 => 0.673,
            32 => 0.697,
            64 => 0.709,
            _ => 0.7213 / (1.0 + 1.079 / m as f64),
        }
    }

    fn memory_usage(&self) -> usize {
        self.buckets.len() + std::mem::size_of::<Self>()
    }

    /// Alias for estimate() method for compatibility
    fn estimate_cardinality(&self) -> u64 {
        self.estimate()
    }
}

impl L0Sampler {
    fn new(max_samples: usize) -> Self {
        Self {
            samples: HashMap::with_capacity(max_samples),
            max_samples,
            threshold: 0.0,
            rng_state: fastrand::u64(..), // Random seed
        }
    }

    fn from_samples(samples: HashMap<u64, f64>, threshold: f64, max_samples: usize) -> Self {
        Self {
            samples,
            max_samples,
            threshold,
            rng_state: fastrand::u64(..),
        }
    }

    fn add(&mut self, item: u64) {
        // Generate priority using hash-based random number
        let priority = self.hash_to_priority(item);

        if self.samples.len() < self.max_samples {
            // Still have space, add directly
            self.samples.insert(item, priority);
            if self.samples.len() == self.max_samples {
                // Now full, set threshold to minimum priority
                self.threshold = self.samples.values().cloned().fold(1.0, f64::min);
            }
        } else if priority > self.threshold {
            // Replace item with minimum priority
            if let Some((&min_item, _)) = self
                .samples
                .iter()
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            {
                self.samples.remove(&min_item);
                self.samples.insert(item, priority);

                // Update threshold
                self.threshold = self.samples.values().cloned().fold(1.0, f64::min);
            }
        }
    }

    fn merge(&mut self, other: &Self) {
        // Combine samples from both samplers
        let mut combined_samples = self.samples.clone();
        combined_samples.extend(other.samples.iter().map(|(&k, &v)| (k, v)));

        // Keep only top samples
        if combined_samples.len() > self.max_samples {
            let mut samples_vec: Vec<_> = combined_samples.into_iter().collect();
            samples_vec.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());
            samples_vec.truncate(self.max_samples);

            self.samples = samples_vec.into_iter().collect();
            self.threshold = self.samples.values().cloned().fold(1.0, f64::min);
        } else {
            self.samples = combined_samples;
            self.threshold = if self.samples.len() == self.max_samples {
                self.samples.values().cloned().fold(1.0, f64::min)
            } else {
                0.0
            };
        }
    }

    fn get_samples_with_weights(&self) -> HashMap<u64, f64> {
        // Convert priorities to abundance estimates
        self.samples
            .iter()
            .map(|(&item, &_priority)| {
                // Abundance is inversely related to priority threshold
                let abundance = if self.threshold > 0.0 {
                    1.0 / self.threshold
                } else {
                    1.0
                };
                (item, abundance)
            })
            .collect()
    }

    fn get_sample_weight(&self, item: u64) -> Option<f64> {
        self.samples.get(&item).map(|&_priority| {
            if self.threshold > 0.0 {
                1.0 / self.threshold
            } else {
                1.0
            }
        })
    }

    fn hash_to_priority(&mut self, item: u64) -> f64 {
        // Use multiplicative hashing for consistent priorities
        self.rng_state = self.rng_state.wrapping_mul(item).wrapping_add(1);
        let hash_val = self.rng_state;

        // Convert to [0, 1] range
        (hash_val as f64) / (u64::MAX as f64)
    }

    fn memory_usage(&self) -> usize {
        self.samples.len() * (std::mem::size_of::<u64>() + std::mem::size_of::<f64>())
            + std::mem::size_of::<Self>()
    }
}

#[derive(Clone)]
pub struct SketchExport {
    hll_buckets: Vec<u8>,
    hll_precision: u8,
    l0_samples: HashMap<u64, f64>,
    l0_threshold: f64,
    config: SketchConfig,
    stats: SketchStats,
}

pub struct ParameterRecommendations {
    pub current_error_rate: f64,
    pub memory_per_kmer: f64,
    pub recommended_hll_precision: u8,
    pub recommended_l0_samples: usize,
    pub memory_usage_ok: bool,
}

/// Parallel processing helper for multiple FASTQ files
pub fn process_multiple_fastq_parallel(
    fastq_files: &[String],
    config: SketchConfig,
) -> Result<HybridAbundanceEstimator> {
    use std::sync::Arc;

    let shared_config = Arc::new(config);

    // Process each file in parallel
    let partial_sketches: Result<Vec<_>> = fastq_files
        .par_iter()
        .map(|file_path| -> Result<HybridAbundanceEstimator> {
            let mut sketch = HybridAbundanceEstimator::new((*shared_config).clone());

            let reader = bio::io::fastq::Reader::from_file(file_path)?;
            for record_result in reader.records() {
                let record = record_result?;
                let sequence = std::str::from_utf8(record.seq())?;
                sketch.add_sequence(sequence)?;
            }

            Ok(sketch)
        })
        .collect();

    let sketches = partial_sketches?;

    // Merge all partial sketches
    let merge_result = HybridAbundanceEstimator::merge_sketches(sketches)?;

    println!(
        "Merged {} sketches:",
        merge_result.merge_stats.sketches_merged
    );
    println!(
        "  Total samples: {} -> {}",
        merge_result.merge_stats.total_samples_before, merge_result.merge_stats.total_samples_after
    );
    println!(
        "  Final unique estimate: {}",
        merge_result.merge_stats.estimated_unique_after
    );

    Ok(merge_result.merged_sketch)
}

/// Integration with existing pipeline
impl L0Sampler {
    /// Upgrade existing L0Sampler to hybrid estimator
    pub fn upgrade_to_hybrid(&self, hll_precision: u8) -> HybridAbundanceEstimator {
        let config = SketchConfig {
            hll_precision,
            max_l0_samples: self.max_samples,
            k: 21, // Default k-mer size
            parallel_processing: true,
            memory_limit_mb: 500, // Default 500MB limit
        };

        let mut hybrid = HybridAbundanceEstimator::new(config);

        // Add existing samples to both sketches
        for &sample in self.samples.keys() {
            hybrid.add_kmer(sample);
        }

        hybrid
    }
}

/// Enhanced sketch output with abundance information
pub fn write_enhanced_sketch(
    sketch: &mut HybridAbundanceEstimator,
    output_path: &str,
) -> Result<()> {
    use std::fs::File;
    use std::io::Write;

    let mut file = File::create(output_path)?;

    // Write header with metadata
    writeln!(file, "# Hybrid Abundance Sketch")?;
    writeln!(
        file,
        "# Estimated unique k-mers: {}",
        sketch.estimate_unique_count()
    )?;
    writeln!(file, "# L0 samples: {}", sketch.stats.l0_samples_count)?;
    writeln!(
        file,
        "# Memory usage: {} bytes",
        sketch.stats.memory_usage_bytes
    )?;
    writeln!(file, "# K-mer size: {}", sketch.config.k)?;
    writeln!(file, "# HLL precision: {}", sketch.config.hll_precision)?;
    writeln!(file)?;

    // Write samples with abundance estimates
    writeln!(file, "kmer_hash\testimated_abundance")?;

    let samples = sketch.get_representative_samples();
    let mut sorted_samples: Vec<_> = samples.into_iter().collect();
    sorted_samples.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());

    for (hash, abundance) in sorted_samples {
        writeln!(file, "{hash}\t{abundance:.6}")?;
    }

    println!("Enhanced sketch written to {output_path}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hyperloglog_basic() {
        let mut hll = HyperLogLog::new(12);

        // Add known set of items
        for i in 0..1000u64 {
            hll.add(i);
        }

        let estimate = hll.estimate();

        // Should be close to 1000 with some error
        assert!(estimate > 800 && estimate < 1200);
        println!("HLL estimate: {estimate} (expected ~1000)");
    }

    #[test]
    fn test_l0_sampler() {
        let mut sampler = L0Sampler::new(100);

        // Add more items than capacity
        for i in 0..1000u64 {
            sampler.add(i);
        }

        // Should retain exactly max_samples
        assert_eq!(sampler.samples.len(), 100);

        // Threshold should be positive
        assert!(sampler.threshold > 0.0);
    }

    #[test]
    fn test_hybrid_estimator() {
        let config = SketchConfig {
            hll_precision: 12,
            max_l0_samples: 1000,
            k: 21,
            parallel_processing: false,
            memory_limit_mb: 100,
        };

        let mut estimator = HybridAbundanceEstimator::new(config);

        // Add test sequence
        let test_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        estimator.add_sequence(test_seq).unwrap();

        let unique_count = estimator.estimate_unique_count();
        let samples = estimator.get_representative_samples();

        assert!(unique_count > 0);
        assert!(!samples.is_empty());

        println!(
            "Unique estimate: {}, Samples: {}",
            unique_count,
            samples.len()
        );
    }

    #[test]
    fn test_sketch_merging() {
        let config = SketchConfig {
            hll_precision: 10,
            max_l0_samples: 100,
            k: 15,
            parallel_processing: false,
            memory_limit_mb: 50,
        };

        // Create multiple sketches
        let mut sketch1 = HybridAbundanceEstimator::new(config.clone());
        let mut sketch2 = HybridAbundanceEstimator::new(config.clone());

        sketch1
            .add_sequence("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
            .unwrap();
        sketch2
            .add_sequence("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
            .unwrap();

        let mut result = HybridAbundanceEstimator::merge_sketches(vec![sketch1, sketch2]).unwrap();

        assert_eq!(result.merge_stats.sketches_merged, 2);
        assert!(result.merged_sketch.estimate_unique_count() > 0);
    }

    #[test]
    fn test_parameter_recommendations() {
        let config = SketchConfig {
            hll_precision: 10,
            max_l0_samples: 1000,
            k: 21,
            parallel_processing: false,
            memory_limit_mb: 100,
        };

        let mut estimator = HybridAbundanceEstimator::new(config);

        // Add enough data to get meaningful recommendations
        for i in 0..10000 {
            estimator.add_kmer(i);
        }

        let recommendations = estimator.recommend_parameters();

        assert!(recommendations.current_error_rate > 0.0);
        assert!(recommendations.memory_per_kmer > 0.0);
        println!(
            "Error rate: {:.4}, Memory per k-mer: {:.2} bytes",
            recommendations.current_error_rate, recommendations.memory_per_kmer
        );
    }

    /// Test HyperLogLog cardinality estimation accuracy within theoretical bounds
    #[test]
    fn test_hyperloglog_cardinality_estimation() {
        let mut hll = HyperLogLog::new(12); // precision=12 gives ~1.04/sqrt(4096) ≈ 1.6% error

        // Add exactly 10,000 unique items
        let true_cardinality = 10000;
        for i in 0..true_cardinality {
            hll.add(i as u64);
        }

        let estimated = hll.estimate_cardinality();
        let relative_error =
            ((estimated as f64 - true_cardinality as f64) / true_cardinality as f64).abs();

        // HyperLogLog should be within ~3% for precision=12 (3 standard deviations)
        assert!(
            relative_error < 0.03,
            "HyperLogLog estimation error {:.3}% exceeds 3% threshold. Estimated: {}, True: {}",
            relative_error * 100.0,
            estimated,
            true_cardinality
        );

        // Verify estimate is reasonable (not zero or way off)
        assert!(estimated > true_cardinality as u64 / 2);
        assert!(estimated < true_cardinality as u64 * 2);
    }

    /// Test L0 sampler respects memory limits and maintains sample quality
    #[test]
    fn test_l0_sampler_memory_bounds() {
        let max_samples = 500;
        let mut sampler = L0Sampler::new(max_samples);

        // Add many more items than max_samples to test memory bounds
        let total_items = 5000;
        for i in 0..total_items {
            sampler.add(i);
        }

        // Verify memory bounds are respected
        assert!(
            sampler.samples.len() <= max_samples,
            "L0Sampler exceeded memory limit: {} samples > {} max",
            sampler.samples.len(),
            max_samples
        );

        // Verify sampler maintained some samples
        assert!(
            sampler.samples.len() > 0,
            "L0Sampler should retain some samples"
        );

        // Test threshold updates - should increase as more items are added
        let initial_threshold = sampler.threshold;

        // Add more items
        for i in total_items..(total_items + 1000) {
            sampler.add(i);
        }

        // Threshold should increase (become more selective)
        assert!(
            sampler.threshold >= initial_threshold,
            "L0Sampler threshold should increase with more items: {} >= {}",
            sampler.threshold,
            initial_threshold
        );

        // Still respect memory bounds
        assert!(sampler.samples.len() <= max_samples);
    }
}
