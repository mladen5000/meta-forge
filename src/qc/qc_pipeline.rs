//! Integrated QC Pipeline
//!
//! Combines quality filtering, adapter trimming, and statistics

use ahash::AHashMap;

use super::qc_stats::FailureReason;
use super::{
    AdapterConfig, AdapterMatch, AdapterTrimmer, QCStats, QualityFilter, QualityFilterConfig,
};
use crate::core::data_structures::CorrectedRead;
use crate::utils::genomic_validator::{GenomicDataValidator, ValidationThresholds};

/// Complete QC pipeline configuration
#[derive(Debug, Clone)]
pub struct QCPipelineConfig {
    /// Enable quality filtering
    pub enable_quality_filter: bool,

    /// Enable adapter trimming
    pub enable_adapter_trimming: bool,

    /// Enable genomic validation
    pub enable_genomic_validation: bool,

    /// Quality filter configuration
    pub quality_config: QualityFilterConfig,

    /// Adapter trimmer configuration
    pub adapter_config: AdapterConfig,

    /// Genomic validation thresholds
    pub validation_thresholds: Option<ValidationThresholds>,

    /// Verbose output
    pub verbose: bool,
}

impl Default for QCPipelineConfig {
    fn default() -> Self {
        Self {
            enable_quality_filter: true,
            enable_adapter_trimming: true,
            enable_genomic_validation: false, // Disabled by default - enable when needed
            quality_config: QualityFilterConfig::default(),
            adapter_config: AdapterConfig::default(),
            validation_thresholds: None, // Uses defaults if None
            verbose: false,              // Disabled for performance
        }
    }
}

/// Complete QC pipeline
pub struct QCPipeline {
    quality_filter: QualityFilter,
    adapter_trimmer: AdapterTrimmer,
    genomic_validator: Option<GenomicDataValidator>,
    config: QCPipelineConfig,
    stats: QCStats,
    reads_processed: usize,
}

impl QCPipeline {
    pub fn new(config: QCPipelineConfig) -> Self {
        let genomic_validator = if config.enable_genomic_validation {
            // Note: GenomicDataValidator uses default thresholds
            // Custom thresholds support could be added via a with_thresholds() method
            Some(GenomicDataValidator::new())
        } else {
            None
        };

        Self {
            quality_filter: QualityFilter::new(config.quality_config.clone()),
            adapter_trimmer: AdapterTrimmer::new(config.adapter_config.clone()),
            genomic_validator,
            config,
            stats: QCStats::new(),
            reads_processed: 0,
        }
    }

    /// Process a single read through the QC pipeline (optimized zero-copy version)
    pub fn process_read(&mut self, read: &CorrectedRead) -> Option<CorrectedRead> {
        self.reads_processed += 1;

        // Work with slices (zero-copy) until we absolutely need to allocate
        let seq = read.corrected.as_str();
        let qual = &read.quality_scores[..];

        // Record input
        self.stats.record_input(seq.len(), qual);

        // Fast path: check average quality first (cheapest check)
        if self.config.enable_quality_filter && !self.quality_filter.passes_quality_threshold(qual)
        {
            self.stats.record_failed(FailureReason::Quality);
            return None;
        }

        // Fast path: check minimum length before any processing
        if seq.len() < self.config.quality_config.min_length {
            self.stats.record_failed(FailureReason::Length);
            return None;
        }

        // Step 0: Genomic validation (if enabled) - work with slices
        if let Some(validator) = &mut self.genomic_validator {
            let validation_result = validator.validate_sequence(seq, Some(qual));
            if !validation_result.passed {
                self.stats.record_failed(FailureReason::Quality);
                return None;
            }
        }

        // Determine trim positions (still zero-copy, just calculating indices)
        let (trim_start, trim_end) = if self.config.enable_quality_filter {
            match self.quality_filter.trim_quality(seq, qual) {
                Some((start, end)) => {
                    let quality_trimmed_bases = seq.len() - (end - start);
                    if quality_trimmed_bases > 0 {
                        self.stats.record_quality_trimming(quality_trimmed_bases);
                    }
                    (start, end)
                }
                None => {
                    self.stats.record_failed(FailureReason::Quality);
                    return None;
                }
            }
        } else {
            (0, seq.len())
        };

        // Check adapter trimming (only if needed)
        let (final_start, final_end) = if self.config.enable_adapter_trimming {
            // Adapter trimming on the slice
            let slice_seq = &seq[trim_start..trim_end];
            let slice_qual = &qual[trim_start..trim_end];
            let (trimmed_seq, _, adapter_match) = self
                .adapter_trimmer
                .trim_adapter_with_quality(slice_seq, slice_qual);

            // Record adapter if found
            if let Some(adapter) = adapter_match {
                let bases_trimmed = slice_seq.len() - trimmed_seq.len();
                self.stats
                    .record_adapter_trimming(adapter.adapter.clone(), bases_trimmed);
            }

            // Calculate new positions relative to original
            if trimmed_seq.len() < slice_seq.len() {
                let adapter_trim = slice_seq.len() - trimmed_seq.len();
                (trim_start, trim_end - adapter_trim)
            } else {
                (trim_start, trim_end)
            }
        } else {
            (trim_start, trim_end)
        };

        // Final length check
        if (final_end - final_start) < self.config.quality_config.min_length {
            self.stats.record_failed(FailureReason::Length);
            return None;
        }

        // Record passed read
        self.stats
            .record_passed(final_end - final_start, &qual[final_start..final_end]);

        // Only NOW do we allocate (when we know the read passes)
        // Zero-copy if no trimming needed
        if final_start == 0 && final_end == seq.len() {
            // No trimming needed - clone once
            Some(read.clone())
        } else {
            // Trimming needed - allocate only trimmed portion
            Some(CorrectedRead {
                id: read.id,
                original: read.original.clone(),
                corrected: seq[final_start..final_end].to_string(),
                corrections: read.corrections.clone(),
                quality_scores: qual[final_start..final_end].to_vec(),
                correction_metadata: read.correction_metadata.clone(),
                kmer_hash_cache: Vec::new(), // Clear cache as sequence changed
            })
        }
    }

    /// Process multiple reads (optimized with rayon parallel processing)
    pub fn process_reads(&mut self, reads: &[CorrectedRead]) -> Vec<CorrectedRead> {
        self.process_reads_with_progress(reads, |_, _| {})
    }

    /// Process multiple reads with progress callback (ultra-fast zero-copy version)
    pub fn process_reads_with_progress<F>(
        &mut self,
        reads: &[CorrectedRead],
        mut progress_callback: F,
    ) -> Vec<CorrectedRead>
    where
        F: FnMut(usize, usize),
    {
        let total = reads.len();

        // Sequential processing with progress updates and full statistics tracking
        // Note: We use sequential instead of parallel because FnMut cannot be shared across threads
        let results: Vec<CorrectedRead> = reads
            .iter()
            .enumerate()
            .filter_map(|(idx, read)| {
                // Use the full process_read method that tracks all statistics
                let result = self.process_read(read);

                // Update progress periodically
                let count = idx + 1;
                if count % 100 == 0 || count == total {
                    progress_callback(count, total);
                }

                result
            })
            .collect();

        // Final progress update
        if total > 0 {
            progress_callback(total, total);
        }

        results
    }

    /// Get QC statistics
    pub fn stats(&mut self) -> QCStats {
        self.stats.finalize();
        self.stats.clone()
    }

    /// Get mutable reference to stats (for updates)
    pub fn stats_mut(&mut self) -> &mut QCStats {
        &mut self.stats
    }

    /// Print comprehensive debug statistics for a read
    fn print_debug_stats(
        &self,
        read: &CorrectedRead,
        original_length: usize,
        quality_trimmed: usize,
        adapter_info: Option<(usize, &AdapterMatch)>,
        final_length: Option<usize>,
        status: &str,
    ) {
        let avg_quality = if !read.quality_scores.is_empty() {
            read.quality_scores
                .iter()
                .map(|&q| (q - 33) as f64)
                .sum::<f64>()
                / read.quality_scores.len() as f64
        } else {
            0.0
        };

        let min_quality = read
            .quality_scores
            .iter()
            .map(|&q| q - 33)
            .min()
            .unwrap_or(0);
        let max_quality = read
            .quality_scores
            .iter()
            .map(|&q| q - 33)
            .max()
            .unwrap_or(0);

        eprintln!("\n╔══════════════════════════════════════════════════════════════════════");
        eprintln!("║ 📊 DEBUG READ STATS #{}", self.reads_processed);
        eprintln!("╠══════════════════════════════════════════════════════════════════════");
        eprintln!("║ Read ID:           {}", read.id);
        eprintln!("║ Status:            {}", status);
        eprintln!("║ ");
        eprintln!("║ LENGTH METRICS:");
        eprintln!("║   Original length:     {} bp", original_length);

        if let Some((adapter_trimmed, adapter)) = adapter_info {
            eprintln!(
                "║   After adapter trim:  {} bp (trimmed {} bp)",
                original_length - adapter_trimmed,
                adapter_trimmed
            );
            eprintln!(
                "║     ↳ Adapter:         {} at position {}",
                &adapter.adapter[..adapter.adapter.len().min(15)],
                adapter.position
            );
            eprintln!(
                "║     ↳ Adapter length:  {} bp (matched with {:.1}% error)",
                adapter.length,
                adapter.error_rate * 100.0
            );
            eprintln!("║     ↳ Mismatches:      {}", adapter.mismatches);
        }

        if quality_trimmed > 0 {
            eprintln!("║   Quality trimmed:     {} bp", quality_trimmed);
        }

        if let Some(final_len) = final_length {
            eprintln!("║   Final length:        {} bp", final_len);
            let total_removed = original_length - final_len;
            let removal_pct = (total_removed as f64 / original_length as f64) * 100.0;
            eprintln!(
                "║   Total removed:       {} bp ({:.1}%)",
                total_removed, removal_pct
            );
        }

        eprintln!("║ ");
        eprintln!("║ QUALITY METRICS:");
        eprintln!("║   Average quality:     Q{:.1}", avg_quality);
        eprintln!("║   Min quality:         Q{}", min_quality);
        eprintln!("║   Max quality:         Q{}", max_quality);
        eprintln!("║   Quality range:       Q{}-Q{}", min_quality, max_quality);

        // GC content
        let gc_count = read
            .corrected
            .chars()
            .filter(|&c| c == 'G' || c == 'C')
            .count();
        let gc_content = if !read.corrected.is_empty() {
            (gc_count as f64 / read.corrected.len() as f64) * 100.0
        } else {
            0.0
        };

        eprintln!("║ ");
        eprintln!("║ SEQUENCE METRICS:");
        eprintln!("║   GC content:          {:.1}%", gc_content);
        eprintln!("║   Corrections made:    {}", read.corrections.len());

        // N bases
        let n_count = read.corrected.chars().filter(|&c| c == 'N').count();
        if n_count > 0 {
            eprintln!(
                "║   N bases:             {} ({:.1}%)",
                n_count,
                (n_count as f64 / read.corrected.len() as f64) * 100.0
            );
        }

        // Base composition
        let a_count = read.corrected.chars().filter(|&c| c == 'A').count();
        let t_count = read.corrected.chars().filter(|&c| c == 'T').count();
        let g_count = read.corrected.chars().filter(|&c| c == 'G').count();
        let c_count = read.corrected.chars().filter(|&c| c == 'C').count();
        eprintln!(
            "║   Base composition:    A={}, T={}, G={}, C={}, N={}",
            a_count, t_count, g_count, c_count, n_count
        );

        // First/Last 20bp preview
        let preview_len = 20.min(read.corrected.len());
        if preview_len > 0 {
            eprintln!("║ ");
            eprintln!("║ SEQUENCE PREVIEW:");
            eprintln!(
                "║   First {}bp:         {}",
                preview_len,
                &read.corrected[..preview_len]
            );
            if read.corrected.len() > preview_len {
                let start = read.corrected.len() - preview_len;
                eprintln!(
                    "║   Last {}bp:          {}",
                    preview_len,
                    &read.corrected[start..]
                );
            }

            // Quality preview
            let qual_preview: Vec<String> = read.quality_scores[..preview_len]
                .iter()
                .map(|&q| format!("Q{}", q - 33))
                .collect();
            eprintln!(
                "║   First {}bp qual:    {}",
                preview_len,
                qual_preview.join(" ")
            );

            if read.quality_scores.len() > preview_len {
                let start = read.quality_scores.len() - preview_len;
                let qual_preview_end: Vec<String> = read.quality_scores[start..]
                    .iter()
                    .map(|&q| format!("Q{}", q - 33))
                    .collect();
                eprintln!(
                    "║   Last {}bp qual:     {}",
                    preview_len,
                    qual_preview_end.join(" ")
                );
            }
        }

        eprintln!("╚══════════════════════════════════════════════════════════════════════\n");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::data_structures::CorrectionMetadata;

    fn create_test_read(id: usize, sequence: &str, quality: Vec<u8>) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: quality,
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.0,
                context_window: 0,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn test_quality_filtering() {
        let config = QCPipelineConfig::default();
        let mut pipeline = QCPipeline::new(config);

        // Good quality read (Q35)
        let good_read = create_test_read(1, "A".repeat(100).as_str(), vec![b'D'; 100]);

        let result = pipeline.process_read(&good_read);
        assert!(result.is_some());

        // Bad quality read (Q5)
        let bad_read = create_test_read(2, "A".repeat(100).as_str(), vec![b'%'; 100]);

        let result = pipeline.process_read(&bad_read);
        assert!(result.is_none());
    }

    #[test]
    fn test_adapter_trimming() {
        let config = QCPipelineConfig::default();
        let mut pipeline = QCPipeline::new(config);

        // Read with adapter (make sequence long enough: 60bp + adapter)
        let base_sequence = "A".repeat(60);
        let sequence = format!(
            "{}{}",
            base_sequence,
            super::super::adapter_trimmer::ILLUMINA_TRUSEQ_R1
        );
        let quality = vec![b'I'; sequence.len()];

        let read = create_test_read(1, &sequence, quality);
        let result = pipeline.process_read(&read);

        assert!(result.is_some());
        let processed = result.unwrap();
        assert_eq!(processed.corrected, base_sequence); // Adapter should be trimmed
        assert!(processed.corrected.len() >= 50); // Should pass length filter
    }

    #[test]
    fn test_quality_and_adapter() {
        let config = QCPipelineConfig::default();
        let mut pipeline = QCPipeline::new(config);

        // Read with low quality at ends and adapter (make long enough: 70bp + adapter)
        let base_sequence = "A".repeat(70);
        let sequence = format!(
            "{}{}",
            base_sequence,
            super::super::adapter_trimmer::ILLUMINA_TRUSEQ_R1
        );
        let mut quality = vec![b'!'; 4]; // Low quality start
        quality.extend(vec![b'I'; 66]); // Good quality middle (70 - 4 = 66)
        quality.extend(vec![
            b'I';
            super::super::adapter_trimmer::ILLUMINA_TRUSEQ_R1.len()
        ]); // Adapter quality

        let read = create_test_read(1, &sequence, quality);
        let result = pipeline.process_read(&read);

        assert!(result.is_some());
        let processed = result.unwrap();
        // Should have trimmed low quality AND adapter
        assert!(processed.corrected.len() < sequence.len());
        assert!(processed.corrected.len() >= 50); // Above min length
        assert!(!processed
            .corrected
            .contains(super::super::adapter_trimmer::ILLUMINA_TRUSEQ_R1)); // No adapter in result
    }

    #[test]
    fn test_statistics() {
        let config = QCPipelineConfig::default();
        let mut pipeline = QCPipeline::new(config);

        let reads = vec![
            create_test_read(1, "A".repeat(100).as_str(), vec![b'I'; 100]), // Good
            create_test_read(2, "A".repeat(100).as_str(), vec![b'#'; 100]), // Bad quality
            create_test_read(3, "A".repeat(30).as_str(), vec![b'I'; 30]),   // Too short
        ];

        let processed = pipeline.process_reads(&reads);
        let stats = pipeline.stats();

        assert_eq!(stats.reads_input, 3);
        assert_eq!(stats.reads_passed, 1);
        assert_eq!(stats.reads_failed, 2);
    }

    #[test]
    fn test_disabled_qc() {
        let config = QCPipelineConfig {
            enable_quality_filter: false,
            enable_adapter_trimming: false,
            ..Default::default()
        };
        let mut pipeline = QCPipeline::new(config);

        // Even bad quality read should pass with QC disabled
        let bad_read = create_test_read(1, "A".repeat(100).as_str(), vec![b'#'; 100]);

        let result = pipeline.process_read(&bad_read);
        assert!(result.is_some());
    }
}
