//! Integrated QC Pipeline
//!
//! Combines quality filtering, adapter trimming, and statistics

use super::qc_stats::FailureReason;
use super::{AdapterConfig, AdapterTrimmer, QCStats, QualityFilter, QualityFilterConfig};
use crate::core::data_structures::CorrectedRead;
use anyhow::Result;

/// Complete QC pipeline configuration
#[derive(Debug, Clone)]
pub struct QCPipelineConfig {
    /// Enable quality filtering
    pub enable_quality_filter: bool,

    /// Enable adapter trimming
    pub enable_adapter_trimming: bool,

    /// Quality filter configuration
    pub quality_config: QualityFilterConfig,

    /// Adapter trimmer configuration
    pub adapter_config: AdapterConfig,

    /// Verbose output
    pub verbose: bool,
}

impl Default for QCPipelineConfig {
    fn default() -> Self {
        Self {
            enable_quality_filter: true,
            enable_adapter_trimming: true,
            quality_config: QualityFilterConfig::default(),
            adapter_config: AdapterConfig::default(),
            verbose: true,
        }
    }
}

/// Complete QC pipeline
pub struct QCPipeline {
    quality_filter: QualityFilter,
    adapter_trimmer: AdapterTrimmer,
    config: QCPipelineConfig,
    stats: QCStats,
}

impl QCPipeline {
    pub fn new(config: QCPipelineConfig) -> Self {
        Self {
            quality_filter: QualityFilter::new(config.quality_config.clone()),
            adapter_trimmer: AdapterTrimmer::new(config.adapter_config.clone()),
            config,
            stats: QCStats::new(),
        }
    }

    /// Process a single read through the QC pipeline
    pub fn process_read(&mut self, read: &CorrectedRead) -> Option<CorrectedRead> {
        let mut sequence = read.corrected.clone();
        let mut quality = read.quality_scores.clone();
        let original_length = sequence.len();

        // Record input
        self.stats.record_input(sequence.len(), &quality);

        // Step 1: Adapter trimming (if enabled)
        if self.config.enable_adapter_trimming {
            let (trimmed_seq, trimmed_qual, adapter_match) = self
                .adapter_trimmer
                .trim_adapter_with_quality(&sequence, &quality);

            if let Some(adapter) = adapter_match {
                let bases_trimmed = sequence.len() - trimmed_seq.len();
                self.stats
                    .record_adapter_trimming(adapter.adapter.clone(), bases_trimmed);

                if self.config.verbose {
                    eprintln!(
                        "  Adapter trimmed from read {}: {} bp ({:.1}% error)",
                        read.id,
                        bases_trimmed,
                        adapter.error_rate * 100.0
                    );
                }
            }

            sequence = trimmed_seq;
            quality = trimmed_qual;
        }

        // Step 2: Quality filtering and trimming (if enabled)
        if self.config.enable_quality_filter {
            // First check average quality
            if !self.quality_filter.passes_quality_threshold(&quality) {
                self.stats.record_failed(FailureReason::Quality);
                return None;
            }

            // Then trim low-quality ends
            if let Some((start, end)) = self.quality_filter.trim_quality(&sequence, &quality) {
                let bases_trimmed = (sequence.len() - (end - start));
                self.stats.record_quality_trimming(bases_trimmed);

                sequence = sequence[start..end].to_string();
                quality = quality[start..end].to_vec();

                if self.config.verbose && bases_trimmed > 0 {
                    eprintln!(
                        "  Quality trimmed from read {}: {} bp",
                        read.id, bases_trimmed
                    );
                }
            } else {
                // Entire read is low quality
                self.stats.record_failed(FailureReason::Quality);
                return None;
            }

            // Final length check
            if !self.quality_filter.passes_length_threshold(sequence.len()) {
                self.stats.record_failed(FailureReason::Length);
                return None;
            }
        }

        // Record passed read
        self.stats.record_passed(sequence.len(), &quality);

        // Create processed read
        Some(CorrectedRead {
            id: read.id,
            original: read.original.clone(),
            corrected: sequence,
            corrections: read.corrections.clone(),
            quality_scores: quality,
            correction_metadata: read.correction_metadata.clone(),
        })
    }

    /// Process multiple reads
    pub fn process_reads(&mut self, reads: &[CorrectedRead]) -> Vec<CorrectedRead> {
        reads
            .iter()
            .filter_map(|read| self.process_read(read))
            .collect()
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
