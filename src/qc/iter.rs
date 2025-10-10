//! Iterator-based QC processing
//!
//! Provides idiomatic iterator chains for read quality control

use super::{QCPipeline, QCPipelineConfig};
use crate::core::data_structures::CorrectedRead;

/// Iterator adapter for QC filtering
pub struct QCFilterIter<I> {
    inner: I,
    pipeline: QCPipeline,
}

impl<I> QCFilterIter<I>
where
    I: Iterator<Item = CorrectedRead>,
{
    pub fn new(iter: I, config: QCPipelineConfig) -> Self {
        Self {
            inner: iter,
            pipeline: QCPipeline::new(config),
        }
    }

    /// Get QC statistics (consumes remaining items if needed)
    pub fn stats(mut self) -> super::QCStats {
        // Process any remaining items to get complete stats
        for _ in self.by_ref() {}
        self.pipeline.stats()
    }
}

impl<I> Iterator for QCFilterIter<I>
where
    I: Iterator<Item = CorrectedRead>,
{
    type Item = CorrectedRead;

    fn next(&mut self) -> Option<Self::Item> {
        // Keep trying until we find a read that passes QC
        for read in self.inner.by_ref() {
            if let Some(filtered) = self.pipeline.process_read(&read) {
                return Some(filtered);
            }
        }
        None
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        // Lower bound is 0 (all could be filtered), upper bound from inner
        let (_, upper) = self.inner.size_hint();
        (0, upper)
    }
}

/// Extension trait for adding QC to iterators
pub trait QCIterExt: Iterator<Item = CorrectedRead> + Sized {
    /// Apply quality control to reads
    fn qc_filter(self, config: QCPipelineConfig) -> QCFilterIter<Self> {
        QCFilterIter::new(self, config)
    }

    /// Apply quality control with default config
    fn qc_default(self) -> QCFilterIter<Self> {
        QCFilterIter::new(self, QCPipelineConfig::default())
    }
}

// Implement for all iterators of CorrectedRead
impl<I> QCIterExt for I where I: Iterator<Item = CorrectedRead> {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::data_structures::CorrectionMetadata;

    fn create_test_read(id: usize, quality: u8) -> CorrectedRead {
        CorrectedRead {
            id,
            original: "A".repeat(100),
            corrected: "A".repeat(100),
            corrections: Vec::new(),
            quality_scores: vec![quality; 100],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.0,
                context_window: 0,
                correction_time_ms: 0,
            },
            kmer_hash_cache: Vec::new(),
        }
    }

    #[test]
    fn test_qc_iter() {
        let reads = vec![
            create_test_read(1, b'I'), // Good (Q40)
            create_test_read(2, b'#'), // Bad (Q2)
            create_test_read(3, b'I'), // Good (Q40)
        ];

        let filtered: Vec<_> = reads.into_iter().qc_default().collect();

        assert_eq!(filtered.len(), 2); // Only 2 good reads
        assert_eq!(filtered[0].id, 1);
        assert_eq!(filtered[1].id, 3);
    }

    #[test]
    fn test_qc_stats() {
        let reads = vec![
            create_test_read(1, b'I'),
            create_test_read(2, b'#'),
            create_test_read(3, b'I'),
        ];

        let stats = reads.into_iter().qc_default().stats();

        assert_eq!(stats.reads_input, 3);
        assert_eq!(stats.reads_passed, 2);
        assert_eq!(stats.reads_failed, 1);
    }
}
