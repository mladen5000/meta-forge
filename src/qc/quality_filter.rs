//! Quality-based read filtering and trimming
//!
//! Implements:
//! - Sliding window quality trimming
//! - Per-base quality filtering
//! - Average quality thresholds
//! - Length-based filtering


/// Quality filtering configuration
#[derive(Debug, Clone)]
pub struct QualityFilterConfig {
    /// Minimum Phred quality score (default: 20 = Q20)
    pub min_quality: u8,

    /// Sliding window size for quality trimming (default: 4)
    pub window_size: usize,

    /// Minimum average quality in sliding window
    pub min_window_quality: f64,

    /// Minimum read length after trimming (default: 50bp)
    pub min_length: usize,

    /// Minimum average quality for entire read
    pub min_avg_quality: f64,

    /// Quality encoding (default: Sanger/Illumina 1.8+ with offset 33)
    pub quality_offset: u8,
}

impl Default for QualityFilterConfig {
    fn default() -> Self {
        Self {
            min_quality: 20,           // Q20 standard
            window_size: 4,            // 4bp sliding window
            min_window_quality: 20.0,  // Q20 average in window
            min_length: 50,            // Minimum 50bp
            min_avg_quality: 25.0,     // Q25 average for whole read
            quality_offset: 33,        // Sanger/Illumina 1.8+ encoding
        }
    }
}

/// Quality filter for reads
pub struct QualityFilter {
    pub config: QualityFilterConfig,
}

impl QualityFilter {
    pub fn new(config: QualityFilterConfig) -> Self {
        Self { config }
    }

    /// Trim low-quality bases from read ends using sliding window
    pub fn trim_quality(&self, sequence: &str, quality: &[u8]) -> Option<(usize, usize)> {
        if sequence.len() != quality.len() {
            return None;
        }

        let length = sequence.len();
        if length < self.config.window_size {
            return None;
        }

        // Convert quality scores from ASCII to numeric
        let qualities: Vec<u8> = quality
            .iter()
            .map(|&q| q.saturating_sub(self.config.quality_offset))
            .collect();

        // Find trim start position (trim from 5' end)
        let start = self.find_trim_start(&qualities);

        // Find trim end position (trim from 3' end)
        let end = self.find_trim_end(&qualities);

        if end <= start || (end - start) < self.config.min_length {
            None
        } else {
            Some((start, end))
        }
    }

    /// Find start position by scanning from 5' end with sliding window
    fn find_trim_start(&self, qualities: &[u8]) -> usize {
        for i in 0..=(qualities.len().saturating_sub(self.config.window_size)) {
            let window = &qualities[i..i + self.config.window_size];
            let avg_quality = window.iter().map(|&q| q as f64).sum::<f64>()
                / self.config.window_size as f64;

            if avg_quality >= self.config.min_window_quality {
                return i;
            }
        }
        qualities.len() // No good region found
    }

    /// Find end position by scanning from 3' end with sliding window
    fn find_trim_end(&self, qualities: &[u8]) -> usize {
        for i in (0..=(qualities.len().saturating_sub(self.config.window_size))).rev() {
            let window = &qualities[i..i + self.config.window_size];
            let avg_quality = window.iter().map(|&q| q as f64).sum::<f64>()
                / self.config.window_size as f64;

            if avg_quality >= self.config.min_window_quality {
                return i + self.config.window_size;
            }
        }
        0 // No good region found
    }

    /// Check if read passes average quality threshold
    pub fn passes_quality_threshold(&self, quality: &[u8]) -> bool {
        if quality.is_empty() {
            return false;
        }

        let qualities: Vec<u8> = quality
            .iter()
            .map(|&q| q.saturating_sub(self.config.quality_offset))
            .collect();

        let avg_quality = qualities.iter().map(|&q| q as f64).sum::<f64>()
            / qualities.len() as f64;

        avg_quality >= self.config.min_avg_quality
    }

    /// Check if read length passes minimum threshold
    pub fn passes_length_threshold(&self, length: usize) -> bool {
        length >= self.config.min_length
    }

    /// Filter and trim a read, returning trimmed sequence and quality
    pub fn filter_read(
        &self,
        sequence: &str,
        quality: &[u8],
    ) -> Option<(String, Vec<u8>)> {
        // First check average quality
        if !self.passes_quality_threshold(quality) {
            return None;
        }

        // Trim based on sliding window
        if let Some((start, end)) = self.trim_quality(sequence, quality) {
            let trimmed_seq = sequence[start..end].to_string();
            let trimmed_qual = quality[start..end].to_vec();

            // Check trimmed length
            if self.passes_length_threshold(trimmed_seq.len()) {
                return Some((trimmed_seq, trimmed_qual));
            }
        }

        None
    }

    /// Calculate quality statistics for a read
    pub fn quality_stats(&self, quality: &[u8]) -> QualityStats {
        if quality.is_empty() {
            return QualityStats::default();
        }

        let qualities: Vec<u8> = quality
            .iter()
            .map(|&q| q.saturating_sub(self.config.quality_offset))
            .collect();

        let min = *qualities.iter().min().unwrap();
        let max = *qualities.iter().max().unwrap();
        let mean = qualities.iter().map(|&q| q as f64).sum::<f64>()
            / qualities.len() as f64;

        // Calculate median
        let mut sorted = qualities.clone();
        sorted.sort_unstable();
        let median = if sorted.len() % 2 == 0 {
            (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) as f64 / 2.0
        } else {
            sorted[sorted.len() / 2] as f64
        };

        // Calculate Q20 and Q30 percentages
        let q20_count = qualities.iter().filter(|&&q| q >= 20).count();
        let q30_count = qualities.iter().filter(|&&q| q >= 30).count();

        QualityStats {
            min_quality: min,
            max_quality: max,
            mean_quality: mean,
            median_quality: median,
            q20_percentage: (q20_count as f64 / qualities.len() as f64) * 100.0,
            q30_percentage: (q30_count as f64 / qualities.len() as f64) * 100.0,
        }
    }
}

/// Quality statistics for a read
#[derive(Debug, Clone, Default)]
pub struct QualityStats {
    pub min_quality: u8,
    pub max_quality: u8,
    pub mean_quality: f64,
    pub median_quality: f64,
    pub q20_percentage: f64,
    pub q30_percentage: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_trimming() {
        let config = QualityFilterConfig {
            min_length: 30, // Lower for test
            ..Default::default()
        };
        let filter = QualityFilter::new(config);

        // Quality scores: ASCII with offset 33
        // Low quality at ends: !!!! (Q0-1), good in middle: IIII (Q40)
        let sequence = format!("AAAA{}GGGG{}TTTT", "T".repeat(30), "C".repeat(30));
        let mut quality = Vec::new();
        quality.extend_from_slice(b"!!!!");
        quality.extend(vec![b'I'; 30]);
        quality.extend_from_slice(b"IIII");
        quality.extend(vec![b'I'; 30]);
        quality.extend_from_slice(b"!!!!");

        let result = filter.trim_quality(&sequence, &quality);
        assert!(result.is_some());

        let (start, end) = result.unwrap();
        // Should trim low quality ends
        assert!(end - start >= 30);
        assert!(!sequence[start..end].starts_with("AAAA"));
        assert!(!sequence[start..end].ends_with("TTTT"));
    }

    #[test]
    fn test_quality_threshold() {
        let config = QualityFilterConfig {
            min_avg_quality: 30.0,
            ..Default::default()
        };
        let filter = QualityFilter::new(config);

        // High quality: J = Q41
        let high_quality = b"JJJJJJJJ";
        assert!(filter.passes_quality_threshold(high_quality));

        // Low quality: # = Q2
        let low_quality = b"########";
        assert!(!filter.passes_quality_threshold(low_quality));
    }

    #[test]
    fn test_length_filtering() {
        let config = QualityFilterConfig {
            min_length: 50,
            ..Default::default()
        };
        let filter = QualityFilter::new(config);

        assert!(filter.passes_length_threshold(100));
        assert!(filter.passes_length_threshold(50));
        assert!(!filter.passes_length_threshold(49));
    }

    #[test]
    fn test_filter_read() {
        let config = QualityFilterConfig::default();
        let filter = QualityFilter::new(config);

        // Good read: 60bp with Q35 quality
        let sequence = "A".repeat(60);
        let quality = vec![b'D'; 60]; // D = Q35

        let result = filter.filter_read(&sequence, &quality);
        assert!(result.is_some());

        let (trimmed_seq, trimmed_qual) = result.unwrap();
        assert!(trimmed_seq.len() >= 50);
        assert_eq!(trimmed_seq.len(), trimmed_qual.len());
    }

    #[test]
    fn test_quality_stats() {
        let config = QualityFilterConfig::default();
        let filter = QualityFilter::new(config);

        // Mixed quality
        let quality = b"####IIIIJJJJ"; // Q2-3, Q40, Q41

        let stats = filter.quality_stats(quality);
        assert!(stats.mean_quality > 20.0);
        assert!(stats.q20_percentage > 50.0);
    }
}
