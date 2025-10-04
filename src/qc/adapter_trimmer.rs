//! Adapter detection and trimming
//!
//! Detects and removes common sequencing adapters:
//! - Illumina TruSeq
//! - Nextera
//! - Custom adapters

use anyhow::Result;

/// Common Illumina adapters
pub const ILLUMINA_TRUSEQ_R1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
pub const ILLUMINA_TRUSEQ_R2: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
pub const ILLUMINA_SMALL_RNA: &str = "TGGAATTCTCGGGTGCCAAGG";
pub const ILLUMINA_SMALL_DNA_1: &str = "CTGTAGGCACCATCAATCGT";
pub const ILLUMINA_SMALL_DNA_2: &str = "CTGTCTCTTATACACATCT";
pub const ILLUMINA_SMALL_DNA_3: &str = "ATGTGTATAAGAGACA";
pub const NEXTERA_R1: &str = "CTGTCTCTTATACACATCT";
pub const NEXTERA_R2: &str = "CTGTCTCTTATACACATCT";

/// Adapter configuration
#[derive(Debug, Clone)]
pub struct AdapterConfig {
    /// List of adapter sequences to detect
    pub adapters: Vec<String>,

    /// Minimum overlap for adapter detection (default: 8bp)
    pub min_overlap: usize,

    /// Maximum error rate for fuzzy matching (default: 0.1 = 10%)
    pub max_error_rate: f64,

    /// Minimum adapter length after trimming
    pub min_adapter_length: usize,
}

impl Default for AdapterConfig {
    fn default() -> Self {
        Self {
            adapters: vec![
                ILLUMINA_TRUSEQ_R1.to_string(),
                ILLUMINA_TRUSEQ_R2.to_string(),
                ILLUMINA_SMALL_RNA.to_string(),
                NEXTERA_R1.to_string(),
                ILLUMINA_SMALL_DNA_1.to_string(),
                ILLUMINA_SMALL_DNA_2.to_string(),
                ILLUMINA_SMALL_DNA_3.to_string(),
            ],
            min_overlap: 8,
            max_error_rate: 0.1,
            min_adapter_length: 5,
        }
    }
}

/// Adapter trimmer
pub struct AdapterTrimmer {
    config: AdapterConfig,
}

impl AdapterTrimmer {
    pub fn new(config: AdapterConfig) -> Self {
        Self { config }
    }

    /// Detect adapter in sequence and return trim position
    /// Optimized: Only search the last 50bp of the read where adapters typically occur
    pub fn detect_adapter(&self, sequence: &str) -> Option<AdapterMatch> {
        let seq_bytes = sequence.as_bytes();

        // Optimization 1: Only search last 50bp where adapters typically are
        let search_start = seq_bytes.len().saturating_sub(50);

        for adapter in &self.config.adapters {
            let adapter_bytes = adapter.as_bytes();

            // Optimization 2: Use Boyer-Moore-like approach - check largest overlaps first
            // and only in the region where adapters appear (3' end)
            for overlap in (self.config.min_overlap..=adapter_bytes.len().min(seq_bytes.len())).rev() {
                // Only search in the likely adapter region (last 50bp)
                let search_end = seq_bytes.len().saturating_sub(overlap) + 1;

                for read_pos in search_start..search_end {
                    let read_segment = &seq_bytes[read_pos..read_pos + overlap];
                    let adapter_segment = &adapter_bytes[0..overlap];

                    // Quick check: if first byte doesn't match and error rate is 0, skip
                    if self.config.max_error_rate == 0.0 && read_segment[0] != adapter_segment[0] {
                        continue;
                    }

                    let mismatches = Self::count_mismatches(read_segment, adapter_segment);
                    let error_rate = mismatches as f64 / overlap as f64;

                    if error_rate <= self.config.max_error_rate {
                        return Some(AdapterMatch {
                            adapter: adapter.clone(),
                            position: read_pos,
                            length: overlap,
                            mismatches,
                            error_rate,
                        });
                    }
                }
            }
        }

        None
    }

    /// Count mismatches between two sequences
    fn count_mismatches(seq1: &[u8], seq2: &[u8]) -> usize {
        seq1.iter().zip(seq2.iter()).filter(|(a, b)| a != b).count()
    }

    /// Trim adapter from sequence
    pub fn trim_adapter(&self, sequence: &str) -> (String, Option<AdapterMatch>) {
        if let Some(adapter_match) = self.detect_adapter(sequence) {
            let trimmed = sequence[0..adapter_match.position].to_string();
            (trimmed, Some(adapter_match))
        } else {
            (sequence.to_string(), None)
        }
    }

    /// Trim adapter from sequence and quality scores
    pub fn trim_adapter_with_quality(
        &self,
        sequence: &str,
        quality: &[u8],
    ) -> (String, Vec<u8>, Option<AdapterMatch>) {
        if let Some(adapter_match) = self.detect_adapter(sequence) {
            let trimmed_seq = sequence[0..adapter_match.position].to_string();
            let trimmed_qual = quality[0..adapter_match.position].to_vec();
            (trimmed_seq, trimmed_qual, Some(adapter_match))
        } else {
            (sequence.to_string(), quality.to_vec(), None)
        }
    }

    /// Detect all adapters in sequence (for reporting)
    pub fn detect_all_adapters(&self, sequence: &str) -> Vec<AdapterMatch> {
        let mut matches = Vec::new();
        let mut remaining_seq = sequence;
        let mut offset = 0;

        loop {
            if let Some(mut adapter_match) = self.detect_adapter(remaining_seq) {
                adapter_match.position += offset;
                matches.push(adapter_match.clone());

                // Continue searching after this match
                offset = adapter_match.position + adapter_match.length;
                if offset >= sequence.len() {
                    break;
                }
                remaining_seq = &sequence[offset..];
            } else {
                break;
            }
        }

        matches
    }
}

/// Adapter match information
#[derive(Debug, Clone)]
pub struct AdapterMatch {
    pub adapter: String,
    pub position: usize,
    pub length: usize,
    pub mismatches: usize,
    pub error_rate: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_adapter_detection() {
        let config = AdapterConfig::default();
        let trimmer = AdapterTrimmer::new(config);

        // Read with TruSeq adapter at end
        let sequence = format!("ATCGATCGATCGATCG{}", ILLUMINA_TRUSEQ_R1);

        let result = trimmer.detect_adapter(&sequence);
        assert!(result.is_some());

        let adapter_match = result.unwrap();
        assert_eq!(adapter_match.position, 16);
        assert_eq!(adapter_match.mismatches, 0);
    }

    #[test]
    fn test_partial_adapter_detection() {
        let config = AdapterConfig {
            min_overlap: 8,
            ..Default::default()
        };
        let trimmer = AdapterTrimmer::new(config);

        // Read with partial TruSeq adapter (first 10bp)
        let sequence = format!("ATCGATCGATCGATCG{}", &ILLUMINA_TRUSEQ_R1[0..10]);

        let result = trimmer.detect_adapter(&sequence);
        assert!(result.is_some());

        let adapter_match = result.unwrap();
        assert_eq!(adapter_match.position, 16);
        assert_eq!(adapter_match.length, 10);
    }

    #[test]
    fn test_fuzzy_adapter_detection() {
        let config = AdapterConfig {
            max_error_rate: 0.15, // Allow 15% errors
            min_overlap: 10,
            ..Default::default()
        };
        let trimmer = AdapterTrimmer::new(config);

        // TruSeq adapter with 1 mismatch: AGATCGGAAG → AGATCGGGAG (A→G at position 8)
        let adapter_with_error = "AGATCGGGAGAGCACACGTCTGAACTCCAGTCA";
        let sequence = format!("ATCGATCGATCGATCG{}", adapter_with_error);

        let result = trimmer.detect_adapter(&sequence);
        // With 15% error rate and 10bp min overlap, 1 error in 10bp = 10% should match
        assert!(
            result.is_some(),
            "Should detect adapter with 1 mismatch in 10bp (10% error rate)"
        );
    }

    #[test]
    fn test_adapter_trimming() {
        let config = AdapterConfig::default();
        let trimmer = AdapterTrimmer::new(config);

        let original = format!("ATCGATCGATCGATCG{}", ILLUMINA_TRUSEQ_R1);
        let (trimmed, adapter_match) = trimmer.trim_adapter(&original);

        assert_eq!(trimmed, "ATCGATCGATCGATCG");
        assert!(adapter_match.is_some());
    }

    #[test]
    fn test_adapter_trimming_with_quality() {
        let config = AdapterConfig::default();
        let trimmer = AdapterTrimmer::new(config);

        let original = format!("ATCGATCGATCGATCG{}", ILLUMINA_TRUSEQ_R1);
        let quality = vec![b'I'; original.len()]; // Q40

        let (trimmed_seq, trimmed_qual, adapter_match) =
            trimmer.trim_adapter_with_quality(&original, &quality);

        assert_eq!(trimmed_seq, "ATCGATCGATCGATCG");
        assert_eq!(trimmed_qual.len(), 16);
        assert!(adapter_match.is_some());
    }

    #[test]
    fn test_no_adapter() {
        let config = AdapterConfig::default();
        let trimmer = AdapterTrimmer::new(config);

        let sequence = "ATCGATCGATCGATCG"; // No adapter

        let result = trimmer.detect_adapter(sequence);
        assert!(result.is_none());

        let (trimmed, adapter_match) = trimmer.trim_adapter(sequence);
        assert_eq!(trimmed, sequence);
        assert!(adapter_match.is_none());
    }

    #[test]
    fn test_multiple_adapters() {
        let config = AdapterConfig::default();
        let trimmer = AdapterTrimmer::new(config);

        // Sequence with multiple adapters (unusual but possible)
        let sequence = format!("{}ATCG{}", ILLUMINA_TRUSEQ_R1, NEXTERA_R1);

        let matches = trimmer.detect_all_adapters(&sequence);
        assert!(matches.len() >= 1); // Should find at least the first one
    }

    #[test]
    fn test_min_overlap_enforcement() {
        let config = AdapterConfig {
            min_overlap: 15, // Require 15bp overlap
            ..Default::default()
        };
        let trimmer = AdapterTrimmer::new(config);

        // Only 10bp of adapter (below minimum)
        let sequence = format!("ATCGATCGATCGATCG{}", &ILLUMINA_TRUSEQ_R1[0..10]);

        let result = trimmer.detect_adapter(&sequence);
        assert!(result.is_none()); // Should not detect with <15bp overlap
    }
}
