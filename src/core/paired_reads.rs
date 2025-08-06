use anyhow::{Result, anyhow, Context};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use rayon::prelude::*;
use ahash::AHashMap;

use crate::core::data_structures::{CorrectedRead, BaseCorrection, ReadPosition, Strand};

/// Enhanced paired-end read support for metagenomic assembly
/// 
/// This module provides comprehensive paired-end read handling including:
/// - Paired read data structures with mate relationship tracking
/// - Support for both interleaved and separate R1/R2 FASTQ files
/// - Insert size estimation and validation
/// - Mate pair relationship validation
/// - Efficient memory layout for large-scale processing

/// A pair of reads (R1 and R2) with their relationship metadata
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadPair {
    /// Forward read (R1)
    pub forward: PairedRead,
    /// Reverse read (R2)
    pub reverse: PairedRead,
    /// Pair-specific metadata
    pub pair_info: PairInfo,
}

/// Individual read that is part of a paired-end set
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedRead {
    /// Unique read identifier
    pub id: usize,
    /// Pair identifier (shared between R1 and R2)
    pub pair_id: String,
    /// Read orientation (Forward/Reverse)
    pub orientation: ReadOrientation,
    /// Original sequence from FASTQ
    pub original: String,
    /// Error-corrected sequence
    pub corrected: String,
    /// Quality scores
    pub quality_scores: Vec<u8>,
    /// Error corrections applied
    pub corrections: Vec<BaseCorrection>,
    /// Read-specific metadata
    pub read_info: ReadInfo,
}

/// Metadata specific to the read pair
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairInfo {
    /// Estimated insert size (distance between reads)
    pub insert_size: Option<usize>,
    /// Insert size confidence
    pub insert_confidence: f64,
    /// Whether reads properly align as a pair
    pub properly_paired: bool,
    /// Library type information
    pub library_type: LibraryType,
    /// Fragment length information
    pub fragment_length: Option<usize>,
}

/// Metadata specific to individual reads
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadInfo {
    /// Read length
    pub length: usize,
    /// Average quality score
    pub avg_quality: f64,
    /// GC content
    pub gc_content: f64,
    /// Complexity score (entropy-based)
    pub complexity: f64,
    /// Whether read passes quality filters
    pub passes_filter: bool,
}

/// Read orientation in the pair
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReadOrientation {
    Forward,  // R1
    Reverse,  // R2
}

/// Library preparation type
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum LibraryType {
    /// Paired-end with forward/reverse orientation
    PairedEnd,
    /// Mate-pair with reverse/forward orientation  
    MatePair,
    /// Single-end (for compatibility)
    SingleEnd,
    /// Unknown library type
    Unknown,
}

/// Statistics for paired-end read processing
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PairedReadStats {
    /// Total number of read pairs processed
    pub total_pairs: usize,
    /// Number of properly paired reads
    pub properly_paired: usize,
    /// Number of singleton reads (missing mate)
    pub singletons: usize,
    /// Insert size statistics
    pub insert_size_stats: InsertSizeStats,
    /// Quality statistics
    pub quality_stats: QualityStats,
    /// Error correction statistics
    pub correction_stats: CorrectionStats,
}

/// Insert size distribution statistics
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct InsertSizeStats {
    pub mean: f64,
    pub median: f64,
    pub std_dev: f64,
    pub min: usize,
    pub max: usize,
    pub distribution: HashMap<usize, usize>,
}

/// Quality score statistics
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct QualityStats {
    pub mean_quality_r1: f64,
    pub mean_quality_r2: f64,
    pub quality_distribution: HashMap<u8, usize>,
    pub low_quality_pairs: usize,
}

/// Error correction statistics
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CorrectionStats {
    pub corrections_applied: usize,
    pub bases_corrected_r1: usize,
    pub bases_corrected_r2: usize,
    pub correction_types: HashMap<String, usize>,
}

impl ReadPair {
    /// Create a new read pair from individual reads
    pub fn new(forward: PairedRead, reverse: PairedRead) -> Result<Self> {
        // Validate that reads are properly paired
        if forward.pair_id != reverse.pair_id {
            return Err(anyhow!(
                "Pair ID mismatch: {} vs {}",
                forward.pair_id,
                reverse.pair_id
            ));
        }

        if forward.orientation == reverse.orientation {
            return Err(anyhow!(
                "Invalid pair: both reads have the same orientation {:?}",
                forward.orientation
            ));
        }

        // Calculate initial pair information
        let pair_info = PairInfo {
            insert_size: None, // Will be estimated later
            insert_confidence: 0.0,
            properly_paired: true,
            library_type: LibraryType::PairedEnd,
            fragment_length: None,
        };

        Ok(Self {
            forward,
            reverse,
            pair_info,
        })
    }

    /// Estimate insert size between paired reads
    pub fn estimate_insert_size(&mut self) -> Result<()> {
        // Simple estimation based on read lengths and expected fragment size
        // In a real implementation, this would use alignment information
        let forward_len = self.forward.read_info.length;
        let reverse_len = self.reverse.read_info.length;
        
        // Estimate based on typical insert sizes for metagenomics
        let estimated_insert = (forward_len + reverse_len + 200).max(300);
        
        self.pair_info.insert_size = Some(estimated_insert);
        self.pair_info.insert_confidence = 0.5; // Placeholder confidence
        
        Ok(())
    }

    /// Validate that the read pair is properly formed
    pub fn validate_pair(&self) -> Result<bool> {
        // Check orientation
        let forward_is_forward = self.forward.orientation == ReadOrientation::Forward;
        let reverse_is_reverse = self.reverse.orientation == ReadOrientation::Reverse;
        
        if !forward_is_forward || !reverse_is_reverse {
            return Ok(false);
        }

        // Check pair IDs match
        if self.forward.pair_id != self.reverse.pair_id {
            return Ok(false);
        }

        // Check read lengths are reasonable
        if self.forward.read_info.length < 50 || self.reverse.read_info.length < 50 {
            return Ok(false);
        }

        // Check quality
        if self.forward.read_info.avg_quality < 20.0 || self.reverse.read_info.avg_quality < 20.0 {
            return Ok(false);
        }

        Ok(true)
    }

    /// Get the combined sequence representation for assembly
    pub fn get_combined_sequence(&self) -> String {
        // For assembly, we might want to represent the pair as a combined sequence
        // with a gap representing the insert
        let gap_size = self.pair_info.insert_size.unwrap_or(300);
        let gap = "N".repeat(gap_size.saturating_sub(self.forward.read_info.length + self.reverse.read_info.length));
        
        format!("{}{}{}", self.forward.corrected, gap, reverse_complement(&self.reverse.corrected))
    }

    /// Extract k-mers considering paired-end information
    pub fn extract_paired_kmers(&self, k: usize) -> Result<Vec<PairedKmer>> {
        let mut kmers = Vec::new();
        
        // Extract k-mers from forward read
        for (i, window) in self.forward.corrected.as_bytes().windows(k).enumerate() {
            let kmer_str = std::str::from_utf8(window)?;
            // Skip k-mers containing ambiguous bases
            if kmer_str.chars().any(|c| matches!(c, 'N' | 'n')) {
                continue;
            }
            kmers.push(PairedKmer {
                sequence: kmer_str.to_string(),
                pair_id: self.forward.pair_id.clone(),
                read_orientation: ReadOrientation::Forward,
                position: i,
                mate_distance: self.pair_info.insert_size,
            });
        }

        // Extract k-mers from reverse read
        for (i, window) in self.reverse.corrected.as_bytes().windows(k).enumerate() {
            let kmer_str = std::str::from_utf8(window)?;
            // Skip k-mers containing ambiguous bases
            if kmer_str.chars().any(|c| matches!(c, 'N' | 'n')) {
                continue;
            }
            kmers.push(PairedKmer {
                sequence: kmer_str.to_string(),
                pair_id: self.reverse.pair_id.clone(),
                read_orientation: ReadOrientation::Reverse,
                position: i,
                mate_distance: self.pair_info.insert_size,
            });
        }

        Ok(kmers)
    }
}

/// K-mer with paired-end context information
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedKmer {
    pub sequence: String,
    pub pair_id: String,
    pub read_orientation: ReadOrientation,
    pub position: usize,
    pub mate_distance: Option<usize>,
}

impl PairedRead {
    /// Create a new paired read
    pub fn new(
        id: usize,
        pair_id: String,
        orientation: ReadOrientation,
        sequence: String,
        quality_scores: Vec<u8>,
    ) -> Self {
        let read_info = ReadInfo {
            length: sequence.len(),
            avg_quality: quality_scores.iter().map(|&q| q as f64).sum::<f64>() / quality_scores.len() as f64,
            gc_content: calculate_gc_content(&sequence),
            complexity: calculate_sequence_complexity(&sequence),
            passes_filter: true, // Initial assumption
        };

        Self {
            id,
            pair_id,
            orientation,
            original: sequence.clone(),
            corrected: sequence,
            quality_scores,
            corrections: Vec::new(),
            read_info,
        }
    }

    /// Apply quality filtering
    pub fn apply_quality_filter(&mut self, min_quality: f64, min_length: usize) {
        self.read_info.passes_filter = 
            self.read_info.avg_quality >= min_quality && 
            self.read_info.length >= min_length;
    }

    /// Convert to regular CorrectedRead for backward compatibility
    pub fn to_corrected_read(&self) -> CorrectedRead {
        CorrectedRead {
            id: self.id,
            original: self.original.clone(),
            corrected: self.corrected.clone(),
            corrections: self.corrections.clone(),
            quality_scores: self.quality_scores.clone(),
        }
    }
}

/// Collection of read pairs with efficient access patterns
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedReadCollection {
    /// All read pairs
    pub pairs: Vec<ReadPair>,
    /// Index from pair_id to pairs vector index
    pub pair_index: AHashMap<String, usize>,
    /// Singleton reads (missing mate)
    pub singletons: Vec<PairedRead>,
    /// Collection statistics
    pub stats: PairedReadStats,
}

impl PairedReadCollection {
    /// Create a new empty collection
    pub fn new() -> Self {
        Self {
            pairs: Vec::new(),
            pair_index: AHashMap::new(),
            singletons: Vec::new(),
            stats: PairedReadStats::default(),
        }
    }

    /// Add a read pair to the collection
    pub fn add_pair(&mut self, pair: ReadPair) -> Result<()> {
        let pair_id = pair.forward.pair_id.clone();
        
        if self.pair_index.contains_key(&pair_id) {
            return Err(anyhow!("Duplicate pair ID: {}", pair_id));
        }

        let index = self.pairs.len();
        self.pair_index.insert(pair_id, index);
        self.pairs.push(pair);
        
        Ok(())
    }

    /// Add a singleton read
    pub fn add_singleton(&mut self, read: PairedRead) {
        self.singletons.push(read);
    }

    /// Get a read pair by pair ID
    pub fn get_pair(&self, pair_id: &str) -> Option<&ReadPair> {
        self.pair_index
            .get(pair_id)
            .and_then(|&index| self.pairs.get(index))
    }

    /// Get mutable reference to a read pair
    pub fn get_pair_mut(&mut self, pair_id: &str) -> Option<&mut ReadPair> {
        if let Some(&index) = self.pair_index.get(pair_id) {
            self.pairs.get_mut(index)
        } else {
            None
        }
    }

    /// Calculate comprehensive statistics
    pub fn calculate_stats(&mut self) {
        self.stats.total_pairs = self.pairs.len();
        self.stats.singletons = self.singletons.len();

        // Calculate properly paired statistics
        self.stats.properly_paired = self.pairs
            .iter()
            .filter(|pair| pair.validate_pair().unwrap_or(false))
            .count();

        // Calculate insert size statistics
        self.calculate_insert_size_stats();
        
        // Calculate quality statistics
        self.calculate_quality_stats();

        // Calculate correction statistics
        self.calculate_correction_stats();
    }

    /// Calculate insert size distribution
    fn calculate_insert_size_stats(&mut self) {
        let insert_sizes: Vec<usize> = self.pairs
            .iter()
            .filter_map(|pair| pair.pair_info.insert_size)
            .collect();

        if insert_sizes.is_empty() {
            return;
        }

        let sum: usize = insert_sizes.iter().sum();
        let mean = sum as f64 / insert_sizes.len() as f64;

        let mut sorted_sizes = insert_sizes.clone();
        sorted_sizes.sort_unstable();
        let median = sorted_sizes[sorted_sizes.len() / 2] as f64;

        let variance = insert_sizes
            .iter()
            .map(|&size| {
                let diff = size as f64 - mean;
                diff * diff
            })
            .sum::<f64>() / insert_sizes.len() as f64;
        
        let std_dev = variance.sqrt();

        // Build distribution histogram
        let mut distribution = HashMap::new();
        for &size in &insert_sizes {
            *distribution.entry(size).or_insert(0) += 1;
        }

        self.stats.insert_size_stats = InsertSizeStats {
            mean,
            median,
            std_dev,
            min: *sorted_sizes.first().unwrap_or(&0),
            max: *sorted_sizes.last().unwrap_or(&0),
            distribution,
        };
    }

    /// Calculate quality score statistics
    fn calculate_quality_stats(&mut self) {
        if self.pairs.is_empty() {
            return;
        }

        let mut r1_qualities = Vec::new();
        let mut r2_qualities = Vec::new();
        let mut quality_dist = HashMap::new();
        let mut low_quality_pairs = 0;

        for pair in &self.pairs {
            let r1_avg = pair.forward.read_info.avg_quality;
            let r2_avg = pair.reverse.read_info.avg_quality;
            
            r1_qualities.push(r1_avg);
            r2_qualities.push(r2_avg);

            // Build quality distribution
            for &qual in &pair.forward.quality_scores {
                *quality_dist.entry(qual).or_insert(0) += 1;
            }
            for &qual in &pair.reverse.quality_scores {
                *quality_dist.entry(qual).or_insert(0) += 1;
            }

            // Count low quality pairs
            if r1_avg < 25.0 || r2_avg < 25.0 {
                low_quality_pairs += 1;
            }
        }

        self.stats.quality_stats = QualityStats {
            mean_quality_r1: r1_qualities.iter().sum::<f64>() / r1_qualities.len() as f64,
            mean_quality_r2: r2_qualities.iter().sum::<f64>() / r2_qualities.len() as f64,
            quality_distribution: quality_dist,
            low_quality_pairs,
        };
    }

    /// Calculate error correction statistics
    fn calculate_correction_stats(&mut self) {
        let mut corrections_applied = 0;
        let mut bases_corrected_r1 = 0;
        let mut bases_corrected_r2 = 0;
        let mut correction_types = HashMap::new();

        for pair in &self.pairs {
            corrections_applied += pair.forward.corrections.len();
            corrections_applied += pair.reverse.corrections.len();

            bases_corrected_r1 += pair.forward.corrections.len();
            bases_corrected_r2 += pair.reverse.corrections.len();

            // Count correction types
            for correction in &pair.forward.corrections {
                let correction_type = format!("{:?}", correction.correction_type);
                *correction_types.entry(correction_type).or_insert(0) += 1;
            }
            for correction in &pair.reverse.corrections {
                let correction_type = format!("{:?}", correction.correction_type);
                *correction_types.entry(correction_type).or_insert(0) += 1;
            }
        }

        self.stats.correction_stats = CorrectionStats {
            corrections_applied,
            bases_corrected_r1,
            bases_corrected_r2,
            correction_types,
        };
    }

    /// Filter collection by quality criteria
    pub fn filter_by_quality(&mut self, min_quality: f64, min_length: usize) {
        self.pairs.retain(|pair| {
            pair.forward.read_info.avg_quality >= min_quality &&
            pair.reverse.read_info.avg_quality >= min_quality &&
            pair.forward.read_info.length >= min_length &&
            pair.reverse.read_info.length >= min_length
        });

        self.singletons.retain(|read| {
            read.read_info.avg_quality >= min_quality &&
            read.read_info.length >= min_length
        });

        // Rebuild index
        self.pair_index.clear();
        for (i, pair) in self.pairs.iter().enumerate() {
            self.pair_index.insert(pair.forward.pair_id.clone(), i);
        }
    }

    /// Extract all k-mers with paired-end context
    pub fn extract_all_paired_kmers(&self, k: usize) -> Result<Vec<PairedKmer>> {
        let mut all_kmers = Vec::new();
        
        // Process pairs in parallel for efficiency
        let pair_kmers: Result<Vec<Vec<PairedKmer>>> = self.pairs
            .par_iter()
            .map(|pair| pair.extract_paired_kmers(k))
            .collect();

        for kmers in pair_kmers? {
            all_kmers.extend(kmers);
        }

        Ok(all_kmers)
    }
}

impl Default for PairedReadCollection {
    fn default() -> Self {
        Self::new()
    }
}

// Utility functions

/// Calculate GC content of a sequence
fn calculate_gc_content(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }

    let gc_count = sequence
        .chars()
        .filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c')
        .count();

    gc_count as f64 / sequence.len() as f64
}

/// Calculate sequence complexity using Shannon entropy
fn calculate_sequence_complexity(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }

    let mut counts = [0u32; 4]; // A, C, G, T
    for byte in sequence.bytes() {
        match byte {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {}
        }
    }

    let total = counts.iter().sum::<u32>() as f64;
    if total == 0.0 {
        return 0.0;
    }

    let entropy = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / total;
            -p * p.log2()
        })
        .sum::<f64>();

    // Normalize entropy (max entropy for DNA is 2.0)
    entropy / 2.0
}

/// Generate reverse complement of a DNA sequence
fn reverse_complement(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            'N' | 'n' => 'N',
            _ => c,
        })
        .collect()
}

/// Extract pair ID from FASTQ header (handles common formats)
pub fn extract_pair_id(header: &str) -> Result<(String, ReadOrientation)> {
    // Remove the '@' prefix if present
    let header = header.strip_prefix('@').unwrap_or(header);
    
    // Handle Illumina format: @SRR123456.1 1:N:0:ACGT or @SRR123456.1/1
    if let Some(space_pos) = header.find(' ') {
        let base_id = header[..space_pos].to_string();
        let suffix = &header[space_pos + 1..];
        
        // Check for /1 or /2 format
        if suffix.starts_with('1') {
            return Ok((base_id, ReadOrientation::Forward));
        } else if suffix.starts_with('2') {
            return Ok((base_id, ReadOrientation::Reverse));
        }
    }
    
    // Handle /1 and /2 suffix format
    if header.ends_with("/1") {
        let base_id = header[..header.len() - 2].to_string();
        return Ok((base_id, ReadOrientation::Forward));
    } else if header.ends_with("/2") {
        let base_id = header[..header.len() - 2].to_string();
        return Ok((base_id, ReadOrientation::Reverse));
    }
    
    // Handle .1 and .2 suffix format
    if header.ends_with(".1") {
        let base_id = header[..header.len() - 2].to_string();
        return Ok((base_id, ReadOrientation::Forward));
    } else if header.ends_with(".2") {
        let base_id = header[..header.len() - 2].to_string();
        return Ok((base_id, ReadOrientation::Reverse));
    }

    // Default: treat as single-end or unknown paired format
    Ok((header.to_string(), ReadOrientation::Forward))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_pair_creation() {
        let forward = PairedRead::new(
            0,
            "pair_1".to_string(),
            ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![30; 12],
        );

        let reverse = PairedRead::new(
            1,
            "pair_1".to_string(),
            ReadOrientation::Reverse,
            "CGATATCGATCG".to_string(),
            vec![30; 12],
        );

        let pair = ReadPair::new(forward, reverse).unwrap();
        assert_eq!(pair.forward.pair_id, pair.reverse.pair_id);
        assert!(pair.validate_pair().unwrap());
    }

    #[test]
    fn test_pair_id_extraction() {
        // Test various FASTQ header formats
        assert_eq!(
            extract_pair_id("@SRR123456.1/1").unwrap(),
            ("SRR123456.1".to_string(), ReadOrientation::Forward)
        );

        assert_eq!(
            extract_pair_id("@SRR123456.2/2").unwrap(),
            ("SRR123456.2".to_string(), ReadOrientation::Reverse)
        );

        assert_eq!(
            extract_pair_id("SRR123456.1 1:N:0:ACGT").unwrap(),
            ("SRR123456.1".to_string(), ReadOrientation::Forward)
        );
    }

    #[test]
    fn test_paired_read_collection() {
        let mut collection = PairedReadCollection::new();

        let forward = PairedRead::new(
            0,
            "pair_1".to_string(),
            ReadOrientation::Forward,
            "ATCGATCGATCGATCGATCG".to_string(),
            vec![30; 20],
        );

        let reverse = PairedRead::new(
            1,
            "pair_1".to_string(),
            ReadOrientation::Reverse,
            "CGATATCGATCGATCGATCG".to_string(),
            vec![30; 20],
        );

        let pair = ReadPair::new(forward, reverse).unwrap();
        collection.add_pair(pair).unwrap();

        assert_eq!(collection.pairs.len(), 1);
        assert!(collection.get_pair("pair_1").is_some());
        
        collection.calculate_stats();
        assert_eq!(collection.stats.total_pairs, 1);
    }

    #[test]
    fn test_insert_size_estimation() {
        let forward = PairedRead::new(
            0,
            "pair_1".to_string(),
            ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![30; 12],
        );

        let reverse = PairedRead::new(
            1,
            "pair_1".to_string(),
            ReadOrientation::Reverse,
            "CGATATCGATCG".to_string(),
            vec![30; 12],
        );

        let mut pair = ReadPair::new(forward, reverse).unwrap();
        pair.estimate_insert_size().unwrap();
        
        assert!(pair.pair_info.insert_size.is_some());
        assert!(pair.pair_info.insert_size.unwrap() > 0);
    }

    #[test]
    fn test_paired_kmer_extraction() {
        let forward = PairedRead::new(
            0,
            "pair_1".to_string(),
            ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![30; 12],
        );

        let reverse = PairedRead::new(
            1,
            "pair_1".to_string(),
            ReadOrientation::Reverse,
            "CGATATCGATCG".to_string(),
            vec![30; 12],
        );

        let pair = ReadPair::new(forward, reverse).unwrap();
        let kmers = pair.extract_paired_kmers(4).unwrap();
        
        assert!(!kmers.is_empty());
        assert!(kmers.iter().any(|k| k.read_orientation == ReadOrientation::Forward));
        assert!(kmers.iter().any(|k| k.read_orientation == ReadOrientation::Reverse));
    }

    #[test]
    fn test_quality_filtering() {
        let mut collection = PairedReadCollection::new();

        // Add high quality pair
        let forward_hq = PairedRead::new(
            0,
            "pair_hq".to_string(),
            ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![35; 12],
        );

        let reverse_hq = PairedRead::new(
            1,
            "pair_hq".to_string(),
            ReadOrientation::Reverse,
            "CGATATCGATCG".to_string(),
            vec![35; 12],
        );

        // Add low quality pair
        let forward_lq = PairedRead::new(
            2,
            "pair_lq".to_string(),
            ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![15; 12],
        );

        let reverse_lq = PairedRead::new(
            3,
            "pair_lq".to_string(),
            ReadOrientation::Reverse,
            "CGATATCGATCG".to_string(),
            vec![15; 12],
        );

        collection.add_pair(ReadPair::new(forward_hq, reverse_hq).unwrap()).unwrap();
        collection.add_pair(ReadPair::new(forward_lq, reverse_lq).unwrap()).unwrap();

        assert_eq!(collection.pairs.len(), 2);

        collection.filter_by_quality(25.0, 10);
        assert_eq!(collection.pairs.len(), 1);
        assert!(collection.get_pair("pair_hq").is_some());
        assert!(collection.get_pair("pair_lq").is_none());
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("atcg"), "cgat");
        assert_eq!(reverse_complement("AAAATTTTGGGGCCCC"), "GGGGCCCCAAAATTTT");
        assert_eq!(reverse_complement("ATCGN"), "NCGAT");
    }

    #[test]
    fn test_gc_content_calculation() {
        assert_eq!(calculate_gc_content("ATCG"), 0.5);
        assert_eq!(calculate_gc_content("AAAA"), 0.0);
        assert_eq!(calculate_gc_content("CCGG"), 1.0);
        assert_eq!(calculate_gc_content(""), 0.0);
    }

    #[test]
    fn test_complexity_calculation() {
        // Uniform sequence (low complexity)
        let complexity1 = calculate_sequence_complexity("AAAAAAA");
        assert!(complexity1 < 0.1);

        // Random sequence (high complexity)
        let complexity2 = calculate_sequence_complexity("ATCGATCG");
        assert!(complexity2 > 0.8);
    }
}