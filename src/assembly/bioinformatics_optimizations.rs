//! Bioinformatics-specific optimizations for assembly graph construction
//! =====================================================================
//! 
//! This module provides specialized optimizations for genomic data processing including:
//! - Bit-packed k-mer representations for memory efficiency
//! - Rolling hash implementations for streaming k-mer processing
//! - SIMD operations for nucleotide sequence processing
//! - Specialized data structures for large-scale genomic datasets

use ahash::{AHashMap, AHashSet};
use anyhow::{Result, anyhow};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Bit-packed k-mer representation for memory-efficient storage
/// Each nucleotide is encoded in 2 bits: A=00, C=01, G=10, T=11
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct BitPackedKmer {
    /// Packed nucleotides (up to 32 nucleotides per u64)
    pub packed_data: Vec<u64>,
    /// Length of the k-mer
    pub k: usize,
    /// Hash value for quick comparison
    pub hash: u64,
    /// Whether this is the canonical representation
    pub is_canonical: bool,
}

impl BitPackedKmer {
    /// Create a new bit-packed k-mer from a DNA sequence
    pub fn new(sequence: &str) -> Result<Self> {
        if sequence.is_empty() {
            return Err(anyhow!("Empty sequence provided"));
        }
        
        let k = sequence.len();
        if k > 1024 {
            return Err(anyhow!("K-mer too long: {} (max 1024)", k));
        }
        
        let sequence_upper = sequence.to_uppercase();
        let rc = Self::reverse_complement(&sequence_upper)?;
        
        let (canonical, is_canonical) = if sequence_upper <= rc {
            (sequence_upper, true)
        } else {
            (rc, false)
        };
        
        let packed_data = Self::pack_sequence(&canonical)?;
        let hash = Self::compute_hash(&packed_data, k);
        
        Ok(Self {
            packed_data,
            k,
            hash,
            is_canonical,
        })
    }
    
    /// Pack a DNA sequence into bit representation
    fn pack_sequence(sequence: &str) -> Result<Vec<u64>> {
        let mut packed = Vec::new();
        let mut current_word = 0u64;
        let mut bits_used = 0;
        
        for nucleotide in sequence.chars() {
            let bits = match nucleotide {
                'A' => 0b00,
                'C' => 0b01,
                'G' => 0b10,
                'T' => 0b11,
                _ => return Err(anyhow!("Invalid nucleotide: {}", nucleotide)),
            };
            
            current_word |= bits << (62 - bits_used);
            bits_used += 2;
            
            if bits_used >= 64 {
                packed.push(current_word);
                current_word = 0;
                bits_used = 0;
            }
        }
        
        if bits_used > 0 {
            packed.push(current_word);
        }
        
        if packed.is_empty() {
            packed.push(0);
        }
        
        Ok(packed)
    }
    
    /// Unpack the bit representation back to a DNA sequence
    pub fn unpack_sequence(&self) -> String {
        let mut sequence = String::with_capacity(self.k);
        let mut remaining = self.k;
        
        for &word in &self.packed_data {
            let nucleotides_in_word = remaining.min(32);
            
            for i in 0..nucleotides_in_word {
                let bits = (word >> (62 - i * 2)) & 0b11;
                let nucleotide = match bits {
                    0b00 => 'A',
                    0b01 => 'C',
                    0b10 => 'G',
                    0b11 => 'T',
                    _ => unreachable!(),
                };
                sequence.push(nucleotide);
            }
            
            remaining -= nucleotides_in_word;
            if remaining == 0 {
                break;
            }
        }
        
        sequence
    }
    
    /// Compute reverse complement of a DNA sequence
    fn reverse_complement(sequence: &str) -> Result<String> {
        sequence
            .chars()
            .rev()
            .map(|c| match c {
                'A' => Ok('T'),
                'T' => Ok('A'),
                'G' => Ok('C'),
                'C' => Ok('G'),
                _ => Err(anyhow!("Invalid nucleotide: {}", c)),
            })
            .collect()
    }
    
    /// Compute hash for packed data
    fn compute_hash(packed_data: &[u64], k: usize) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        let mut hasher = DefaultHasher::new();
        packed_data.hash(&mut hasher);
        k.hash(&mut hasher);
        hasher.finish()
    }
    
    /// Get memory usage in bytes
    pub fn memory_usage(&self) -> usize {
        self.packed_data.len() * std::mem::size_of::<u64>() + std::mem::size_of::<Self>()
    }
}

/// Rolling hash implementation for efficient k-mer processing
#[derive(Debug, Clone)]
pub struct RollingHash {
    /// Current hash value
    hash: u64,
    /// K-mer size
    k: usize,
    /// Base for polynomial hash
    base: u64,
    /// Precomputed power for efficient rolling
    base_power: u64,
    /// Current k-mer window
    window: Vec<u8>,
}

impl RollingHash {
    /// Create a new rolling hash for given k-mer size
    pub fn new(k: usize) -> Self {
        let base = 4u64; // 4 possible nucleotides
        let base_power = base.pow((k as u32).saturating_sub(1));
        
        Self {
            hash: 0,
            k,
            base,
            base_power,
            window: Vec::with_capacity(k),
        }
    }
    
    /// Add a nucleotide and compute rolling hash
    pub fn push(&mut self, nucleotide: char) -> Result<Option<u64>> {
        let encoded = Self::encode_nucleotide(nucleotide)?;
        
        if self.window.len() < self.k {
            // Building initial window
            self.window.push(encoded);
            self.hash = self.hash * self.base + encoded as u64;
            
            if self.window.len() == self.k {
                Ok(Some(self.hash))
            } else {
                Ok(None)
            }
        } else {
            // Rolling update
            let outgoing = self.window[0] as u64;
            self.hash = (self.hash - outgoing * self.base_power) * self.base + encoded as u64;
            
            // Shift window
            self.window.remove(0);
            self.window.push(encoded);
            
            Ok(Some(self.hash))
        }
    }
    
    /// Encode nucleotide to number
    fn encode_nucleotide(nucleotide: char) -> Result<u8> {
        match nucleotide.to_ascii_uppercase() {
            'A' => Ok(0),
            'C' => Ok(1),
            'G' => Ok(2),
            'T' => Ok(3),
            _ => Err(anyhow!("Invalid nucleotide: {}", nucleotide)),
        }
    }
    
    /// Get current k-mer as string
    pub fn current_kmer(&self) -> String {
        if self.window.len() != self.k {
            return String::new();
        }
        
        self.window
            .iter()
            .map(|&encoded| match encoded {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => 'N',
            })
            .collect()
    }
}

/// Streaming k-mer processor for memory-efficient large file processing
pub struct StreamingKmerProcessor {
    k: usize,
    rolling_hash: RollingHash,
    canonical_kmers: AHashSet<u64>,
    kmer_counts: AHashMap<u64, u32>,
    processed_sequences: usize,
    total_kmers: usize,
}

impl StreamingKmerProcessor {
    /// Create new streaming processor
    pub fn new(k: usize) -> Self {
        Self {
            k,
            rolling_hash: RollingHash::new(k),
            canonical_kmers: AHashSet::new(),
            kmer_counts: AHashMap::new(),
            processed_sequences: 0,
            total_kmers: 0,
        }
    }
    
    /// Process a sequence and extract k-mers
    pub fn process_sequence(&mut self, sequence: &str) -> Result<Vec<(String, u64)>> {
        let mut kmers = Vec::new();
        self.rolling_hash = RollingHash::new(self.k); // Reset for new sequence
        
        for nucleotide in sequence.chars() {
            if let Some(_hash) = self.rolling_hash.push(nucleotide)? {
                let kmer_seq = self.rolling_hash.current_kmer();
                
                // Create canonical k-mer
                let bit_kmer = BitPackedKmer::new(&kmer_seq)?;
                let canonical_hash = bit_kmer.hash;
                
                // Update counts
                *self.kmer_counts.entry(canonical_hash).or_insert(0) += 1;
                self.canonical_kmers.insert(canonical_hash);
                
                kmers.push((kmer_seq, canonical_hash));
                self.total_kmers += 1;
            }
        }
        
        self.processed_sequences += 1;
        Ok(kmers)
    }
    
    /// Get statistics about processed data
    pub fn get_statistics(&self) -> ProcessingStatistics {
        ProcessingStatistics {
            unique_kmers: self.canonical_kmers.len(),
            total_kmers: self.total_kmers,
            processed_sequences: self.processed_sequences,
            average_kmer_frequency: if self.canonical_kmers.is_empty() {
                0.0
            } else {
                self.total_kmers as f64 / self.canonical_kmers.len() as f64
            },
            memory_usage_bytes: self.estimate_memory_usage(),
        }
    }
    
    /// Estimate memory usage
    fn estimate_memory_usage(&self) -> usize {
        let kmers_size = self.canonical_kmers.len() * std::mem::size_of::<u64>();
        let counts_size = self.kmer_counts.len() * (std::mem::size_of::<u64>() + std::mem::size_of::<u32>());
        let rolling_hash_size = std::mem::size_of::<RollingHash>();
        
        kmers_size + counts_size + rolling_hash_size
    }
    
    /// Get frequent k-mers above threshold
    pub fn get_frequent_kmers(&self, min_frequency: u32) -> Vec<(u64, u32)> {
        self.kmer_counts
            .iter()
            .filter(|(_, &count)| count >= min_frequency)
            .map(|(&hash, &count)| (hash, count))
            .collect()
    }
}

/// Statistics for k-mer processing
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingStatistics {
    pub unique_kmers: usize,
    pub total_kmers: usize,
    pub processed_sequences: usize,
    pub average_kmer_frequency: f64,
    pub memory_usage_bytes: usize,
}

/// SIMD-optimized nucleotide operations (placeholder for future SIMD implementation)
pub struct SimdNucleotideOps;

impl SimdNucleotideOps {
    /// Count nucleotides in parallel using SIMD when available
    pub fn count_nucleotides(sequence: &str) -> [usize; 4] {
        // Parallel counting using rayon for now
        // Future: implement with SIMD intrinsics for better performance
        let chunks: Vec<_> = sequence.as_bytes().par_chunks(1024).collect();
        
        chunks
            .par_iter()
            .map(|chunk| {
                let mut counts = [0usize; 4];
                for &byte in *chunk {
                    match byte.to_ascii_uppercase() {
                        b'A' => counts[0] += 1,
                        b'C' => counts[1] += 1,
                        b'G' => counts[2] += 1,
                        b'T' => counts[3] += 1,
                        _ => {} // Ignore invalid characters
                    }
                }
                counts
            })
            .reduce(
                || [0usize; 4],
                |mut acc, counts| {
                    for i in 0..4 {
                        acc[i] += counts[i];
                    }
                    acc
                },
            )
    }
    
    /// Calculate GC content efficiently
    pub fn gc_content(sequence: &str) -> f64 {
        let [a_count, c_count, g_count, t_count] = Self::count_nucleotides(sequence);
        let total = a_count + c_count + g_count + t_count;
        
        if total == 0 {
            0.0
        } else {
            (c_count + g_count) as f64 / total as f64
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_packed_kmer() {
        let sequence = "ATCGATCG";
        let kmer = BitPackedKmer::new(sequence).unwrap();
        let unpacked = kmer.unpack_sequence();
        assert_eq!(unpacked, sequence);
    }

    #[test]
    fn test_rolling_hash() {
        let mut rolling_hash = RollingHash::new(3);
        let sequence = "ATCGATCG";
        
        let mut hashes = Vec::new();
        for nucleotide in sequence.chars() {
            if let Some(hash) = rolling_hash.push(nucleotide).unwrap() {
                hashes.push((rolling_hash.current_kmer(), hash));
            }
        }
        
        assert_eq!(hashes.len(), sequence.len() - 3 + 1);
    }

    #[test]
    fn test_streaming_processor() {
        let mut processor = StreamingKmerProcessor::new(4);
        let kmers = processor.process_sequence("ATCGATCGATCG").unwrap();
        
        assert!(!kmers.is_empty());
        let stats = processor.get_statistics();
        assert!(stats.unique_kmers > 0);
        assert!(stats.total_kmers > 0);
    }

    #[test]
    fn test_simd_ops() {
        let sequence = "ATCGATCGATCG";
        let counts = SimdNucleotideOps::count_nucleotides(sequence);
        let gc = SimdNucleotideOps::gc_content(sequence);
        
        assert_eq!(counts[0] + counts[1] + counts[2] + counts[3], sequence.len());
        assert!(gc >= 0.0 && gc <= 1.0);
    }
}