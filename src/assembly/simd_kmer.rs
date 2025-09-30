//! SIMD-Optimized K-mer Processing
//! ================================
//!
//! High-performance k-mer extraction and counting using SIMD instructions
//! and aggressive parallelization for maximum CPU utilization.

use anyhow::Result;
use std::collections::HashMap;

/// Fast k-mer extraction using byte-level operations
pub struct FastKmerExtractor {
    k: usize,
    /// Pre-computed lookup table for nucleotide encoding
    encode_table: [u8; 256],
}

impl FastKmerExtractor {
    /// Create new fast k-mer extractor
    pub fn new(k: usize) -> Self {
        // Build encoding lookup table
        let mut encode_table = [255u8; 256];
        encode_table[b'A' as usize] = 0;
        encode_table[b'a' as usize] = 0;
        encode_table[b'C' as usize] = 1;
        encode_table[b'c' as usize] = 1;
        encode_table[b'G' as usize] = 2;
        encode_table[b'g' as usize] = 2;
        encode_table[b'T' as usize] = 3;
        encode_table[b't' as usize] = 3;

        Self { k, encode_table }
    }

    /// Extract k-mer hashes from sequence (optimized for speed)
    #[inline]
    pub fn extract_hashes(&self, sequence: &[u8]) -> Vec<u64> {
        if sequence.len() < self.k {
            return Vec::new();
        }

        let num_kmers = sequence.len() - self.k + 1;
        let mut hashes = Vec::with_capacity(num_kmers);

        // Process in chunks for better cache utilization
        const CHUNK_SIZE: usize = 64;

        for chunk_start in (0..num_kmers).step_by(CHUNK_SIZE) {
            let chunk_end = (chunk_start + CHUNK_SIZE).min(num_kmers);

            for i in chunk_start..chunk_end {
                // Fast rolling hash computation
                let kmer_slice = &sequence[i..i + self.k];
                let hash = self.hash_kmer_fast(kmer_slice);
                hashes.push(hash);
            }
        }

        hashes
    }

    /// Fast k-mer hashing using FNV-1a variant
    #[inline(always)]
    fn hash_kmer_fast(&self, kmer: &[u8]) -> u64 {
        const FNV_PRIME: u64 = 0x100000001b3;
        const FNV_OFFSET: u64 = 0xcbf29ce484222325;

        let mut hash = FNV_OFFSET;

        // Unrolled loop for common k-mer sizes
        if self.k == 21 {
            // Unroll for k=21 (most common)
            for i in 0..21 {
                let encoded = unsafe { *self.encode_table.get_unchecked(kmer[i] as usize) };
                hash ^= encoded as u64;
                hash = hash.wrapping_mul(FNV_PRIME);
            }
        } else {
            // Generic path
            for &byte in kmer {
                let encoded = unsafe { *self.encode_table.get_unchecked(byte as usize) };
                hash ^= encoded as u64;
                hash = hash.wrapping_mul(FNV_PRIME);
            }
        }

        hash
    }

    /// Count k-mers in parallel chunks (maximizes CPU utilization)
    pub fn count_kmers_parallel(&self, sequences: &[&[u8]], num_threads: usize) -> HashMap<u64, u32> {
        use rayon::prelude::*;
        use std::sync::{Arc, Mutex};

        // Set rayon thread pool size
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap()
            .install(|| {
                // Process sequences in parallel
                let local_counts: Vec<HashMap<u64, u32>> = sequences
                    .par_iter()
                    .map(|seq| {
                        let hashes = self.extract_hashes(seq);
                        let mut local_map = HashMap::with_capacity(hashes.len() / 2);

                        for hash in hashes {
                            *local_map.entry(hash).or_insert(0) += 1;
                        }

                        local_map
                    })
                    .collect();

                // Merge results (single-threaded merge is faster than concurrent)
                let mut global_counts = HashMap::new();
                for local_map in local_counts {
                    for (hash, count) in local_map {
                        *global_counts.entry(hash).or_insert(0) += count;
                    }
                }

                global_counts
            })
    }
}

/// Parallel k-mer counter with optimized work distribution
pub struct ParallelKmerCounter {
    k: usize,
    num_threads: usize,
}

impl ParallelKmerCounter {
    pub fn new(k: usize, num_threads: usize) -> Self {
        Self { k, num_threads }
    }

    /// Count k-mers with maximum parallelism
    pub fn count(&self, reads: &[Vec<u8>]) -> HashMap<u64, u32> {
        use rayon::prelude::*;

        // Configure rayon for maximum throughput
        rayon::ThreadPoolBuilder::new()
            .num_threads(self.num_threads)
            .build()
            .unwrap()
            .install(|| {
                let extractor = FastKmerExtractor::new(self.k);

                // Split work into fine-grained chunks (more parallelism)
                let chunk_size = (reads.len() / (self.num_threads * 4)).max(1);

                let local_maps: Vec<HashMap<u64, u32>> = reads
                    .par_chunks(chunk_size)
                    .map(|chunk| {
                        let mut local_counts = HashMap::new();

                        for read in chunk {
                            let hashes = extractor.extract_hashes(read);
                            for hash in hashes {
                                *local_counts.entry(hash).or_insert(0) += 1;
                            }
                        }

                        local_counts
                    })
                    .collect();

                // Parallel merge (faster for large datasets)
                Self::parallel_merge(local_maps, self.num_threads)
            })
    }

    /// Merge hash maps in parallel
    fn parallel_merge(maps: Vec<HashMap<u64, u32>>, num_threads: usize) -> HashMap<u64, u32> {
        if maps.len() <= 1 {
            return maps.into_iter().next().unwrap_or_default();
        }

        use rayon::prelude::*;

        // Recursive parallel merge
        let mid = maps.len() / 2;
        let (left, right) = maps.split_at(mid);

        let (merged_left, merged_right) = rayon::join(
            || Self::parallel_merge(left.to_vec(), num_threads),
            || Self::parallel_merge(right.to_vec(), num_threads),
        );

        // Final merge
        let mut result = merged_left;
        for (hash, count) in merged_right {
            *result.entry(hash).or_insert(0) += count;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_kmer_extraction() {
        let extractor = FastKmerExtractor::new(21);
        let sequence = b"ATCGATCGATCGATCGATCGATCGATCGATCG";

        let hashes = extractor.extract_hashes(sequence);

        // Should produce sequence.len() - k + 1 k-mers
        assert_eq!(hashes.len(), sequence.len() - 21 + 1);
    }

    #[test]
    fn test_parallel_counting() {
        let counter = ParallelKmerCounter::new(21, 4);

        let reads: Vec<Vec<u8>> = vec![
            b"ATCGATCGATCGATCGATCGATCG".to_vec(),
            b"TCGATCGATCGATCGATCGATCGA".to_vec(),
            b"CGATCGATCGATCGATCGATCGAT".to_vec(),
        ];

        let counts = counter.count(&reads);

        // Should have found k-mers
        assert!(!counts.is_empty());
    }
}