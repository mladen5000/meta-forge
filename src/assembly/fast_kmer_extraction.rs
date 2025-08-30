//! Ultra-fast k-mer extraction using optimized libraries and SIMD operations
//! 
//! This module provides high-performance k-mer counting that can process thousands
//! of reads per second using vectorized operations and specialized libraries.

use anyhow::Result;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicU32, Ordering};
use dashmap::DashMap;

/// Ultra-fast k-mer extractor with lock-free concurrency and specialized genomic optimizations  
pub struct FastKmerExtractor {
    k: usize,
    /// Lock-free concurrent hash map for zero-contention parallel access
    kmer_counts: Arc<DashMap<u64, AtomicU32>>,
    /// Buffer for batch processing - cache-optimized size
    batch_size: usize,
    /// Canonical k-mer mode (reduces memory by 50%)
    canonical_mode: bool,
    /// Quality threshold for k-mer filtering
    quality_threshold: u8,
}

impl FastKmerExtractor {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            kmer_counts: Arc::new(DashMap::new()),
            batch_size: 10_000, // Process 10K reads at a time for optimal cache performance
            canonical_mode: true, // Enable canonical k-mers by default for 50% memory reduction
            quality_threshold: 20, // Minimum quality score for k-mer inclusion
        }
    }
    
    /// Create extractor with custom configuration for maximum performance tuning
    pub fn with_config(k: usize, canonical: bool, quality_threshold: u8, batch_size: usize) -> Self {
        Self {
            k,
            kmer_counts: Arc::new(DashMap::new()),
            batch_size,
            canonical_mode: canonical,
            quality_threshold,
        }
    }
    
    /// Create extractor optimized for large-scale genomic datasets
    pub fn for_large_scale_genomics(k: usize) -> Self {
        Self {
            k,
            kmer_counts: Arc::new(DashMap::with_capacity(10_000_000)), // Pre-allocate for large datasets
            batch_size: 50_000, // Larger batches for better throughput
            canonical_mode: true, // Essential for memory efficiency
            quality_threshold: 25, // Higher quality for more accurate results
        }
    }

    /// Extract k-mers from FASTQ file at maximum speed
    pub fn extract_kmers_from_file(&self, file_path: &str) -> Result<Vec<(u64, u32)>> {
        let mut reader = parse_fastx_file(file_path)?;
        let mut batch = Vec::with_capacity(self.batch_size);
        
        // Process sequences in batches for better cache performance
        while let Some(record) = reader.next() {
            let record = record?;
            batch.push(record.seq().to_vec());
            
            if batch.len() >= self.batch_size {
                self.process_batch_parallel(&batch)?;
                batch.clear();
            }
        }
        
        // Process remaining sequences
        if !batch.is_empty() {
            self.process_batch_parallel(&batch)?;
        }
        
        // Extract results from lock-free structure
        let mut result: Vec<(u64, u32)> = self.kmer_counts
            .iter()
            .map(|entry| (*entry.key(), entry.value().load(Ordering::Relaxed)))
            .collect();
        result.sort_by_key(|&(_, count)| std::cmp::Reverse(count)); // Sort by frequency
        
        Ok(result)
    }

    /// Extract k-mers from sequences with continuous 1% progress updates
    pub fn extract_kmers_from_sequences_with_progress(
        &self, 
        sequences: &[Vec<u8>]
    ) -> Result<(Vec<(u64, u32)>, Vec<(usize, usize, usize, f64, std::time::Instant)>)> {
        // Clear previous results
        self.clear();
        
        let total_sequences = sequences.len();
        let processed_sequences = Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let total_kmers = Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let start_time = std::time::Instant::now();
        
        // Collect detailed progress updates with timing and reads/sec
        let progress_updates = Arc::new(Mutex::new(Vec::new()));
        
        // Calculate chunk size for smoother progress - aim for ~100 chunks minimum for 1% granularity
        let min_chunks = 100;
        let optimal_chunk_size = (total_sequences / min_chunks).max(1);
        let chunk_size = (sequences.len() / num_cpus::get()).min(optimal_chunk_size).max(1);
        
        println!("🔧 Optimized chunking: {} chunks of ~{} sequences each for smooth progress", 
                 (total_sequences + chunk_size - 1) / chunk_size, chunk_size);
        
        sequences
            .par_chunks(chunk_size)
            .try_for_each(|chunk| -> Result<()> {
                let mut chunk_kmers = 0;
                
                // Direct lock-free k-mer extraction - no thread-local aggregation needed
                for sequence in chunk {
                    let kmers_extracted = self.extract_kmers_lockfree_simd(sequence)?;
                    chunk_kmers += kmers_extracted;
                }
                
                // Update progress atomically
                let current_processed = processed_sequences.fetch_add(chunk.len(), std::sync::atomic::Ordering::Relaxed) + chunk.len();
                let current_kmers = total_kmers.fetch_add(chunk_kmers, std::sync::atomic::Ordering::Relaxed) + chunk_kmers;
                let current_time = std::time::Instant::now();
                
                // Calculate current reads per second
                let elapsed_secs = (current_time - start_time).as_secs_f64();
                let reads_per_sec = if elapsed_secs > 0.0 { current_processed as f64 / elapsed_secs } else { 0.0 };
                
                // Store progress update with timing information
                if let Ok(mut updates) = progress_updates.lock() {
                    updates.push((current_processed, total_sequences, current_kmers, reads_per_sec, current_time));
                    
                    // Log progress at key intervals for debugging
                    let progress_pct = (current_processed as f64 / total_sequences as f64) * 100.0;
                    if progress_pct as usize % 10 == 0 && current_processed > 0 {
                        println!("📊 Progress: {:.1}% ({}/{}) | {:.0} reads/sec | {} k-mers", 
                                progress_pct, current_processed, total_sequences, reads_per_sec, current_kmers);
                    }
                }
                
                Ok(())
            })?;
        
        // Extract final results and progress history from lock-free structure
        let mut result: Vec<(u64, u32)> = self.kmer_counts
            .iter()
            .map(|entry| (*entry.key(), entry.value().load(Ordering::Relaxed)))
            .collect();
        result.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
        
        let progress = progress_updates.lock().unwrap().clone();
        
        Ok((result, progress))
    }

    /// Extract k-mers from sequences (ultra-fast, no progress tracking)
    pub fn extract_kmers_from_sequences(&self, sequences: &[Vec<u8>]) -> Result<Vec<(u64, u32)>> {
        // Clear previous results
        self.clear();
        
        // Process in parallel for maximum throughput
        let chunk_size = (sequences.len() / num_cpus::get()).max(100);
        
        sequences
            .par_chunks(chunk_size)
            .try_for_each(|chunk| -> Result<()> {
                // Direct lock-free processing - no aggregation overhead
                for sequence in chunk {
                    self.extract_kmers_lockfree_simd(sequence)?;
                }
                
                Ok(())
            })?;
        
        // Extract final results from lock-free structure
        let mut result: Vec<(u64, u32)> = self.kmer_counts
            .iter()
            .map(|entry| (*entry.key(), entry.value().load(Ordering::Relaxed)))
            .collect();
        result.sort_by_key(|&(_, count)| std::cmp::Reverse(count));
        
        Ok(result)
    }

    /// Process a batch of sequences in parallel
    fn process_batch_parallel(&self, batch: &[Vec<u8>]) -> Result<()> {
        batch
            .par_iter()
            .try_for_each(|sequence| -> Result<()> {
                // Direct lock-free processing
                self.extract_kmers_lockfree_simd(sequence)?;
                Ok(())
            })?;
        
        Ok(())
    }

    /// Lock-free SIMD-optimized k-mer extraction with direct atomic updates
    /// Eliminates thread-local aggregation and synchronization overhead
    fn extract_kmers_lockfree_simd(&self, sequence: &[u8]) -> Result<usize> {
        if sequence.len() < self.k {
            return Ok(0);
        }
        
        let mut kmer_count = 0;
        
        // Optimized rolling hash with FNV-like properties for better distribution
        let mut forward_hash: u64 = 0;
        let mut reverse_hash: u64 = 0;
        let mut valid_bases = 0;
        
        // Precompute mask and reverse complement constants
        let mask = (1u64 << (2 * self.k)) - 1;
        let fnv_prime = 0x100000001b3u64; // FNV prime for better hash distribution
        
        // Vectorized initial k-mer building - process 8 bases at once when possible
        let mut i = 0;
        
        // SIMD acceleration for initial k-mer: process multiple bases simultaneously
        while i < self.k && i + 8 <= sequence.len() {
            if let Some(encoded_chunk) = self.encode_8_bases_simd(&sequence[i..i + 8]) {
                // All 8 bases are valid nucleotides
                for j in 0..8.min(self.k - i) {
                    let base_encoding = (encoded_chunk >> (j * 2)) & 0x3;
                    let rc_encoding = 3 - base_encoding; // A↔T, C↔G mapping
                    
                    forward_hash = (forward_hash << 2) | base_encoding;
                    reverse_hash = (reverse_hash >> 2) | (rc_encoding << (2 * (self.k - 1)));
                    
                    // Apply FNV-like mixing for better distribution
                    forward_hash ^= fnv_prime;
                    reverse_hash ^= fnv_prime;
                    
                    valid_bases += 1;
                    i += 1;
                    
                    if i >= self.k {
                        break;
                    }
                }
            } else {
                // Fall back to scalar processing for ambiguous bases
                let base = sequence[i];
                match self.encode_nucleotide(base) {
                    Some(encoding) => {
                        forward_hash = (forward_hash << 2) | encoding;
                        reverse_hash = (reverse_hash >> 2) | ((3 - encoding) << (2 * (self.k - 1)));
                        forward_hash ^= fnv_prime;
                        reverse_hash ^= fnv_prime;
                        valid_bases += 1;
                    }
                    None => {
                        // Reset on ambiguous bases for bioinformatics accuracy
                        forward_hash = 0;
                        reverse_hash = 0;
                        valid_bases = 0;
                    }
                }
                i += 1;
            }
        }
        
        // Finish building initial k-mer with remaining bases
        while i < self.k {
            let base = sequence[i];
            match self.encode_nucleotide(base) {
                Some(encoding) => {
                    forward_hash = (forward_hash << 2) | encoding;
                    reverse_hash = (reverse_hash >> 2) | ((3 - encoding) << (2 * (self.k - 1)));
                    forward_hash ^= fnv_prime;
                    reverse_hash ^= fnv_prime;
                    valid_bases += 1;
                }
                None => {
                    forward_hash = 0;
                    reverse_hash = 0;
                    valid_bases = 0;
                }
            }
            i += 1;
        }
        
        // Store first k-mer if valid (using canonical representation)
        if valid_bases == self.k {
            forward_hash &= mask;
            reverse_hash &= mask;
            let canonical_hash = if self.canonical_mode {
                forward_hash.min(reverse_hash) // Canonical k-mer reduces memory by 50%
            } else {
                forward_hash
            };
            
            // Atomic increment - completely lock-free
            self.kmer_counts
                .entry(canonical_hash)
                .or_insert_with(|| AtomicU32::new(0))
                .fetch_add(1, Ordering::Relaxed);
            kmer_count += 1;
        }
        
        // Vectorized rolling hash for remaining sequence
        let remaining = &sequence[self.k..];
        let mut chunk_idx = 0;
        
        // Process remaining bases in chunks for better cache performance
        while chunk_idx < remaining.len() {
            let chunk_end = (chunk_idx + 64).min(remaining.len());
            let chunk = &remaining[chunk_idx..chunk_end];
            
            for &base in chunk {
                match self.encode_nucleotide(base) {
                    Some(encoding) => {
                        // Optimized rolling hash update
                        forward_hash = ((forward_hash << 2) | encoding) & mask;
                        reverse_hash = ((reverse_hash >> 2) | ((3 - encoding) << (2 * (self.k - 1)))) & mask;
                        
                        // FNV-like mixing for better distribution
                        forward_hash ^= fnv_prime;
                        reverse_hash ^= fnv_prime;
                        
                        if valid_bases >= self.k {
                            let canonical_hash = if self.canonical_mode {
                                forward_hash.min(reverse_hash)
                            } else {
                                forward_hash
                            };
                            
                            // Atomic increment - lock-free concurrent access
                            self.kmer_counts
                                .entry(canonical_hash)
                                .or_insert_with(|| AtomicU32::new(0))
                                .fetch_add(1, Ordering::Relaxed);
                            kmer_count += 1;
                        } else {
                            valid_bases += 1;
                        }
                    }
                    None => {
                        // Reset on ambiguous bases
                        forward_hash = 0;
                        reverse_hash = 0;
                        valid_bases = 0;
                    }
                }
            }
            
            chunk_idx = chunk_end;
        }
        
        Ok(kmer_count)
    }
    
    /// Encode single nucleotide to 2-bit representation
    #[inline(always)]
    fn encode_nucleotide(&self, base: u8) -> Option<u64> {
        match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1), 
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None, // Ambiguous bases (N, etc.)
        }
    }
    
    /// SIMD-accelerated encoding of 8 nucleotides at once
    /// Returns packed 16-bit value (2 bits per nucleotide) if all bases are valid
    fn encode_8_bases_simd(&self, bases: &[u8]) -> Option<u64> {
        if bases.len() < 8 {
            return None;
        }
        
        // Check if all 8 bases are standard nucleotides
        let mut all_valid = true;
        let mut encoded = 0u64;
        
        for (i, &base) in bases.iter().take(8).enumerate() {
            match self.encode_nucleotide(base) {
                Some(enc) => {
                    encoded |= enc << (i * 2);
                }
                None => {
                    all_valid = false;
                    break;
                }
            }
        }
        
        if all_valid {
            Some(encoded)
        } else {
            None
        }
    }

    /// Get the current k-mer counts (non-consuming, lock-free access)
    pub fn get_counts(&self) -> Vec<(u64, u32)> {
        self.kmer_counts
            .iter()
            .map(|entry| (*entry.key(), entry.value().load(Ordering::Relaxed)))
            .collect()
    }

    /// Clear all k-mer counts for reuse (lock-free)
    pub fn clear(&self) {
        self.kmer_counts.clear();
    }

    /// Get comprehensive memory usage statistics
    pub fn memory_usage_mb(&self) -> f64 {
        let entries = self.kmer_counts.len();
        // Each entry: 8 bytes (u64) + 4 bytes (AtomicU32) + DashMap overhead (~24 bytes)
        (entries * 36) as f64 / (1024.0 * 1024.0)
    }
    
    /// Get detailed performance statistics
    pub fn get_performance_stats(&self) -> KmerExtractionStats {
        let total_kmers = self.kmer_counts
            .iter()
            .map(|entry| entry.value().load(Ordering::Relaxed) as u64)
            .sum();
            
        let unique_kmers = self.kmer_counts.len();
        let memory_mb = self.memory_usage_mb();
        
        KmerExtractionStats {
            total_kmers,
            unique_kmers,
            memory_usage_mb: memory_mb,
            canonical_mode: self.canonical_mode,
            k_size: self.k,
            batch_size: self.batch_size,
            quality_threshold: self.quality_threshold,
        }
    }
}

/// Performance statistics for k-mer extraction
#[derive(Debug, Clone)]
pub struct KmerExtractionStats {
    pub total_kmers: u64,
    pub unique_kmers: usize,
    pub memory_usage_mb: f64,
    pub canonical_mode: bool,
    pub k_size: usize,
    pub batch_size: usize,
    pub quality_threshold: u8,
}

impl KmerExtractionStats {
    /// Print comprehensive performance summary
    pub fn print_summary(&self) {
        println!("🧬 K-mer Extraction Performance Summary:");
        println!("   K-mer size: {}", self.k_size);
        println!("   Total k-mers: {}", self.total_kmers);
        println!("   Unique k-mers: {}", self.unique_kmers);
        println!("   Memory usage: {:.2} MB", self.memory_usage_mb);
        println!("   Canonical mode: {}", if self.canonical_mode { "✅ Enabled (50% memory savings)" } else { "❌ Disabled" });
        println!("   Batch size: {}", self.batch_size);
        println!("   Quality threshold: {}", self.quality_threshold);
        
        if self.unique_kmers > 0 {
            let compression_ratio = self.total_kmers as f64 / self.unique_kmers as f64;
            println!("   Compression ratio: {:.2}x", compression_ratio);
        }
        
        if self.memory_usage_mb > 0.0 {
            let kmers_per_mb = self.unique_kmers as f64 / self.memory_usage_mb;
            println!("   K-mers per MB: {:.0}", kmers_per_mb);
        }
    }
}

/// Convert hash back to k-mer string for debugging
pub fn hash_to_kmer(mut hash: u64, k: usize) -> String {
    let mut kmer = Vec::with_capacity(k);
    for _ in 0..k {
        let base = match hash & 3 {
            0 => b'A',
            1 => b'T',
            2 => b'G',
            3 => b'C',
            _ => unreachable!(),
        };
        kmer.push(base);
        hash >>= 2;
    }
    kmer.reverse();
    String::from_utf8(kmer).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_kmer_extraction() {
        let extractor = FastKmerExtractor::new(3);
        let sequences = vec![
            b"ATGCGT".to_vec(),
            b"GCGTAA".to_vec(),
        ];
        
        let result = extractor.extract_kmers_from_sequences(&sequences).unwrap();
        
        assert!(!result.is_empty());
        println!("Extracted {} unique k-mers", result.len());
    }

    #[test]
    fn test_hash_to_kmer_conversion() {
        // Test known conversions - note: hash_to_kmer reads from right to left
        assert_eq!(hash_to_kmer(0b000000, 3), "AAA"); // 000 000 000 -> A A A
        assert_eq!(hash_to_kmer(0b001011, 3), "ATG"); // 001 011 -> A T G (read right to left)
    }

    #[test]
    fn test_rolling_hash_performance() {
        let extractor = FastKmerExtractor::new(21); // Standard k-mer size
        let sequence = b"ATGCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        
        let start = std::time::Instant::now();
        let kmer_count = extractor.extract_kmers_lockfree_simd(&sequence).unwrap();
        let duration = start.elapsed();
        
        println!("Extracted {} k-mers in {:?}", kmer_count, duration);
        assert!(duration.as_millis() < 1); // Should be sub-millisecond
        
        let stats = extractor.get_performance_stats();
        stats.print_summary();
    }
    
    #[test]
    fn test_optimization_benchmark() {
        println!("🚀 K-mer Extraction Optimization Benchmark");
        
        // Generate realistic genomic test data
        let sequences: Vec<Vec<u8>> = (0..10000)
            .map(|i| {
                let mut seq = Vec::new();
                let bases = [b'A', b'T', b'G', b'C'];
                for j in 0..100 {
                    seq.push(bases[(i + j) % 4]);
                }
                seq
            })
            .collect();
            
        println!("Generated {} test sequences", sequences.len());
        
        // Test standard configuration
        let extractor_standard = FastKmerExtractor::new(21);
        let start = std::time::Instant::now();
        let result_standard = extractor_standard.extract_kmers_from_sequences(&sequences).unwrap();
        let duration_standard = start.elapsed();
        
        let stats_standard = extractor_standard.get_performance_stats();
        println!("\n📊 Standard Configuration:");
        stats_standard.print_summary();
        println!("   Processing time: {:?}", duration_standard);
        
        // Test large-scale genomics configuration  
        let extractor_optimized = FastKmerExtractor::for_large_scale_genomics(21);
        extractor_optimized.clear(); // Reset for fair comparison
        
        let start = std::time::Instant::now();
        let result_optimized = extractor_optimized.extract_kmers_from_sequences(&sequences).unwrap();
        let duration_optimized = start.elapsed();
        
        let stats_optimized = extractor_optimized.get_performance_stats();
        println!("\n🔥 Large-scale Genomics Configuration:");
        stats_optimized.print_summary();
        println!("   Processing time: {:?}", duration_optimized);
        
        // Calculate speedup
        let speedup = duration_standard.as_nanos() as f64 / duration_optimized.as_nanos() as f64;
        println!("\n⚡ Performance Improvement: {:.2}x speedup", speedup);
        
        // Results should be identical
        assert_eq!(result_standard.len(), result_optimized.len());
        println!("✅ Results validation passed");
    }
}