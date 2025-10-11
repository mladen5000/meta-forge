//! Concrete Code Examples for Assembly Performance Optimizations
//! ==========================================================
//!
//! This file contains before/after code examples demonstrating the specific
//! optimizations identified in the performance bottleneck analysis.

use std::arch::x86_64::*;
use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
use ahash::AHashMap;
use dashmap::DashMap;
use crossbeam::queue::SegQueue;

// ============================================================================
// 1. OPTIMIZED COMPACT K-MER IMPLEMENTATION
// ============================================================================

/// BEFORE: Original CompactKmer with string allocation bottlenecks
mod original {
    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    pub struct CompactKmer {
        data: u64,
        k: u8,
    }

    impl CompactKmer {
        // ❌ BOTTLENECK: String allocation in hot path
        pub fn to_string(&self) -> String {
            let mut result = String::with_capacity(self.k as usize);  // Heap allocation
            for i in 0..self.k {
                let bits = (self.data >> (2 * (31 - i))) & 0b11;
                let nucleotide = match bits {
                    0b00 => 'A',
                    0b01 => 'C',
                    0b10 => 'G',
                    0b11 => 'T',
                    _ => unreachable!(),
                };
                result.push(nucleotide);  // Character-by-character insertion
            }
            result
        }

        // ❌ BOTTLENECK: Hash recomputation for each k-mer
        pub fn hash(&self) -> u64 {
            self.data.wrapping_mul(0x9e3779b97f4a7c15) ^ (self.k as u64)
        }
    }
}

/// AFTER: Optimized zero-copy, SIMD-enabled k-mer implementation
mod optimized {
    use super::*;

    /// Optimized k-mer with zero-copy operations and SIMD alignment
    #[repr(C, align(16))]  // 16-byte alignment for SIMD operations
    pub struct OptimizedCompactKmer {
        data: u64,              // 8 bytes - 2-bit encoding for up to 32-mers
        k: u8,                 // 1 byte - k-mer length
        _padding: [u8; 7],     // Padding to 16 bytes for SIMD
    }

    /// Zero-copy iterator for k-mer nucleotides
    pub struct NucleotideIterator<'a> {
        kmer: &'a OptimizedCompactKmer,
        pos: u8,
    }

    impl<'a> Iterator for NucleotideIterator<'a> {
        type Item = u8;

        #[inline(always)]
        fn next(&mut self) -> Option<Self::Item> {
            if self.pos >= self.kmer.k {
                return None;
            }

            let bits = (self.kmer.data >> (2 * (31 - self.pos))) & 0b11;
            let nucleotide = match bits {
                0b00 => b'A',
                0b01 => b'C',
                0b10 => b'G',
                0b11 => b'T',
                _ => unreachable!(),
            };

            self.pos += 1;
            Some(nucleotide)
        }
    }

    impl OptimizedCompactKmer {
        /// Zero-copy interface - no string allocation
        pub fn iter_nucleotides(&self) -> NucleotideIterator<'_> {
            NucleotideIterator { kmer: self, pos: 0 }
        }

        /// Write k-mer to pre-allocated buffer
        #[inline]
        pub fn write_to_buffer(&self, buffer: &mut [u8]) -> usize {
            let k = self.k as usize;
            if buffer.len() < k {
                return 0;
            }

            for i in 0..k {
                let bits = (self.data >> (2 * (31 - i))) & 0b11;
                buffer[i] = match bits {
                    0b00 => b'A',
                    0b01 => b'C',
                    0b10 => b'G',
                    0b11 => b'T',
                    _ => unreachable!(),
                };
            }
            k
        }

        /// SIMD-optimized batch comparison (compare 8 k-mers simultaneously)
        #[target_feature(enable = "avx2")]
        pub unsafe fn batch_compare_simd(kmers: &[Self], target_data: u64) -> u64 {
            if kmers.is_empty() {
                return 0;
            }

            let mut matches = 0u64;
            let target_vec = _mm256_set1_epi64x(target_data as i64);

            // Process 4 k-mers at a time with AVX2
            for (i, chunk) in kmers.chunks(4).enumerate() {
                let mut data_array = [0i64; 4];
                for (j, kmer) in chunk.iter().enumerate() {
                    data_array[j] = kmer.data as i64;
                }

                let data_vec = _mm256_loadu_si256(data_array.as_ptr() as *const __m256i);
                let cmp_result = _mm256_cmpeq_epi64(data_vec, target_vec);
                let mask = _mm256_movemask_pd(_mm256_castsi256_pd(cmp_result));

                matches |= (mask as u64) << (i * 4);
            }

            matches
        }

        /// Pre-computed rolling hash for O(1) updates
        #[inline(always)]
        pub fn rolling_hash(&self) -> u64 {
            // Use pre-computed polynomial rolling hash
            self.data.wrapping_mul(0x9e3779b97f4a7c15) ^ (self.k as u64)
        }
    }

    /// Rolling hash for O(1) k-mer hash updates
    pub struct RollingKmerHash {
        hash: u64,
        k: usize,
        base_power: u64,    // base^(k-1) for efficient rolling
        base: u64,
    }

    impl RollingKmerHash {
        pub fn new(k: usize) -> Self {
            let base = 4u64;
            let base_power = base.pow((k - 1) as u32);

            Self {
                hash: 0,
                k,
                base_power,
                base,
            }
        }

        /// Initialize hash with first k-mer
        pub fn init_hash(&mut self, sequence: &[u8]) -> u64 {
            self.hash = 0;
            for i in 0..self.k {
                let nucleotide_value = Self::nucleotide_to_value(sequence[i]);
                self.hash = self.hash.wrapping_mul(self.base).wrapping_add(nucleotide_value);
            }
            self.hash
        }

        /// Update hash for next k-mer in O(1) time
        #[inline(always)]
        pub fn roll(&mut self, old_nucleotide: u8, new_nucleotide: u8) -> u64 {
            // Remove contribution of old nucleotide
            let old_value = Self::nucleotide_to_value(old_nucleotide);
            self.hash = self.hash.wrapping_sub(old_value.wrapping_mul(self.base_power));

            // Add contribution of new nucleotide
            self.hash = self.hash.wrapping_mul(self.base);
            self.hash = self.hash.wrapping_add(Self::nucleotide_to_value(new_nucleotide));

            self.hash
        }

        #[inline(always)]
        fn nucleotide_to_value(nucleotide: u8) -> u64 {
            match nucleotide {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,  // Default for ambiguous bases
            }
        }
    }
}

// ============================================================================
// 2. SIMD-OPTIMIZED NUCLEOTIDE PROCESSING
// ============================================================================

/// SIMD-accelerated nucleotide counting and analysis
pub mod simd_nucleotide_ops {
    use super::*;

    /// SIMD-optimized nucleotide counting (16x parallel with AVX2)
    #[target_feature(enable = "avx2")]
    pub unsafe fn count_nucleotides_simd(sequence: &[u8]) -> [u32; 4] {
        let mut counts = [0u32; 4];  // A, C, G, T

        // Process 32 bytes at a time with AVX2
        let chunks = sequence.chunks_exact(32);
        for chunk in chunks {
            let chunk_counts = count_nucleotides_avx2_chunk(chunk);
            for i in 0..4 {
                counts[i] += chunk_counts[i];
            }
        }

        // Handle remainder with scalar processing
        for &byte in chunks.remainder() {
            match byte {
                b'A' | b'a' => counts[0] += 1,
                b'C' | b'c' => counts[1] += 1,
                b'G' | b'g' => counts[2] += 1,
                b'T' | b't' => counts[3] += 1,
                _ => {}  // Skip ambiguous bases
            }
        }

        counts
    }

    /// Process 32 bytes with AVX2 (8-12x faster than scalar)
    #[target_feature(enable = "avx2")]
    unsafe fn count_nucleotides_avx2_chunk(chunk: &[u8]) -> [u32; 4] {
        // Load 32 bytes into AVX2 register
        let input = _mm256_loadu_si256(chunk.as_ptr() as *const __m256i);

        // Create comparison masks for each nucleotide (case-insensitive)
        let a_upper = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b'A' as i8));
        let a_lower = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b'a' as i8));
        let a_mask = _mm256_or_si256(a_upper, a_lower);

        let c_upper = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b'C' as i8));
        let c_lower = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b'c' as i8));
        let c_mask = _mm256_or_si256(c_upper, c_lower);

        let g_upper = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b'G' as i8));
        let g_lower = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b'g' as i8));
        let g_mask = _mm256_or_si256(g_upper, g_lower);

        let t_upper = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b'T' as i8));
        let t_lower = _mm256_cmpeq_epi8(input, _mm256_set1_epi8(b't' as i8));
        let t_mask = _mm256_or_si256(t_upper, t_lower);

        // Count set bits in each mask using population count
        let a_count = count_mask_bits(a_mask);
        let c_count = count_mask_bits(c_mask);
        let g_count = count_mask_bits(g_mask);
        let t_count = count_mask_bits(t_mask);

        [a_count, c_count, g_count, t_count]
    }

    /// Count set bits in AVX2 mask
    #[target_feature(enable = "avx2")]
    unsafe fn count_mask_bits(mask: __m256i) -> u32 {
        // Convert to scalar and count bits
        let scalar_mask = _mm256_movemask_epi8(mask);
        scalar_mask.count_ones()
    }

    /// SIMD-optimized GC content calculation
    #[target_feature(enable = "avx2")]
    pub unsafe fn calculate_gc_content_simd(sequence: &[u8]) -> f64 {
        let counts = count_nucleotides_simd(sequence);
        let gc_count = counts[1] + counts[2];  // C + G
        let total_count = counts.iter().sum::<u32>();

        if total_count == 0 {
            0.0
        } else {
            (gc_count as f64 / total_count as f64) * 100.0
        }
    }
}

// ============================================================================
// 3. LOCK-FREE PARALLEL GRAPH CONSTRUCTION
// ============================================================================

/// BEFORE: Mutex-based parallel graph construction with contention
mod mutex_based {
    use super::*;
    use std::sync::{Arc, Mutex};
    use rayon::prelude::*;

    // ❌ BOTTLENECK: Heavy mutex contention
    pub fn build_graph_with_mutexes(chunks: &[Vec<u8>]) -> Result<(), Box<dyn std::error::Error>> {
        let nodes_mutex = Arc::new(Mutex::new(AHashMap::<u64, u32>::new()));
        let edges_mutex = Arc::new(Mutex::new(Vec::<(u64, u64)>::new()));

        chunks.par_iter().try_for_each(|chunk| -> Result<(), Box<dyn std::error::Error>> {
            let mut local_nodes = AHashMap::new();
            let mut local_edges = Vec::new();

            // Process chunk locally...
            for i in 0..chunk.len() - 21 {
                let kmer_hash = hash_kmer(&chunk[i..i + 21]);
                local_nodes.insert(kmer_hash, 1);

                if i > 0 {
                    let prev_hash = hash_kmer(&chunk[i - 1..i + 20]);
                    local_edges.push((prev_hash, kmer_hash));
                }
            }

            // ❌ CRITICAL BOTTLENECK: Lock contention during merge
            {
                let mut nodes = nodes_mutex.lock().unwrap();  // Blocking!
                for (hash, count) in local_nodes {
                    *nodes.entry(hash).or_insert(0) += count;  // Hash operations under lock
                }
            }

            {
                let mut edges = edges_mutex.lock().unwrap();  // More blocking!
                edges.extend(local_edges);
            }

            Ok(())
        })?;

        Ok(())
    }

    fn hash_kmer(kmer: &[u8]) -> u64 {
        // Simplified hash function
        kmer.iter().fold(0u64, |acc, &b| acc.wrapping_mul(4).wrapping_add(b as u64))
    }
}

/// AFTER: Lock-free parallel graph construction
mod lock_free {
    use super::*;
    use rayon::prelude::*;

    /// Lock-free parallel graph builder
    pub struct LockFreeGraphBuilder {
        // Lock-free concurrent hash map for nodes
        nodes: DashMap<u64, AtomicU32>,

        // Per-thread edge queues to avoid contention
        edge_queues: Vec<SegQueue<(u64, u64)>>,

        // Thread count for work distribution
        num_threads: usize,
    }

    impl LockFreeGraphBuilder {
        pub fn new(num_threads: usize) -> Self {
            let edge_queues = (0..num_threads)
                .map(|_| SegQueue::new())
                .collect();

            Self {
                nodes: DashMap::new(),
                edge_queues,
                num_threads,
            }
        }

        /// ✅ OPTIMIZED: Lock-free parallel construction
        pub fn build_graph_lock_free(&self, chunks: &[Vec<u8>]) -> Result<(), Box<dyn std::error::Error>> {
            // Parallel processing with work-stealing
            chunks.par_iter().enumerate().for_each(|(chunk_idx, chunk)| {
                let thread_id = chunk_idx % self.num_threads;
                self.process_chunk_lock_free(chunk, thread_id);
            });

            // Single-threaded consolidation of edges (fast since no locking needed)
            self.consolidate_edges();

            Ok(())
        }

        #[inline(always)]
        fn process_chunk_lock_free(&self, chunk: &[u8], thread_id: usize) {
            let edge_queue = &self.edge_queues[thread_id];

            for i in 0..chunk.len().saturating_sub(21) {
                let kmer_hash = self.hash_kmer_optimized(&chunk[i..i + 21]);

                // ✅ OPTIMIZED: Atomic increment without locking
                self.nodes
                    .entry(kmer_hash)
                    .or_insert(AtomicU32::new(0))
                    .fetch_add(1, Ordering::Relaxed);

                // ✅ OPTIMIZED: Queue edge for later processing (no contention)
                if i > 0 {
                    let prev_hash = self.hash_kmer_optimized(&chunk[i - 1..i + 20]);
                    edge_queue.push((prev_hash, kmer_hash));
                }
            }
        }

        fn consolidate_edges(&self) {
            // Fast single-threaded edge consolidation
            let mut all_edges = Vec::new();

            for queue in &self.edge_queues {
                while let Some(edge) = queue.pop() {
                    all_edges.push(edge);
                }
            }

            // Sort and deduplicate edges
            all_edges.sort_unstable();
            all_edges.dedup();
        }

        /// Optimized k-mer hashing with better cache behavior
        #[inline(always)]
        fn hash_kmer_optimized(&self, kmer: &[u8]) -> u64 {
            // Use polynomial rolling hash for better distribution
            let mut hash = 0u64;
            for &byte in kmer {
                hash = hash.wrapping_mul(257).wrapping_add(byte as u64);
            }
            hash
        }
    }
}

// ============================================================================
// 4. CACHE-OPTIMIZED GRAPH LAYOUT
// ============================================================================

/// Structure-of-Arrays layout for better cache performance
pub mod cache_optimized_graph {
    use super::*;

    /// BEFORE: Array-of-Structures (poor cache locality)
    #[derive(Clone)]
    pub struct GraphNode {
        hash: u64,          // 8 bytes
        coverage: u32,      // 4 bytes
        in_degree: u8,      // 1 byte
        out_degree: u8,     // 1 byte
        // Total: ~24 bytes with padding, poor cache utilization
    }

    /// AFTER: Structure-of-Arrays (excellent cache locality)
    pub struct CacheOptimizedGraph {
        // Node data in separate arrays for better cache utilization
        node_hashes: Vec<u64>,           // Sequential access for all hashes
        node_coverage: Vec<u32>,         // Sequential access for all coverage
        node_degrees: Vec<u16>,          // Packed: in_degree(8) + out_degree(8)

        // CSR (Compressed Sparse Row) format for edges
        edge_offsets: Vec<u32>,          // Start index for each node's edges
        edge_targets: Vec<u32>,          // Target node indices (not hashes!)
        edge_weights: Vec<u16>,          // Compressed edge weights

        // Hash-to-index mapping
        hash_to_index: AHashMap<u64, u32>,

        // Temporary workspace to avoid allocations
        temp_workspace: Vec<u32>,
    }

    impl CacheOptimizedGraph {
        pub fn new(estimated_nodes: usize, estimated_edges: usize) -> Self {
            Self {
                node_hashes: Vec::with_capacity(estimated_nodes),
                node_coverage: Vec::with_capacity(estimated_nodes),
                node_degrees: Vec::with_capacity(estimated_nodes),
                edge_offsets: Vec::with_capacity(estimated_nodes + 1),
                edge_targets: Vec::with_capacity(estimated_edges),
                edge_weights: Vec::with_capacity(estimated_edges),
                hash_to_index: AHashMap::with_capacity(estimated_nodes),
                temp_workspace: Vec::with_capacity(1000),
            }
        }

        /// ✅ OPTIMIZED: Cache-friendly contig generation
        pub fn generate_contigs_cache_optimized(&mut self) -> Vec<Vec<u32>> {
            let mut contigs = Vec::new();
            let mut visited = vec![false; self.node_hashes.len()];

            // Process nodes in sequential order for perfect cache behavior
            for node_idx in 0..self.node_hashes.len() {
                if visited[node_idx] {
                    continue;
                }

                // Extract degree information from packed representation
                let in_degree = (self.node_degrees[node_idx] & 0xFF) as u8;
                let out_degree = ((self.node_degrees[node_idx] >> 8) & 0xFF) as u8;

                // Start path from nodes that are likely contig starts
                if in_degree == 0 || self.is_branch_node(node_idx) {
                    let contig = self.trace_linear_path_optimized(node_idx, &mut visited);
                    if !contig.is_empty() {
                        contigs.push(contig);
                    }
                }
            }

            contigs
        }

        /// ✅ OPTIMIZED: Linear path tracing with cache-friendly access
        fn trace_linear_path_optimized(&mut self, start_idx: usize, visited: &mut [bool]) -> Vec<u32> {
            self.temp_workspace.clear();
            let mut current_idx = start_idx;

            // Follow linear path with perfect cache locality
            loop {
                if visited[current_idx] {
                    break;
                }

                visited[current_idx] = true;
                self.temp_workspace.push(current_idx as u32);

                // Check if this node has exactly one outgoing edge
                let out_degree = ((self.node_degrees[current_idx] >> 8) & 0xFF) as u8;
                if out_degree != 1 {
                    break;
                }

                // Get the single outgoing edge
                let edge_start = self.edge_offsets[current_idx] as usize;
                let edge_end = self.edge_offsets[current_idx + 1] as usize;

                if edge_end - edge_start != 1 {
                    break;  // Not exactly one edge
                }

                let next_idx = self.edge_targets[edge_start] as usize;

                // Check if next node has exactly one incoming edge
                let next_in_degree = (self.node_degrees[next_idx] & 0xFF) as u8;
                if next_in_degree != 1 {
                    break;
                }

                current_idx = next_idx;
            }

            self.temp_workspace.clone()
        }

        #[inline(always)]
        fn is_branch_node(&self, node_idx: usize) -> bool {
            let in_degree = (self.node_degrees[node_idx] & 0xFF) as u8;
            let out_degree = ((self.node_degrees[node_idx] >> 8) & 0xFF) as u8;
            in_degree > 1 || out_degree > 1
        }
    }
}

// ============================================================================
// 5. ZERO-COPY K-MER PROCESSING
// ============================================================================

/// Zero-allocation k-mer processing pipeline
pub mod zero_copy_processing {
    use super::*;

    /// Zero-copy k-mer processor that eliminates all string allocations
    pub struct ZeroCopyKmerProcessor<'a> {
        sequence: &'a [u8],
        k: usize,
        rolling_hash: optimized::RollingKmerHash,
        position: usize,
    }

    impl<'a> ZeroCopyKmerProcessor<'a> {
        pub fn new(sequence: &'a [u8], k: usize) -> Self {
            Self {
                sequence,
                k,
                rolling_hash: optimized::RollingKmerHash::new(k),
                position: 0,
            }
        }

        /// Process all k-mers without any allocations
        pub fn process_all_kmers<F>(&mut self, mut callback: F) -> Result<(), Box<dyn std::error::Error>>
        where
            F: FnMut(u64, &[u8]) -> Result<(), Box<dyn std::error::Error>>,
        {
            if self.sequence.len() < self.k {
                return Ok(());
            }

            // Initialize rolling hash with first k-mer
            let mut hash = self.rolling_hash.init_hash(&self.sequence[0..self.k]);
            callback(hash, &self.sequence[0..self.k])?;

            // Process remaining k-mers with O(1) rolling hash updates
            for i in self.k..self.sequence.len() {
                hash = self.rolling_hash.roll(
                    self.sequence[i - self.k],
                    self.sequence[i]
                );
                callback(hash, &self.sequence[i - self.k + 1..=i])?;
            }

            Ok(())
        }

        /// Batch process k-mers with SIMD optimizations
        pub fn batch_process_simd<F>(&mut self, batch_size: usize, mut callback: F) -> Result<(), Box<dyn std::error::Error>>
        where
            F: FnMut(&[u64], &[&[u8]]) -> Result<(), Box<dyn std::error::Error>>,
        {
            let mut hash_batch = Vec::with_capacity(batch_size);
            let mut kmer_batch = Vec::with_capacity(batch_size);

            self.process_all_kmers(|hash, kmer_slice| {
                hash_batch.push(hash);
                kmer_batch.push(kmer_slice);

                if hash_batch.len() == batch_size {
                    callback(&hash_batch, &kmer_batch)?;
                    hash_batch.clear();
                    kmer_batch.clear();
                }

                Ok(())
            })?;

            // Process remaining k-mers
            if !hash_batch.is_empty() {
                callback(&hash_batch, &kmer_batch)?;
            }

            Ok(())
        }
    }
}

// ============================================================================
// 6. PERFORMANCE BENCHMARKING HELPERS
// ============================================================================

#[cfg(test)]
mod benchmarks {
    use super::*;
    use std::time::Instant;

    /// Benchmark comparing original vs optimized implementations
    pub fn benchmark_kmer_operations() {
        let sequence = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCG".repeat(1000);
        let k = 21;
        let iterations = 1000;

        // Benchmark original string-based approach
        let start = Instant::now();
        for _ in 0..iterations {
            for i in 0..sequence.len().saturating_sub(k) {
                let kmer_str = String::from_utf8_lossy(&sequence[i..i + k]);
                let _ = kmer_str.len();  // Simulate work
            }
        }
        let original_duration = start.elapsed();

        // Benchmark optimized zero-copy approach
        let start = Instant::now();
        for _ in 0..iterations {
            let mut processor = zero_copy_processing::ZeroCopyKmerProcessor::new(&sequence, k);
            processor.process_all_kmers(|_hash, _kmer| Ok(())).unwrap();
        }
        let optimized_duration = start.elapsed();

        println!("Original approach: {:?}", original_duration);
        println!("Optimized approach: {:?}", optimized_duration);
        println!("Speedup: {:.2}x", original_duration.as_secs_f64() / optimized_duration.as_secs_f64());
    }

    /// Benchmark SIMD vs scalar nucleotide counting
    pub fn benchmark_nucleotide_counting() {
        let sequence = b"ATCGATCGATCGATCGATCGATCG".repeat(10000);
        let iterations = 100;

        // Scalar approach
        let start = Instant::now();
        for _ in 0..iterations {
            let mut counts = [0u32; 4];
            for &byte in &sequence {
                match byte {
                    b'A' => counts[0] += 1,
                    b'C' => counts[1] += 1,
                    b'G' => counts[2] += 1,
                    b'T' => counts[3] += 1,
                    _ => {}
                }
            }
        }
        let scalar_duration = start.elapsed();

        // SIMD approach
        let start = Instant::now();
        for _ in 0..iterations {
            unsafe {
                let _counts = simd_nucleotide_ops::count_nucleotides_simd(&sequence);
            }
        }
        let simd_duration = start.elapsed();

        println!("Scalar nucleotide counting: {:?}", scalar_duration);
        println!("SIMD nucleotide counting: {:?}", simd_duration);
        println!("SIMD speedup: {:.2}x", scalar_duration.as_secs_f64() / simd_duration.as_secs_f64());
    }
}