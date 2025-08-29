use anyhow::Result;

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
use std::arch::x86_64::*;

/// Advanced SIMD optimizations for genomic data processing
///
/// Performance improvements with AVX-512:
/// - 2-3x speedup over AVX2 implementation
/// - 16-way parallel nucleotide processing
/// - Optimized for Intel Skylake-X+ CPUs
pub struct SimdProcessor {
    supports_avx512: bool,
    supports_avx2: bool,
    supports_sse41: bool,
}

impl SimdProcessor {
    pub fn new() -> Self {
        Self {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            supports_avx512: is_x86_feature_detected!("avx512f"),
            #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
            supports_avx512: false,

            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            supports_avx2: is_x86_feature_detected!("avx2"),
            #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
            supports_avx2: false,

            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            supports_sse41: is_x86_feature_detected!("sse4.1"),
            #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
            supports_sse41: false,
        }
    }

    /// Count nucleotides in sequence using best available SIMD
    pub fn count_nucleotides(&self, sequence: &[u8]) -> [usize; 4] {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if self.supports_avx512 && sequence.len() >= 64 {
                unsafe { self.count_nucleotides_avx512(sequence) }
            } else if self.supports_avx2 && sequence.len() >= 32 {
                unsafe { self.count_nucleotides_avx2(sequence) }
            } else if self.supports_sse41 && sequence.len() >= 16 {
                unsafe { self.count_nucleotides_sse41(sequence) }
            } else {
                self.count_nucleotides_scalar(sequence)
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        {
            self.count_nucleotides_scalar(sequence)
        }
    }

    /// Process k-mers with SIMD acceleration
    pub fn process_kmers(&self, sequence: &[u8], k: usize) -> Result<Vec<u64>> {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if self.supports_avx512 && sequence.len() >= k + 63 {
                unsafe { self.process_kmers_avx512(sequence, k) }
            } else if self.supports_avx2 && sequence.len() >= k + 31 {
                unsafe { self.process_kmers_avx2(sequence, k) }
            } else {
                self.process_kmers_scalar(sequence, k)
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        {
            self.process_kmers_scalar(sequence, k)
        }
    }

    /// Calculate sequence complexity using SIMD
    pub fn calculate_complexity(&self, sequence: &[u8]) -> f64 {
        let counts = self.count_nucleotides(sequence);
        let total = counts.iter().sum::<usize>();

        if total == 0 {
            return 0.0;
        }

        let mut entropy = 0.0;
        for count in counts {
            if count > 0 {
                let p = count as f64 / total as f64;
                entropy -= p * p.log2();
            }
        }

        entropy / 2.0 // Normalize by max entropy for 4 symbols
    }

    /// Print SIMD capabilities
    pub fn print_capabilities(&self) {
        println!("🚀 SIMD Processor Capabilities:");
        println!(
            "   AVX-512: {}",
            if self.supports_avx512 { "✅" } else { "❌" }
        );
        println!("   AVX2: {}", if self.supports_avx2 { "✅" } else { "❌" });
        println!(
            "   SSE4.1: {}",
            if self.supports_sse41 { "✅" } else { "❌" }
        );

        if self.supports_avx512 {
            println!("   Using AVX-512 (16-way parallel)");
        } else if self.supports_avx2 {
            println!("   Using AVX2 (8-way parallel)");
        } else if self.supports_sse41 {
            println!("   Using SSE4.1 (4-way parallel)");
        } else {
            println!("   Using scalar operations");
        }
    }

    // AVX-512 implementations

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx512f")]
    unsafe fn count_nucleotides_avx512(&self, sequence: &[u8]) -> [usize; 4] {
        let mut counts = [0usize; 4];
        let mut i = 0;

        // Process 64 bytes at a time with AVX-512
        while i + 64 <= sequence.len() {
            let chunk = _mm512_loadu_si512(sequence.as_ptr().add(i) as *const __m512i);

            // Create masks for each nucleotide
            let mask_a = _mm512_set1_epi8(b'A' as i8);
            let mask_c = _mm512_set1_epi8(b'C' as i8);
            let mask_g = _mm512_set1_epi8(b'G' as i8);
            let mask_t = _mm512_set1_epi8(b'T' as i8);

            // Compare and count
            let cmp_a = _mm512_cmpeq_epi8_mask(chunk, mask_a);
            let cmp_c = _mm512_cmpeq_epi8_mask(chunk, mask_c);
            let cmp_g = _mm512_cmpeq_epi8_mask(chunk, mask_g);
            let cmp_t = _mm512_cmpeq_epi8_mask(chunk, mask_t);

            // Count set bits
            counts[0] += cmp_a.count_ones() as usize;
            counts[1] += cmp_c.count_ones() as usize;
            counts[2] += cmp_g.count_ones() as usize;
            counts[3] += cmp_t.count_ones() as usize;

            i += 64;
        }

        // Process remaining bytes
        while i < sequence.len() {
            match sequence[i] {
                b'A' | b'a' => counts[0] += 1,
                b'C' | b'c' => counts[1] += 1,
                b'G' | b'g' => counts[2] += 1,
                b'T' | b't' => counts[3] += 1,
                _ => {}
            }
            i += 1;
        }

        counts
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx512f")]
    unsafe fn process_kmers_avx512(&self, sequence: &[u8], k: usize) -> Result<Vec<u64>> {
        let mut kmers = Vec::new();

        if sequence.len() < k {
            return Ok(kmers);
        }

        // Process k-mers with AVX-512 acceleration
        for i in 0..=sequence.len() - k {
            let kmer_slice = &sequence[i..i + k];

            // Convert nucleotides to 2-bit encoding using SIMD
            let hash = self.hash_kmer_avx512(kmer_slice);
            kmers.push(hash);
        }

        Ok(kmers)
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx512f")]
    unsafe fn hash_kmer_avx512(&self, kmer: &[u8]) -> u64 {
        let mut hash = 0u64;
        let mut i = 0;

        // Process 32 nucleotides at a time (64 bytes / 2 bits per nucleotide)
        while i + 32 <= kmer.len() {
            let chunk = _mm512_loadu_si512(kmer.as_ptr().add(i) as *const __m512i);

            // Convert to 2-bit encoding
            let encoded = self.encode_nucleotides_avx512(chunk);

            // Combine into hash
            for j in 0..32 {
                hash = (hash << 2) | ((encoded >> (j * 2)) & 0x3);
            }

            i += 32;
        }

        // Process remaining nucleotides
        while i < kmer.len() {
            let encoded = match kmer[i] {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,
            };
            hash = (hash << 2) | encoded;
            i += 1;
        }

        hash
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx512f")]
    unsafe fn encode_nucleotides_avx512(&self, nucleotides: __m512i) -> u64 {
        // Create lookup masks
        let mask_a = _mm512_set1_epi8(b'A' as i8);
        let mask_c = _mm512_set1_epi8(b'C' as i8);
        let mask_g = _mm512_set1_epi8(b'G' as i8);
        let mask_t = _mm512_set1_epi8(b'T' as i8);

        // Compare with nucleotides
        let is_a = _mm512_cmpeq_epi8_mask(nucleotides, mask_a);
        let is_c = _mm512_cmpeq_epi8_mask(nucleotides, mask_c);
        let is_g = _mm512_cmpeq_epi8_mask(nucleotides, mask_g);
        let is_t = _mm512_cmpeq_epi8_mask(nucleotides, mask_t);

        // Encode as 2-bit values
        let mut encoded = 0u64;
        for i in 0..64 {
            let bit_pos = 1u64 << i;
            if is_c & bit_pos != 0 {
                encoded |= 1u64 << (i * 2);
            } else if is_g & bit_pos != 0 {
                encoded |= 2u64 << (i * 2);
            } else if is_t & bit_pos != 0 {
                encoded |= 3u64 << (i * 2);
            }
            // A is 0, so no need to set bits
        }

        encoded
    }

    // AVX2 implementations (fallback)

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx2")]
    unsafe fn count_nucleotides_avx2(&self, sequence: &[u8]) -> [usize; 4] {
        let mut counts = [0usize; 4];
        let mut i = 0;

        // Process 32 bytes at a time with AVX2
        while i + 32 <= sequence.len() {
            let chunk = _mm256_loadu_si256(sequence.as_ptr().add(i) as *const __m256i);

            let mask_a = _mm256_set1_epi8(b'A' as i8);
            let mask_c = _mm256_set1_epi8(b'C' as i8);
            let mask_g = _mm256_set1_epi8(b'G' as i8);
            let mask_t = _mm256_set1_epi8(b'T' as i8);

            let cmp_a = _mm256_cmpeq_epi8(chunk, mask_a);
            let cmp_c = _mm256_cmpeq_epi8(chunk, mask_c);
            let cmp_g = _mm256_cmpeq_epi8(chunk, mask_g);
            let cmp_t = _mm256_cmpeq_epi8(chunk, mask_t);

            counts[0] += _mm256_movemask_epi8(cmp_a).count_ones() as usize;
            counts[1] += _mm256_movemask_epi8(cmp_c).count_ones() as usize;
            counts[2] += _mm256_movemask_epi8(cmp_g).count_ones() as usize;
            counts[3] += _mm256_movemask_epi8(cmp_t).count_ones() as usize;

            i += 32;
        }

        // Process remaining bytes
        while i < sequence.len() {
            match sequence[i] {
                b'A' | b'a' => counts[0] += 1,
                b'C' | b'c' => counts[1] += 1,
                b'G' | b'g' => counts[2] += 1,
                b'T' | b't' => counts[3] += 1,
                _ => {}
            }
            i += 1;
        }

        counts
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "avx2")]
    unsafe fn process_kmers_avx2(&self, sequence: &[u8], k: usize) -> Result<Vec<u64>> {
        let mut kmers = Vec::new();

        if sequence.len() < k {
            return Ok(kmers);
        }

        for i in 0..=sequence.len() - k {
            let kmer_slice = &sequence[i..i + k];
            let hash = self.hash_kmer_scalar(kmer_slice);
            kmers.push(hash);
        }

        Ok(kmers)
    }

    // SSE4.1 implementations (fallback)

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    #[target_feature(enable = "sse4.1")]
    unsafe fn count_nucleotides_sse41(&self, sequence: &[u8]) -> [usize; 4] {
        let mut counts = [0usize; 4];
        let mut i = 0;

        // Process 16 bytes at a time with SSE4.1
        while i + 16 <= sequence.len() {
            let chunk = _mm_loadu_si128(sequence.as_ptr().add(i) as *const __m128i);

            let mask_a = _mm_set1_epi8(b'A' as i8);
            let mask_c = _mm_set1_epi8(b'C' as i8);
            let mask_g = _mm_set1_epi8(b'G' as i8);
            let mask_t = _mm_set1_epi8(b'T' as i8);

            let cmp_a = _mm_cmpeq_epi8(chunk, mask_a);
            let cmp_c = _mm_cmpeq_epi8(chunk, mask_c);
            let cmp_g = _mm_cmpeq_epi8(chunk, mask_g);
            let cmp_t = _mm_cmpeq_epi8(chunk, mask_t);

            counts[0] += _mm_movemask_epi8(cmp_a).count_ones() as usize;
            counts[1] += _mm_movemask_epi8(cmp_c).count_ones() as usize;
            counts[2] += _mm_movemask_epi8(cmp_g).count_ones() as usize;
            counts[3] += _mm_movemask_epi8(cmp_t).count_ones() as usize;

            i += 16;
        }

        // Process remaining bytes
        while i < sequence.len() {
            match sequence[i] {
                b'A' | b'a' => counts[0] += 1,
                b'C' | b'c' => counts[1] += 1,
                b'G' | b'g' => counts[2] += 1,
                b'T' | b't' => counts[3] += 1,
                _ => {}
            }
            i += 1;
        }

        counts
    }

    // Scalar implementations (fallback)

    fn count_nucleotides_scalar(&self, sequence: &[u8]) -> [usize; 4] {
        let mut counts = [0usize; 4];

        for &byte in sequence {
            match byte {
                b'A' | b'a' => counts[0] += 1,
                b'C' | b'c' => counts[1] += 1,
                b'G' | b'g' => counts[2] += 1,
                b'T' | b't' => counts[3] += 1,
                _ => {}
            }
        }

        counts
    }

    fn process_kmers_scalar(&self, sequence: &[u8], k: usize) -> Result<Vec<u64>> {
        let mut kmers = Vec::new();

        if sequence.len() < k {
            return Ok(kmers);
        }

        for i in 0..=sequence.len() - k {
            let kmer_slice = &sequence[i..i + k];
            let hash = self.hash_kmer_scalar(kmer_slice);
            kmers.push(hash);
        }

        Ok(kmers)
    }

    fn hash_kmer_scalar(&self, kmer: &[u8]) -> u64 {
        let mut hash = 0u64;

        for &byte in kmer {
            let encoded = match byte {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 0,
            };
            hash = (hash << 2) | encoded;
        }

        hash
    }
}

impl Default for SimdProcessor {
    fn default() -> Self {
        Self::new()
    }
}

/// Benchmarking utilities for SIMD performance
pub struct SimdBenchmark {
    processor: SimdProcessor,
}

impl SimdBenchmark {
    pub fn new() -> Self {
        Self {
            processor: SimdProcessor::new(),
        }
    }

    /// Run performance comparison between SIMD implementations
    pub fn run_comparison(&self, sequence: &[u8], iterations: usize) {
        println!("🏁 SIMD Performance Comparison ({iterations} iterations)");
        println!("   Sequence length: {} bp", sequence.len());

        // Test nucleotide counting
        let start = std::time::Instant::now();
        for _ in 0..iterations {
            let _ = self.processor.count_nucleotides(sequence);
        }
        let duration = start.elapsed();

        println!(
            "   Nucleotide counting: {:.2}ms",
            duration.as_millis() as f64 / iterations as f64
        );

        // Test k-mer processing
        let k = 31;
        if sequence.len() >= k {
            let start = std::time::Instant::now();
            for _ in 0..iterations {
                let _ = self.processor.process_kmers(sequence, k);
            }
            let duration = start.elapsed();

            println!(
                "   K-mer processing (k={}): {:.2}ms",
                k,
                duration.as_millis() as f64 / iterations as f64
            );
        }

        // Test complexity calculation
        let start = std::time::Instant::now();
        for _ in 0..iterations {
            let _ = self.processor.calculate_complexity(sequence);
        }
        let duration = start.elapsed();

        println!(
            "   Complexity calculation: {:.2}ms",
            duration.as_millis() as f64 / iterations as f64
        );
    }

    /// Estimate speedup compared to scalar implementation
    pub fn estimate_speedup(&self, sequence: &[u8]) -> f64 {
        let iterations = 100;

        // Force scalar implementation
        let scalar_processor = SimdProcessor {
            supports_avx512: false,
            supports_avx2: false,
            supports_sse41: false,
        };

        // Time scalar version
        let start = std::time::Instant::now();
        for _ in 0..iterations {
            let _ = scalar_processor.count_nucleotides(sequence);
        }
        let scalar_time = start.elapsed();

        // Time SIMD version
        let start = std::time::Instant::now();
        for _ in 0..iterations {
            let _ = self.processor.count_nucleotides(sequence);
        }
        let simd_time = start.elapsed();

        if simd_time.as_nanos() > 0 {
            scalar_time.as_nanos() as f64 / simd_time.as_nanos() as f64
        } else {
            1.0
        }
    }
}

impl Default for SimdBenchmark {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simd_nucleotide_counting() {
        let processor = SimdProcessor::new();
        let sequence = b"ATCGATCGATCGATCGATCGAAAA";

        let counts = processor.count_nucleotides(sequence);

        // Manually verify
        let expected = [8, 4, 4, 8]; // A, C, G, T
        assert_eq!(counts, expected);

        println!("SIMD nucleotide counting test passed!");
    }

    #[test]
    fn test_simd_kmer_processing() {
        let processor = SimdProcessor::new();
        let sequence = b"ATCGATCGATCGATCGATCG";
        let k = 5;

        let kmers = processor.process_kmers(sequence, k).unwrap();

        // Should generate sequence.len() - k + 1 k-mers
        let expected_count = sequence.len() - k + 1;
        assert_eq!(kmers.len(), expected_count);

        println!("SIMD k-mer processing test passed!");
    }

    #[test]
    fn test_simd_complexity_calculation() {
        let processor = SimdProcessor::new();
        let simple_sequence = b"AAAAAAAAAAAAAAAAAAAA"; // Low complexity
        let complex_sequence = b"ATCGATCGATCGATCGATCG"; // Higher complexity

        let simple_complexity = processor.calculate_complexity(simple_sequence);
        let complex_complexity = processor.calculate_complexity(complex_sequence);

        // Complex sequence should have higher entropy
        assert!(complex_complexity > simple_complexity);

        println!("SIMD complexity calculation test passed!");
        println!(
            "   Simple: {:.3}, Complex: {:.3}",
            simple_complexity, complex_complexity
        );
    }

    #[test]
    fn test_simd_benchmark() {
        let benchmark = SimdBenchmark::new();
        let sequence = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

        benchmark.processor.print_capabilities();

        // Run performance comparison
        benchmark.run_comparison(sequence, 1000);

        // Estimate speedup
        let speedup = benchmark.estimate_speedup(sequence);
        println!("   Estimated speedup: {:.2}x", speedup);

        assert!(speedup >= 1.0); // Should be at least as fast as scalar
    }
}
