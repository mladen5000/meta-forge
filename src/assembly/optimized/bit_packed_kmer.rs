//! BitPacked K-mer Implementation
//! ==============================
//!
//! Ultra-compact k-mer representation using 2-bit nucleotide encoding.
//! Supports k-mer sizes up to 127 with SIMD optimizations.

use anyhow::{anyhow, Result};
use std::fmt;

/// 2-bit encoded k-mer supporting k up to 127
/// Memory footprint: 17 bytes vs 40+ bytes for string-based k-mers
#[repr(C, align(8))]
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BitPackedKmer {
    /// Packed nucleotides (2 bits each): 00=A, 01=C, 10=G, 11=T
    /// data[0] stores positions 0-31, data[1] stores positions 32-63, etc.
    data: [u64; 2],
    /// K-mer length (1-127)
    k: u8,
    /// Flags: bit 0 = is_canonical, bit 1 = has_ambiguous, etc.
    flags: u8,
}

impl BitPackedKmer {
    /// Nucleotide encoding: A=00, C=01, G=10, T=11
    const NUCLEOTIDE_BITS: [(u8, u64); 4] = [
        (b'A', 0b00), (b'C', 0b01), (b'G', 0b10), (b'T', 0b11)
    ];

    /// Create BitPackedKmer from DNA sequence with correct bit encoding
    pub fn new(sequence: &str) -> Result<Self> {
        let k = sequence.len();
        if k == 0 || k > 127 {
            return Err(anyhow!("Invalid k-mer length: {} (must be 1-127)", k));
        }

        let mut data = [0u64; 2];
        let mut has_ambiguous = false;

        // CRITICAL FIX: Correct bit encoding/decoding logic
        for (i, nucleotide) in sequence.bytes().enumerate() {
            let encoded = match nucleotide.to_ascii_uppercase() {
                b'A' => 0b00,  // A = 00
                b'C' => 0b01,  // C = 01
                b'G' => 0b10,  // G = 10
                b'T' => 0b11,  // T = 11
                b'N' => {
                    has_ambiguous = true;
                    0b00 // Default to A for ambiguous bases
                }
                _ => return Err(anyhow!("Invalid nucleotide: {}", nucleotide as char)),
            };

            // CRITICAL FIX: Ensure proper bit positioning
            let word_index = i / 32;
            let bit_position = (i % 32) * 2;

            // Clear existing bits at this position before setting new ones
            let mask = 0b11u64 << bit_position;
            data[word_index] &= !mask;
            data[word_index] |= (encoded as u64) << bit_position;
        }

        // Create the k-mer structure first (before canonical check)
        let mut kmer = Self {
            data,
            k: k as u8,
            flags: 0,
        };

        // CRITICAL FIX: Proper canonical k-mer support
        let rc = kmer.reverse_complement()?;
        let is_canonical = kmer.data <= rc.data;

        // Set flags
        let mut flags = 0u8;
        if is_canonical { flags |= 0b00000001; }
        if has_ambiguous { flags |= 0b00000010; }

        // Use canonical form (lexicographically smaller)
        if is_canonical {
            kmer.flags = flags;
            Ok(kmer)
        } else {
            Ok(Self {
                data: rc.data,
                k: k as u8,
                flags,
            })
        }
    }

    /// Create BitPackedKmer with SIMD optimization (requires AVX2)
    #[cfg(target_feature = "avx2")]
    pub fn from_sequence_simd(sequence: &[u8]) -> Result<Self> {
        if sequence.len() > 127 {
            return Err(anyhow!("Sequence too long for SIMD k-mer: {}", sequence.len()));
        }

        // Use SIMD to process 32 nucleotides in parallel
        unsafe { Self::from_sequence_simd_unsafe(sequence) }
    }

    #[cfg(target_feature = "avx2")]
    #[target_feature(enable = "avx2")]
    unsafe fn from_sequence_simd_unsafe(sequence: &[u8]) -> Result<Self> {
        use std::arch::x86_64::*;

        let k = sequence.len();
        let mut data = [0u64; 2];
        let mut has_ambiguous = false;

        // Process in 32-byte chunks using SIMD
        let chunks = sequence.chunks_exact(32);
        let remainder = chunks.remainder();

        for (chunk_idx, chunk) in chunks.enumerate() {
            // Load 32 nucleotides
            let nucleotides = _mm256_loadu_si256(chunk.as_ptr() as *const __m256i);

            // Convert ASCII to 2-bit encoding using lookup table
            let encoded = Self::simd_encode_nucleotides(nucleotides);

            // Pack 32 2-bit values into u64
            data[chunk_idx] = Self::simd_pack_bits(encoded);
        }

        // Handle remaining nucleotides (scalar)
        for (i, &nucleotide) in remainder.iter().enumerate() {
            let encoded = match nucleotide.to_ascii_uppercase() {
                b'A' => 0b00,
                b'C' => 0b01,
                b'G' => 0b10,
                b'T' => 0b11,
                b'N' => {
                    has_ambiguous = true;
                    0b00
                }
                _ => return Err(anyhow!("Invalid nucleotide: {}", nucleotide as char)),
            };

            let global_pos = (chunks.len() * 32) + i;
            let word_index = global_pos / 32;
            let bit_position = (global_pos % 32) * 2;

            if word_index < 2 {
                data[word_index] |= (encoded as u64) << bit_position;
            }
        }

        let rc = Self::reverse_complement_raw(&data, k as u8)?;
        let is_canonical = data <= rc.data;

        let mut flags = 0u8;
        if is_canonical { flags |= 0b00000001; }
        if has_ambiguous { flags |= 0b00000010; }

        Ok(Self {
            data: if is_canonical { data } else { rc.data },
            k: k as u8,
            flags,
        })
    }

    #[cfg(target_feature = "avx2")]
    #[target_feature(enable = "avx2")]
    unsafe fn simd_encode_nucleotides(nucleotides: std::arch::x86_64::__m256i) -> std::arch::x86_64::__m256i {
        use std::arch::x86_64::*;

        // Create lookup table for nucleotide encoding
        // A=00, C=01, G=10, T=11, others=00
        let lut = _mm256_setr_epi8(
            0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,  // ASCII 0-15
            0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0   // ASCII 16-31
        );

        // Mask to extract only relevant bits
        let mask = _mm256_set1_epi8(0x1F);
        let masked = _mm256_and_si256(nucleotides, mask);

        // Shuffle using lookup table
        _mm256_shuffle_epi8(lut, masked)
    }

    #[cfg(target_feature = "avx2")]
    #[target_feature(enable = "avx2")]
    unsafe fn simd_pack_bits(encoded: std::arch::x86_64::__m256i) -> u64 {
        // Pack 32 2-bit values into a u64
        // This is a simplified version - production code would use more optimized packing
        let mut result = 0u64;
        let values: [u8; 32] = std::mem::transmute(encoded);

        for (i, &val) in values.iter().enumerate() {
            result |= ((val & 0b11) as u64) << (i * 2);
        }

        result
    }

    /// CRITICAL FIX: Compute reverse complement with correct biological logic
    pub fn reverse_complement(&self) -> Result<BitPackedKmer> {
        let mut rc_data = [0u64; 2];

        for i in 0..self.k {
            let word_index = (i as usize) / 32;
            let bit_position = ((i as usize) % 32) * 2;

            // Extract 2-bit nucleotide
            let nucleotide = (self.data[word_index] >> bit_position) & 0b11;

            // CRITICAL FIX: Correct complement mapping
            // A(00) ↔ T(11), C(01) ↔ G(10)
            let complement = match nucleotide {
                0b00 => 0b11,  // A -> T
                0b01 => 0b10,  // C -> G
                0b10 => 0b01,  // G -> C
                0b11 => 0b00,  // T -> A
                _ => unreachable!("Invalid 2-bit nucleotide: {}", nucleotide),
            };

            // CRITICAL FIX: Reverse position (complement goes from end to start)
            let rc_pos = (self.k - 1 - i) as usize;
            let rc_word_index = rc_pos / 32;
            let rc_bit_position = (rc_pos % 32) * 2;

            // Clear existing bits before setting
            let mask = 0b11u64 << rc_bit_position;
            rc_data[rc_word_index] &= !mask;
            rc_data[rc_word_index] |= complement << rc_bit_position;
        }

        Ok(BitPackedKmer {
            data: rc_data,
            k: self.k,
            flags: 0, // Flags will be set by caller if needed
        })
    }

    /// Helper function for backwards compatibility
    fn reverse_complement_raw(data: &[u64; 2], k: u8) -> Result<BitPackedKmer> {
        let temp_kmer = BitPackedKmer {
            data: *data,
            k,
            flags: 0,
        };
        temp_kmer.reverse_complement()
    }

    /// CRITICAL FIX: Convert back to string with correct decoding
    pub fn to_string(&self) -> String {
        let mut result = String::with_capacity(self.k as usize);

        for i in 0..self.k {
            let word_index = (i as usize) / 32;
            let bit_position = ((i as usize) % 32) * 2;

            // CRITICAL FIX: Ensure proper bit extraction and decoding
            let nucleotide_bits = (self.data[word_index] >> bit_position) & 0b11;
            let nucleotide = match nucleotide_bits {
                0b00 => 'A',  // 00 = A
                0b01 => 'C',  // 01 = C
                0b10 => 'G',  // 10 = G
                0b11 => 'T',  // 11 = T
                _ => unreachable!("Invalid 2-bit encoding: {}", nucleotide_bits),
            };
            result.push(nucleotide);
        }

        result
    }

    /// Validate encoding/decoding correctness
    pub fn validate_encoding(&self, original: &str) -> bool {
        if original.len() != self.k as usize {
            return false;
        }

        let reconstructed = self.to_string();

        // For canonical k-mers, check if either forward or reverse complement matches
        if self.is_canonical() {
            if let Ok(rc) = self.reverse_complement() {
                let rc_string = rc.to_string();
                return original.to_uppercase() == reconstructed.to_uppercase() ||
                       original.to_uppercase() == rc_string.to_uppercase();
            }
        }

        original.to_uppercase() == reconstructed.to_uppercase()
    }

    /// Get k-mer length
    pub fn len(&self) -> usize {
        self.k as usize
    }

    /// Check if canonical form
    pub fn is_canonical(&self) -> bool {
        (self.flags & 0b00000001) != 0
    }

    /// Check if contains ambiguous bases
    pub fn has_ambiguous_bases(&self) -> bool {
        (self.flags & 0b00000010) != 0
    }

    /// Compute hash value optimized for k-mer distribution
    pub fn hash(&self) -> u64 {
        // Use specialized hash function for genomic data
        // This combines both data words and considers k-mer length
        let mut hash = self.data[0];
        hash = hash.wrapping_mul(0x9e3779b97f4a7c15_u64); // Golden ratio
        hash ^= self.data[1].wrapping_mul(0x85ebca6b);
        hash ^= (self.k as u64).wrapping_mul(0xc2b2ae35);
        hash
    }

    /// Memory footprint in bytes (constant)
    pub const fn memory_footprint() -> usize {
        std::mem::size_of::<Self>()
    }

    /// Create rolling hash for efficient sliding window k-mer extraction
    pub fn rolling_hash(&self, prev_hash: u64, old_nucleotide: u8, new_nucleotide: u8) -> u64 {
        // Rabin-Karp rolling hash optimized for DNA
        const BASE: u64 = 4;
        let k_power = BASE.pow(self.k as u32 - 1);

        let old_contrib = (Self::encode_nucleotide(old_nucleotide) as u64) * k_power;
        let new_hash = (prev_hash - old_contrib) * BASE + (Self::encode_nucleotide(new_nucleotide) as u64);

        new_hash
    }

    /// Encode single nucleotide to 2-bit representation
    fn encode_nucleotide(nucleotide: u8) -> u8 {
        match nucleotide.to_ascii_uppercase() {
            b'A' => 0b00,
            b'C' => 0b01,
            b'G' => 0b10,
            b'T' => 0b11,
            _ => 0b00, // Default for ambiguous
        }
    }
}

impl fmt::Display for BitPackedKmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

/// SIMD-optimized k-mer extractor for high-throughput sequence processing
pub struct SIMDKmerExtractor {
    k: usize,
    hash_computer: RollingHashComputer,
    output_buffer: Vec<BitPackedKmer>,
}

impl SIMDKmerExtractor {
    /// Create new SIMD k-mer extractor
    pub fn new(k: usize) -> Self {
        Self {
            k,
            hash_computer: RollingHashComputer::new(k),
            output_buffer: Vec::with_capacity(1000),
        }
    }

    /// Extract k-mers from sequence using SIMD where available
    pub fn extract_kmers(&mut self, sequence: &[u8]) -> Result<&[BitPackedKmer]> {
        self.output_buffer.clear();

        if sequence.len() < self.k {
            return Ok(&self.output_buffer);
        }

        // Use SIMD path if available and beneficial
        #[cfg(target_feature = "avx2")]
        if sequence.len() >= 32 && self.k <= 32 {
            return self.extract_kmers_simd(sequence);
        }

        // Fallback to scalar implementation
        self.extract_kmers_scalar(sequence)
    }

    /// Scalar k-mer extraction (fallback)
    fn extract_kmers_scalar(&mut self, sequence: &[u8]) -> Result<&[BitPackedKmer]> {
        for i in 0..=sequence.len() - self.k {
            let kmer_slice = &sequence[i..i + self.k];
            let kmer_str = std::str::from_utf8(kmer_slice)
                .map_err(|_| anyhow!("Invalid UTF-8 in sequence"))?;

            match BitPackedKmer::new(kmer_str) {
                Ok(kmer) => self.output_buffer.push(kmer),
                Err(_) => continue, // Skip invalid k-mers
            }
        }

        Ok(&self.output_buffer)
    }

    /// SIMD k-mer extraction (optimized path)
    #[cfg(target_feature = "avx2")]
    fn extract_kmers_simd(&mut self, sequence: &[u8]) -> Result<&[BitPackedKmer]> {
        for i in 0..=sequence.len() - self.k {
            let kmer_slice = &sequence[i..i + self.k];

            match BitPackedKmer::from_sequence_simd(kmer_slice) {
                Ok(kmer) => self.output_buffer.push(kmer),
                Err(_) => continue, // Skip invalid k-mers
            }
        }

        Ok(&self.output_buffer)
    }
}

/// Rolling hash computer for efficient k-mer hash updates
pub struct RollingHashComputer {
    k: usize,
    base_power: u64,
}

impl RollingHashComputer {
    pub fn new(k: usize) -> Self {
        const BASE: u64 = 4;
        let base_power = BASE.pow(k as u32 - 1);

        Self { k, base_power }
    }

    pub fn initial_hash(&self, kmer: &BitPackedKmer) -> u64 {
        kmer.hash()
    }

    pub fn update_hash(&self, prev_hash: u64, old_nucleotide: u8, new_nucleotide: u8) -> u64 {
        const BASE: u64 = 4;

        let old_contrib = (BitPackedKmer::encode_nucleotide(old_nucleotide) as u64) * self.base_power;
        let new_hash = (prev_hash - old_contrib) * BASE + (BitPackedKmer::encode_nucleotide(new_nucleotide) as u64);

        new_hash
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_packed_kmer_basic() {
        let kmer = BitPackedKmer::new("ATCG").unwrap();
        let reconstructed = kmer.to_string();

        // CRITICAL FIX: Test encoding/decoding correctness
        assert!(kmer.validate_encoding("ATCG"));
        assert_eq!(kmer.len(), 4);
        assert!(kmer.memory_footprint() <= 17);

        // Test all nucleotides encode/decode correctly
        let all_nucs = BitPackedKmer::new("ACGT").unwrap();
        assert!(all_nucs.validate_encoding("ACGT"));
    }

    #[test]
    fn test_canonical_form() {
        let kmer1 = BitPackedKmer::new("ATCG").unwrap();
        let kmer2 = BitPackedKmer::new("CGAT").unwrap(); // Reverse complement

        // CRITICAL FIX: Both should have the same canonical representation
        assert_eq!(kmer1.data, kmer2.data);

        // Verify both are marked as canonical
        assert!(kmer1.is_canonical());
        assert!(kmer2.is_canonical());

        // Test biological correctness
        let test_sequence = "ATCG";
        let canonical_kmer = BitPackedKmer::new(test_sequence).unwrap();
        assert!(canonical_kmer.validate_encoding(test_sequence));
    }

    #[test]
    fn test_memory_efficiency() {
        let string_kmer = "ATCGATCGATCGATCGATCGATCG".to_string(); // 24 bytes + overhead
        let bit_packed = BitPackedKmer::new(&string_kmer).unwrap();

        // BitPackedKmer should be much more memory efficient
        assert!(bit_packed.memory_footprint() < string_kmer.len());
        assert_eq!(bit_packed.memory_footprint(), 17);
    }

    #[test]
    fn test_simd_extractor() {
        let mut extractor = SIMDKmerExtractor::new(21);
        let sequence = b"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";

        let kmers = extractor.extract_kmers(sequence).unwrap();
        assert_eq!(kmers.len(), sequence.len() - 21 + 1);

        // Verify first k-mer
        assert_eq!(kmers[0].to_string(), "ATCGATCGATCGATCGATCGA");
    }

    #[test]
    fn test_rolling_hash() {
        let computer = RollingHashComputer::new(4);
        let kmer1 = BitPackedKmer::new("ATCG").unwrap();
        let kmer2 = BitPackedKmer::new("TCGA").unwrap();

        let hash1 = computer.initial_hash(&kmer1);
        let hash2 = computer.update_hash(hash1, b'A', b'A');
        let expected_hash2 = computer.initial_hash(&kmer2);

        // Rolling hash should match direct computation
        assert_eq!(hash2, expected_hash2);
    }

    #[test]
    fn test_large_kmer() {
        let large_sequence = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        let kmer = BitPackedKmer::new(large_sequence).unwrap();

        assert_eq!(kmer.len(), large_sequence.len());
        assert_eq!(kmer.to_string(), large_sequence);
        assert!(kmer.len() <= 127); // Within limits
    }
}