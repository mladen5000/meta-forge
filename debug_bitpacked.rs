// Copy of the BitPackedKmer implementation for testing
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BitPackedKmer {
    pub packed_data: Vec<u64>,
    pub k: usize,
    pub hash: u64,
}

impl BitPackedKmer {
    pub fn new(seq: &str) -> Result<Self, String> {
        if seq.is_empty() {
            return Err("Empty k-mer".to_string());
        }
        let k = seq.len();
        if k > 1024 {
            return Err("k-mer too long (>1024)".to_string());
        }
        let mut packed = Vec::with_capacity(k.div_ceil(32));
        let mut word = 0u64;
        let mut used = 0;
        for c in seq.chars() {
            let bits = match c.to_ascii_uppercase() {
                'A' => 0b00,
                'C' => 0b01,
                'G' => 0b10,
                'T' => 0b11,
                x => return Err(format!("Invalid nt {}", x)),
            } as u64;
            word |= bits << (62 - used);
            used += 2;
            if used >= 64 {
                packed.push(word);
                word = 0;
                used = 0;
            }
        }
        if used > 0 {
            packed.push(word);
        }
        
        let mut h = std::collections::hash_map::DefaultHasher::new();
        packed.hash(&mut h);
        k.hash(&mut h);
        let hash = h.finish();
        
        Ok(Self {
            packed_data: packed,
            k,
            hash,
        })
    }
}

fn main() {
    let k1 = BitPackedKmer::new("ATCG").unwrap();
    let k2 = BitPackedKmer::new("CGAT").unwrap();
    
    println!("ATCG: packed={:?}, hash={}", k1.packed_data, k1.hash);
    println!("CGAT: packed={:?}, hash={}", k2.packed_data, k2.hash);
    println!("Same packed: {}", k1.packed_data == k2.packed_data);
    println!("Same hash: {}", k1.hash == k2.hash);
}
