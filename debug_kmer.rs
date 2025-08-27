use std::hash::{Hash, Hasher};

fn main() {
    // Test ATCG vs CGAT packing
    let seq1 = "ATCG";
    let seq2 = "CGAT";
    
    println!("Testing {} vs {}", seq1, seq2);
    
    // Manual bit packing like in BitPackedKmer
    let mut packed1 = Vec::new();
    let mut packed2 = Vec::new();
    let mut word1 = 0u64;
    let mut word2 = 0u64;
    let mut used = 0;
    
    for c in seq1.chars() {
        let bits = match c {
            'A' => 0b00,
            'C' => 0b01,
            'G' => 0b10,
            'T' => 0b11,
            _ => panic!("Invalid nucleotide"),
        } as u64;
        word1 |= bits << (62 - used);
        used += 2;
    }
    packed1.push(word1);
    
    used = 0;
    for c in seq2.chars() {
        let bits = match c {
            'A' => 0b00,
            'C' => 0b01,
            'G' => 0b10,
            'T' => 0b11,
            _ => panic!("Invalid nucleotide"),
        } as u64;
        word2 |= bits << (62 - used);
        used += 2;
    }
    packed2.push(word2);
    
    println!("Packed1: {:064b} ({})", word1, word1);
    println!("Packed2: {:064b} ({})", word2, word2);
    
    // Hash them
    let mut h1 = std::collections::hash_map::DefaultHasher::new();
    packed1.hash(&mut h1);
    4usize.hash(&mut h1);
    let hash1 = h1.finish();
    
    let mut h2 = std::collections::hash_map::DefaultHasher::new();
    packed2.hash(&mut h2);
    4usize.hash(&mut h2);
    let hash2 = h2.finish();
    
    println!("Hash1: {}", hash1);
    println!("Hash2: {}", hash2);
    println!("Equal: {}", hash1 == hash2);
}
