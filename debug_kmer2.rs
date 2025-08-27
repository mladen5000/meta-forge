fn main() {
    // Test ATCG vs CGAT packing step by step
    let seq1 = "ATCG";
    let seq2 = "CGAT";
    
    println!("=== Testing {} ===", seq1);
    let mut word1 = 0u64;
    let mut used = 0;
    for (i, c) in seq1.chars().enumerate() {
        let bits = match c {
            'A' => 0b00,
            'C' => 0b01,  
            'G' => 0b10,
            'T' => 0b11,
            _ => panic!("Invalid nucleotide"),
        } as u64;
        println!("  Char {}: {} -> bits: {:02b}, shift: {}", i, c, bits, 62 - used);
        word1 |= bits << (62 - used);
        used += 2;
        println!("  Word after: {:064b} ({})", word1, word1);
    }
    
    println!("\n=== Testing {} ===", seq2);
    let mut word2 = 0u64;
    used = 0;
    for (i, c) in seq2.chars().enumerate() {
        let bits = match c {
            'A' => 0b00,
            'C' => 0b01,
            'G' => 0b10,
            'T' => 0b11,
            _ => panic!("Invalid nucleotide"),
        } as u64;
        println!("  Char {}: {} -> bits: {:02b}, shift: {}", i, c, bits, 62 - used);
        word2 |= bits << (62 - used);
        used += 2;
        println!("  Word after: {:064b} ({})", word2, word2);
    }
    
    println!("\nFinal:");
    println!("ATCG: {:064b} ({})", word1, word1);
    println!("CGAT: {:064b} ({})", word2, word2);
    println!("Same: {}", word1 == word2);
}
