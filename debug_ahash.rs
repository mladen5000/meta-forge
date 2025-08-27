// Test with ahash to match the real implementation
use std::hash::{Hash, Hasher};

fn main() {
    let data1 = vec![3891110078048108544u64];
    let data2 = vec![7133701809754865664u64];
    
    println!("Testing with ahash:");
    
    // Test with ahash
    let mut h1 = ahash::AHasher::default();
    data1.hash(&mut h1);
    4usize.hash(&mut h1);
    let hash1 = h1.finish();
    
    let mut h2 = ahash::AHasher::default();
    data2.hash(&mut h2);
    4usize.hash(&mut h2);
    let hash2 = h2.finish();
    
    println!("data1 hash: {}", hash1);
    println!("data2 hash: {}", hash2);
    println!("Same: {}", hash1 == hash2);
    
    // Test if ahash is deterministic
    let mut h3 = ahash::AHasher::default();
    data1.hash(&mut h3);
    4usize.hash(&mut h3);
    let hash3 = h3.finish();
    
    println!("data1 hash again: {}", hash3);
    println!("Deterministic: {}", hash1 == hash3);
    
    // Now test the exact same data (which shouldn't happen)
    let mut h4 = ahash::AHasher::default();
    vec![3891110078048108544u64].hash(&mut h4);
    4usize.hash(&mut h4);
    let hash4 = h4.finish();
    
    let mut h5 = ahash::AHasher::default();
    vec![3891110078048108544u64].hash(&mut h5);
    4usize.hash(&mut h5);
    let hash5 = h5.finish();
    
    println!("Same data hashes: {} vs {}", hash4, hash5);
    println!("Same data equal: {}", hash4 == hash5);
}
