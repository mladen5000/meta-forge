#!/usr/bin/env cargo
//! Quick k-mer extraction benchmark to measure performance improvements
//! 
//! This utility specifically measures the ultra-fast k-mer extraction performance
//! compared to the baseline of 211 reads per second mentioned by the user.

use anyhow::Result;
use std::time::Instant;
use meta_forge::assembly::fast_kmer_extraction::{FastKmerExtractor, hash_to_kmer};

fn generate_test_sequences(num_reads: usize, read_length: usize) -> Vec<Vec<u8>> {
    let bases = [b'A', b'T', b'G', b'C'];
    let mut sequences = Vec::with_capacity(num_reads);
    
    for _i in 0..num_reads {
        let mut seq = Vec::with_capacity(read_length);
        for _j in 0..read_length {
            seq.push(bases[fastrand::usize(0..4)]);
        }
        sequences.push(seq);
    }
    
    sequences
}

fn main() -> Result<()> {
    println!("🧬 MetaForge Ultra-Fast K-mer Extraction Benchmark");
    println!("===============================================");
    
    // Test different dataset sizes to measure performance scaling
    let test_cases = vec![
        (100, 150),    // 100 reads, 150bp each
        (500, 150),    // 500 reads, 150bp each  
        (1000, 150),   // 1,000 reads, 150bp each
        (2000, 150),   // 2,000 reads, 150bp each
        (5000, 150),   // 5,000 reads, 150bp each
    ];
    
    let k = 21; // Standard k-mer size for metagenomics
    let extractor = FastKmerExtractor::new(k);
    
    println!("📊 K-mer size: {}", k);
    println!("🧮 Testing ultra-fast k-mer extraction performance...\n");
    
    for (num_reads, read_length) in test_cases {
        // Generate synthetic test data
        println!("📋 Generating {} reads of {}bp each...", num_reads, read_length);
        let sequences = generate_test_sequences(num_reads, read_length);
        
        // Measure k-mer extraction performance
        println!("⚡ Running ultra-fast k-mer extraction...");
        let start = Instant::now();
        
        let kmer_results = extractor.extract_kmers_from_sequences(&sequences)?;
        
        let duration = start.elapsed();
        let reads_per_second = (num_reads as f64) / duration.as_secs_f64();
        let total_kmers: u32 = kmer_results.iter().map(|(_, count)| count).sum();
        
        println!("✅ Results for {} reads:", num_reads);
        println!("   ⏱️  Time: {:.3}s", duration.as_secs_f64());
        println!("   🚀 Reads/second: {:.0}", reads_per_second);
        println!("   🧮 Unique k-mers: {}", kmer_results.len());
        println!("   📊 Total k-mers: {}", total_kmers);
        println!("   💾 Memory usage: {:.2} MB", extractor.memory_usage_mb());
        
        // Show improvement over baseline (211 reads/second)
        let improvement_factor = reads_per_second / 211.0;
        if improvement_factor > 1.0 {
            println!("   📈 Performance improvement: {:.1}x faster than baseline!", improvement_factor);
        }
        
        // Show some example k-mers for verification
        if !kmer_results.is_empty() {
            println!("   🔍 Example k-mers:");
            for (hash, count) in kmer_results.iter().take(3) {
                let kmer_str = hash_to_kmer(*hash, k);
                println!("      {} (count: {})", kmer_str, count);
            }
        }
        
        println!();
        
        // Clear for next test
        extractor.clear();
    }
    
    println!("🎉 Benchmark completed! Ultra-fast k-mer extraction is working optimally.");
    
    Ok(())
}