//! Comprehensive Assembly Accuracy Test Suite
//! ==========================================
//!
//! Tests for assembly quality metrics including N50, N90, coverage distribution,
//! contiguity, and accuracy validation against known reference sequences.
//! This suite ensures the laptop assembly optimization maintains high quality.

use meta_forge::assembly::{LaptopAssembler, LaptopConfig, LaptopAssemblyGraph};
use meta_forge::core::data_structures::{CorrectedRead, AssemblyStats, Contig, ContigType, CorrectionMetadata};
use meta_forge::assembly::adaptive_k::{AdaptiveKSelector, AdaptiveKConfig};
use anyhow::Result;
use std::collections::HashMap;
use std::time::Instant;

/// Test data generator for reproducible assembly validation
#[derive(Debug, Clone)]
pub struct TestDataGenerator {
    reference_sequence: String,
    read_length: usize,
    coverage_depth: f64,
    error_rate: f64,
}

impl TestDataGenerator {
    pub fn new(reference: &str, read_len: usize, coverage: f64, error_rate: f64) -> Self {
        Self {
            reference_sequence: reference.to_string(),
            read_length: read_len,
            coverage_depth: coverage,
            error_rate,
        }
    }

    /// Generate synthetic reads with controlled parameters
    pub fn generate_reads(&self) -> Vec<CorrectedRead> {
        let num_reads = ((self.reference_sequence.len() as f64 * self.coverage_depth) / self.read_length as f64) as usize;
        let mut reads = Vec::new();
        
        for i in 0..num_reads {
            let start_pos = fastrand::usize(0..=(self.reference_sequence.len().saturating_sub(self.read_length)));
            let end_pos = (start_pos + self.read_length).min(self.reference_sequence.len());
            
            if end_pos <= start_pos {
                continue;
            }
            
            let mut sequence = self.reference_sequence[start_pos..end_pos].to_string();
            
            // Introduce random errors based on error rate
            if self.error_rate > 0.0 {
                sequence = self.introduce_errors(&sequence);
            }
            
            reads.push(CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; self.read_length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test_generator".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
            });
        }
        
        reads
    }
    
    fn introduce_errors(&self, sequence: &str) -> String {
        let mut result = String::new();
        let bases = ['A', 'C', 'G', 'T'];
        
        for c in sequence.chars() {
            if fastrand::f64() < self.error_rate {
                // Introduce random substitution
                let random_base = bases[fastrand::usize(0..4)];
                result.push(random_base);
            } else {
                result.push(c);
            }
        }
        
        result
    }
}

/// Assembly quality metrics calculator
#[derive(Debug, Clone)]
pub struct AssemblyQualityMetrics {
    pub n50: usize,
    pub n90: usize,
    pub total_length: usize,
    pub num_contigs: usize,
    pub largest_contig: usize,
    pub coverage_mean: f64,
    pub coverage_std: f64,
    pub gc_content: f64,
    pub gaps_count: usize,
    pub contiguity_score: f64,
}

impl AssemblyQualityMetrics {
    pub fn calculate(contigs: &[Contig]) -> Self {
        if contigs.is_empty() {
            return Self::default();
        }
        
        let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
        lengths.sort_by(|a, b| b.cmp(a)); // Sort descending
        
        let total_length: usize = lengths.iter().sum();
        let half_length = total_length / 2;
        let ninety_percent_length = (total_length as f64 * 0.9) as usize;
        
        // Calculate N50
        let mut cumulative_length = 0;
        let mut n50 = 0;
        for &length in &lengths {
            cumulative_length += length;
            if cumulative_length >= half_length {
                n50 = length;
                break;
            }
        }
        
        // Calculate N90
        cumulative_length = 0;
        let mut n90 = 0;
        for &length in &lengths {
            cumulative_length += length;
            if cumulative_length >= ninety_percent_length {
                n90 = length;
                break;
            }
        }
        
        // Coverage statistics
        let coverages: Vec<f64> = contigs.iter().map(|c| c.coverage).collect();
        let coverage_mean = coverages.iter().sum::<f64>() / coverages.len() as f64;
        let coverage_variance = coverages.iter()
            .map(|c| (c - coverage_mean).powi(2))
            .sum::<f64>() / coverages.len() as f64;
        let coverage_std = coverage_variance.sqrt();
        
        // GC content calculation
        let total_gc = contigs.iter()
            .map(|c| c.sequence.chars().filter(|&ch| ch == 'G' || ch == 'C').count())
            .sum::<usize>();
        let total_bases = contigs.iter().map(|c| c.sequence.len()).sum::<usize>();
        let gc_content = if total_bases > 0 { total_gc as f64 / total_bases as f64 } else { 0.0 };
        
        // Gaps count (N characters)
        let gaps_count = contigs.iter()
            .map(|c| c.sequence.chars().filter(|&ch| ch == 'N').count())
            .sum::<usize>();
        
        // Contiguity score (weighted by contig length)
        let contiguity_score = if lengths.len() > 1 {
            let largest = lengths[0] as f64;
            let weighted_sum: f64 = lengths.iter()
                .map(|&len| (len as f64).powi(2))
                .sum();
            weighted_sum / (total_length as f64 * largest)
        } else {
            1.0
        };
        
        Self {
            n50,
            n90,
            total_length,
            num_contigs: contigs.len(),
            largest_contig: lengths.first().copied().unwrap_or(0),
            coverage_mean,
            coverage_std,
            gc_content,
            gaps_count,
            contiguity_score,
        }
    }
}

impl Default for AssemblyQualityMetrics {
    fn default() -> Self {
        Self {
            n50: 0,
            n90: 0,
            total_length: 0,
            num_contigs: 0,
            largest_contig: 0,
            coverage_mean: 0.0,
            coverage_std: 0.0,
            gc_content: 0.0,
            gaps_count: 0,
            contiguity_score: 0.0,
        }
    }
}

/// Test suite for assembly accuracy validation
#[test]
fn test_assembly_accuracy_synthetic_data() {
    // Test with known reference sequence
    let reference = "ATCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGT";
    let generator = TestDataGenerator::new(reference, 50, 10.0, 0.01); // 10x coverage, 1% error
    let reads = generator.generate_reads();
    
    println!("Generated {} reads for {}bp reference", reads.len(), reference.len());
    
    let config = LaptopConfig::medium_memory();
    let assembler = LaptopAssembler::new(config);
    
    let start_time = Instant::now();
    let contigs = assembler.assemble(&reads).unwrap();
    let assembly_time = start_time.elapsed();
    
    let metrics = AssemblyQualityMetrics::calculate(&contigs);
    
    println!("Assembly Accuracy Results:");
    println!("  Time: {:?}", assembly_time);
    println!("  Contigs: {}", metrics.num_contigs);
    println!("  Total length: {}bp", metrics.total_length);
    println!("  N50: {}bp", metrics.n50);
    println!("  N90: {}bp", metrics.n90);
    println!("  Largest contig: {}bp", metrics.largest_contig);
    println!("  Coverage mean: {:.2}", metrics.coverage_mean);
    println!("  Coverage std: {:.2}", metrics.coverage_std);
    println!("  GC content: {:.2}%", metrics.gc_content * 100.0);
    println!("  Contiguity score: {:.3}", metrics.contiguity_score);
    
    // Quality assertions
    assert!(!contigs.is_empty(), "Assembly should produce contigs");
    assert!(metrics.total_length >= reference.len() / 2, "Assembly should recover at least 50% of reference length");
    assert!(metrics.n50 > 20, "N50 should be reasonable for synthetic data");
    assert!(metrics.coverage_mean > 1.0, "Mean coverage should be positive");
    assert!(metrics.gc_content >= 0.0 && metrics.gc_content <= 1.0, "GC content should be valid percentage");
}

#[test]
fn test_assembly_accuracy_high_coverage() {
    // Test with high coverage to validate assembly completeness
    let reference = "ATGCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTAG";
    let generator = TestDataGenerator::new(reference, 40, 25.0, 0.005); // 25x coverage, 0.5% error
    let reads = generator.generate_reads();
    
    let config = LaptopConfig::high_memory();
    let assembler = LaptopAssembler::new(config);
    
    let contigs = assembler.assemble(&reads).unwrap();
    let metrics = AssemblyQualityMetrics::calculate(&contigs);
    
    println!("High Coverage Assembly:");
    println!("  Reference: {}bp", reference.len());
    println!("  Assembled: {}bp", metrics.total_length);
    println!("  Recovery rate: {:.1}%", (metrics.total_length as f64 / reference.len() as f64) * 100.0);
    println!("  N50: {}bp", metrics.n50);
    
    // High coverage should improve assembly quality
    assert!(metrics.total_length >= reference.len() * 8 / 10, "High coverage should recover ≥80% of reference");
    assert!(metrics.n50 >= 30, "High coverage should improve N50");
    assert!(metrics.coverage_mean >= 5.0, "High coverage should maintain good coverage depth");
}

#[test]
fn test_assembly_accuracy_low_coverage() {
    // Test with low coverage to validate graceful degradation
    let reference = "ATCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGTCGATCGATCGAAGT";
    let generator = TestDataGenerator::new(reference, 30, 3.0, 0.02); // 3x coverage, 2% error
    let reads = generator.generate_reads();
    
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);
    
    let contigs = assembler.assemble(&reads).unwrap();
    let metrics = AssemblyQualityMetrics::calculate(&contigs);
    
    println!("Low Coverage Assembly:");
    println!("  Contigs: {}", metrics.num_contigs);
    println!("  Recovery rate: {:.1}%", (metrics.total_length as f64 / reference.len() as f64) * 100.0);
    
    // Low coverage should still produce reasonable results
    assert!(!contigs.is_empty(), "Low coverage should still produce some contigs");
    assert!(metrics.total_length >= reference.len() / 4, "Should recover at least 25% with low coverage");
    assert!(metrics.n50 > 0, "N50 should be positive even with low coverage");
}

#[test]
fn test_assembly_accuracy_repetitive_regions() {
    // Test with repetitive sequences that challenge assembly
    let repeat_unit = "ATCGATCG";
    let reference = format!(
        "GCTAGCTA{}{}{}{}TACGTACG", 
        repeat_unit, repeat_unit, repeat_unit, repeat_unit
    );
    
    let generator = TestDataGenerator::new(&reference, 25, 15.0, 0.01);
    let reads = generator.generate_reads();
    
    let config = LaptopConfig::medium_memory();
    let assembler = LaptopAssembler::new(config);
    
    let contigs = assembler.assemble(&reads).unwrap();
    let metrics = AssemblyQualityMetrics::calculate(&contigs);
    
    println!("Repetitive Region Assembly:");
    println!("  Contigs: {}", metrics.num_contigs);
    println!("  Largest contig: {}bp", metrics.largest_contig);
    println!("  Contiguity score: {:.3}", metrics.contiguity_score);
    
    // Repetitive regions may fragment assembly
    assert!(!contigs.is_empty(), "Should handle repetitive regions");
    assert!(metrics.contiguity_score >= 0.0, "Contiguity score should be valid");
}

#[test]
fn test_assembly_n50_n90_calculation() {
    // Test N50/N90 calculations with known contig lengths
    let test_contigs = vec![
        create_test_contig(0, "A".repeat(100), 5.0),
        create_test_contig(1, "T".repeat(80), 4.0),
        create_test_contig(2, "G".repeat(60), 3.0),
        create_test_contig(3, "C".repeat(40), 2.0),
        create_test_contig(4, "A".repeat(20), 1.0),
    ];
    
    let metrics = AssemblyQualityMetrics::calculate(&test_contigs);
    
    // Total: 300bp, sorted lengths: [100, 80, 60, 40, 20]
    // Cumulative: [100, 180, 240, 280, 300]
    // N50 (≥150bp): 80bp
    // N90 (≥270bp): 20bp
    
    assert_eq!(metrics.total_length, 300);
    assert_eq!(metrics.n50, 80);
    assert_eq!(metrics.n90, 20);
    assert_eq!(metrics.largest_contig, 100);
    assert_eq!(metrics.num_contigs, 5);
    
    println!("N50/N90 Validation: N50={}, N90={}", metrics.n50, metrics.n90);
}

#[test]
fn test_assembly_accuracy_edge_cases() {
    let config = LaptopConfig::low_memory();
    let assembler = LaptopAssembler::new(config);
    
    // Test empty reads
    let empty_result = assembler.assemble(&vec![]);
    assert!(empty_result.is_ok(), "Should handle empty input gracefully");
    
    // Test single read
    let single_read = vec![create_test_read(0, "ATCGATCGATCG")];
    let single_result = assembler.assemble(&single_read).unwrap();
    let single_metrics = AssemblyQualityMetrics::calculate(&single_result);
    
    println!("Single read assembly: {} contigs", single_metrics.num_contigs);
    
    // Test identical reads
    let identical_reads = vec![
        create_test_read(0, "ATCGATCGATCG"),
        create_test_read(1, "ATCGATCGATCG"),
        create_test_read(2, "ATCGATCGATCG"),
    ];
    let identical_result = assembler.assemble(&identical_reads).unwrap();
    let identical_metrics = AssemblyQualityMetrics::calculate(&identical_result);
    
    println!("Identical reads assembly: coverage_mean={:.1}", identical_metrics.coverage_mean);
    assert!(identical_metrics.coverage_mean >= 1.0, "Identical reads should increase coverage");
}

/// Helper function to create test reads
fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
    CorrectedRead {
        id,
        original: sequence.to_string(),
        corrected: sequence.to_string(),
        corrections: Vec::new(),
        quality_scores: vec![30; sequence.len()],
        correction_metadata: CorrectionMetadata {
            algorithm: "test".to_string(),
            confidence_threshold: 0.95,
            context_window: 5,
            correction_time_ms: 0,
        },
    }
}

/// Helper function to create test contigs
fn create_test_contig(id: usize, sequence: String, coverage: f64) -> Contig {
    Contig {
        id,
        sequence: sequence.clone(),
        length: sequence.len(),
        coverage,
        contig_type: ContigType::Primary,
    }
}
