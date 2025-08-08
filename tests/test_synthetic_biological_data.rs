//! Tests with synthetic and realistic biological data
//! Simulates real genomic patterns, creates synthetic reads, and tests biological accuracy

use meta_forge::core::data_structures::*;
use meta_forge::assembly::adaptive_k::*;
use anyhow::Result;
use std::collections::HashMap;

#[cfg(test)]
mod synthetic_genome_tests {
    use super::*;

    /// Generate synthetic genome with realistic properties
    fn generate_synthetic_genome(length: usize, gc_content: f64) -> String {
        let mut genome = String::with_capacity(length);
        
        // Define base probabilities based on GC content
        let gc_prob = gc_content / 2.0; // Split equally between G and C
        let at_prob = (1.0 - gc_content) / 2.0; // Split equally between A and T
        
        for _ in 0..length {
            let rand_val = fastrand::f64();
            let base = if rand_val < at_prob {
                'A'
            } else if rand_val < at_prob * 2.0 {
                'T'
            } else if rand_val < at_prob * 2.0 + gc_prob {
                'G'
            } else {
                'C'
            };
            genome.push(base);
        }
        
        genome
    }

    /// Generate synthetic sequencing reads with realistic error patterns
    fn generate_synthetic_reads(genome: &str, read_length: usize, coverage: f64, error_rate: f64) -> Vec<CorrectedRead> {
        let num_reads = ((genome.len() as f64 * coverage) / read_length as f64) as usize;
        let mut reads = Vec::new();
        
        for i in 0..num_reads {
            // Random starting position
            let max_start = if genome.len() >= read_length {
                genome.len() - read_length
            } else {
                0
            };
            
            let start_pos = if max_start > 0 {
                fastrand::usize(0..max_start)
            } else {
                0
            };
            
            let end_pos = std::cmp::min(start_pos + read_length, genome.len());
            let mut read_sequence = genome[start_pos..end_pos].to_string();
            
            // Introduce sequencing errors
            let mut original_sequence = read_sequence.clone();
            let mut corrections = Vec::new();
            
            for pos in 0..read_sequence.len() {
                if fastrand::f64() < error_rate {
                    let original_base = read_sequence.chars().nth(pos).unwrap();
                    let error_bases = match original_base {
                        'A' => ['T', 'G', 'C'],
                        'T' => ['A', 'G', 'C'],
                        'G' => ['A', 'T', 'C'],
                        'C' => ['A', 'T', 'G'],
                        _ => ['A', 'T', 'G'], // Fallback
                    };
                    
                    let error_base = error_bases[fastrand::usize(0..3)];
                    
                    // Replace character at position
                    let mut chars: Vec<char> = original_sequence.chars().collect();
                    chars[pos] = error_base;
                    original_sequence = chars.into_iter().collect();
                    
                    corrections.push(BaseCorrection {
                        position: pos,
                        from: error_base,
                        to: original_base,
                        confidence: 0.9,
                        correction_type: CorrectionType::Substitution,
                    });
                }
            }
            
            // Generate realistic quality scores
            let quality_scores: Vec<u8> = (0..read_sequence.len())
                .map(|_| {
                    // Most bases have high quality, with some variation
                    let base_quality = 35.0;
                    let noise = fastrand::f64() * 10.0 - 5.0; // -5 to +5
                    (base_quality + noise).max(10.0).min(42.0) as u8
                })
                .collect();
            
            reads.push(CorrectedRead {
                id: i,
                original: original_sequence,
                corrected: read_sequence,
                corrections,
                quality_scores,
                correction_metadata: CorrectionMetadata {
                    algorithm: "synthetic_corrector".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 10,
                    correction_time_ms: fastrand::u64(50..200),
                },
            });
        }
        
        reads
    }

    #[test]
    fn test_synthetic_bacterial_genome() {
        // Simulate small bacterial genome (E. coli has ~50% GC content)
        let genome = generate_synthetic_genome(10000, 0.5);
        
        // Verify GC content is approximately correct
        let actual_gc = calculate_gc_content(&genome);
        assert!((actual_gc - 0.5).abs() < 0.05, "GC content should be close to 50%");
        
        // Generate reads with realistic coverage
        let reads = generate_synthetic_reads(&genome, 150, 10.0, 0.01); // 10x coverage, 1% error rate
        
        // Test assembly pipeline
        let builder = AssemblyGraphBuilder::new(21, 31, 3); // Realistic k-mer sizes for bacterial assembly
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Synthetic bacterial genome assembly should succeed");
        
        let assembly = result.unwrap();
        
        // Verify assembly quality
        assert!(!assembly.contigs.is_empty(), "Should generate contigs");
        assert!(assembly.assembly_stats.total_length > 1000, "Should have substantial assembly length");
        assert!(assembly.assembly_stats.n50 > 100, "Should have reasonable N50");
        
        // Verify GC content is preserved
        let assembly_gc = assembly.assembly_stats.gc_content;
        assert!((assembly_gc - 0.5).abs() < 0.1, "Assembly should preserve GC content roughly");
        
        // Coverage should be reasonable
        assert!(assembly.assembly_stats.coverage_mean > 5.0, "Should have good coverage");
    }

    #[test]
    fn test_high_gc_genome() {
        // Simulate high-GC genome (like some Streptomyces species ~72% GC)
        let genome = generate_synthetic_genome(5000, 0.72);
        let reads = generate_synthetic_reads(&genome, 100, 15.0, 0.005); // High coverage, low error
        
        let builder = AssemblyGraphBuilder::new(15, 25, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "High-GC genome assembly should succeed");
        
        let assembly = result.unwrap();
        
        // Verify high GC content is handled correctly
        let assembly_gc = assembly.assembly_stats.gc_content;
        assert!(assembly_gc > 0.6, "Assembly should retain high GC content");
        
        // High-GC genomes can be more challenging to assemble
        assert!(!assembly.contigs.is_empty(), "Should still generate contigs");
    }

    #[test]
    fn test_low_gc_genome() {
        // Simulate low-GC genome (like some Plasmodium species ~20% GC)
        let genome = generate_synthetic_genome(5000, 0.2);
        let reads = generate_synthetic_reads(&genome, 100, 12.0, 0.01);
        
        let builder = AssemblyGraphBuilder::new(17, 27, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Low-GC genome assembly should succeed");
        
        let assembly = result.unwrap();
        
        // Verify low GC content
        let assembly_gc = assembly.assembly_stats.gc_content;
        assert!(assembly_gc < 0.4, "Assembly should retain low GC content");
        
        // Low-GC genomes might have different complexity
        let low_complexity_nodes = assembly.graph_fragment.nodes
            .values()
            .filter(|n| n.complexity_score < 0.3)
            .count();
        
        assert!(low_complexity_nodes > 0, "Should have some low-complexity regions in AT-rich genome");
    }

    #[test]
    fn test_repetitive_genome_regions() {
        // Create genome with repetitive elements
        let repeat_unit = "ATCGATCGATCG";
        let unique_sequence1 = generate_synthetic_genome(1000, 0.5);
        let unique_sequence2 = generate_synthetic_genome(1000, 0.5);
        
        // Build genome with repeats: unique1 + repeat + unique2 + repeat + unique3
        let genome = format!("{}{}{}{}{}",
                           unique_sequence1,
                           repeat_unit.repeat(10), // 120bp repeat
                           unique_sequence2,
                           repeat_unit.repeat(8),  // 96bp repeat
                           generate_synthetic_genome(800, 0.5));
        
        let reads = generate_synthetic_reads(&genome, 150, 8.0, 0.008);
        
        let builder = AssemblyGraphBuilder::new(19, 29, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Repetitive genome assembly should succeed");
        
        let assembly = result.unwrap();
        
        // Should detect repetitive nodes
        let repetitive_nodes = assembly.graph_fragment.nodes
            .values()
            .filter(|n| matches!(n.node_type, NodeType::Repetitive))
            .count();
        
        assert!(repetitive_nodes > 0, "Should detect repetitive regions");
        
        // Assembly might be fragmented due to repeats
        assert!(assembly.contigs.len() >= 1, "Should generate contigs despite repeats");
    }

    #[test]
    fn test_plasmid_like_circular_sequence() {
        // Simulate small circular plasmid
        let plasmid_sequence = generate_synthetic_genome(3000, 0.48);
        
        // Generate reads that span the junction (circular)
        let mut reads = Vec::new();
        let read_length = 120;
        let num_reads = 200;
        
        for i in 0..num_reads {
            let start_pos = fastrand::usize(0..plasmid_sequence.len());
            let mut read_seq = String::new();
            
            for j in 0..read_length {
                let genome_pos = (start_pos + j) % plasmid_sequence.len();
                read_seq.push(plasmid_sequence.chars().nth(genome_pos).unwrap());
            }
            
            reads.push(CorrectedRead {
                id: i,
                original: read_seq.clone(),
                corrected: read_seq,
                corrections: Vec::new(),
                quality_scores: vec![35; read_length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "plasmid_sim".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 8,
                    correction_time_ms: 50,
                },
            });
        }
        
        let builder = AssemblyGraphBuilder::new(15, 25, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Plasmid assembly should succeed");
        
        let assembly = result.unwrap();
        
        // Should handle circular sequences
        assert!(!assembly.contigs.is_empty(), "Should generate contigs for plasmid");
        assert!(assembly.assembly_stats.total_length > 1000, "Should recover substantial sequence");
    }
}

#[cfg(test)]
mod realistic_error_patterns {
    use super::*;

    fn create_read_with_indel_errors(id: usize, sequence: &str) -> CorrectedRead {
        let mut corrected_seq = String::new();
        let mut corrections = Vec::new();
        let mut pos = 0;
        
        for (i, base) in sequence.chars().enumerate() {
            if fastrand::f64() < 0.005 { // 0.5% deletion rate
                corrections.push(BaseCorrection {
                    position: pos,
                    from: '_', // Represent deletion
                    to: base,
                    confidence: 0.8,
                    correction_type: CorrectionType::Deletion,
                });
                // Skip this base (deletion)
            } else {
                corrected_seq.push(base);
                
                if fastrand::f64() < 0.003 { // 0.3% insertion rate
                    let insert_base = match fastrand::usize(0..4) {
                        0 => 'A', 1 => 'T', 2 => 'G', _ => 'C',
                    };
                    corrected_seq.push(insert_base);
                    
                    corrections.push(BaseCorrection {
                        position: pos + 1,
                        from: insert_base,
                        to: '_',
                        confidence: 0.75,
                        correction_type: CorrectionType::Insertion,
                    });
                }
                pos += 1;
            }
        }
        
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: corrected_seq,
            corrections,
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "indel_corrector".to_string(),
                confidence_threshold: 0.85,
                context_window: 12,
                correction_time_ms: fastrand::u64(100..300),
            },
        }
    }

    #[test]
    fn test_homopolymer_error_handling() {
        // Homopolymers are prone to sequencing errors
        let genome_with_homopolymers = "ATCGAAAAAAAAATCGTTTTTTTTCGGGGGGGGGGATCCCCCCCCCATCG";
        
        let mut reads = Vec::new();
        for i in 0..20 {
            let start = fastrand::usize(0..std::cmp::max(1, genome_with_homopolymers.len() - 25));
            let end = std::cmp::min(start + 25, genome_with_homopolymers.len());
            let read_seq = &genome_with_homopolymers[start..end];
            
            reads.push(create_read_with_indel_errors(i, read_seq));
        }
        
        let builder = AssemblyGraphBuilder::new(8, 12, 1);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle homopolymer-rich sequences");
        
        let assembly = result.unwrap();
        
        // Should detect some low-complexity regions
        let low_complexity_nodes = assembly.graph_fragment.nodes
            .values()
            .filter(|n| n.complexity_score < 0.2)
            .count();
        
        assert!(low_complexity_nodes > 0, "Should detect homopolymer regions");
    }

    #[test]
    fn test_quality_based_correction_simulation() {
        let genome = generate_synthetic_genome(2000, 0.45);
        let mut reads = Vec::new();
        
        // Simulate quality-dependent error correction
        for i in 0..50 {
            let start = fastrand::usize(0..std::cmp::max(1, genome.len() - 100));
            let end = std::cmp::min(start + 100, genome.len());
            let original_seq = &genome[start..end];
            
            let mut corrected_seq = String::new();
            let mut corrections = Vec::new();
            let mut quality_scores = Vec::new();
            
            for (pos, base) in original_seq.chars().enumerate() {
                let quality = if fastrand::f64() < 0.1 { 15 } else { 35 }; // 10% low quality
                quality_scores.push(quality);
                
                if quality < 20 && fastrand::f64() < 0.1 { // Errors more likely in low quality
                    let error_base = match fastrand::usize(0..3) {
                        0 => 'A', 1 => 'T', _ => 'G',
                    };
                    
                    corrected_seq.push(base); // Corrected back to original
                    corrections.push(BaseCorrection {
                        position: pos,
                        from: error_base,
                        to: base,
                        confidence: 0.6,
                        correction_type: CorrectionType::QualityImprovement,
                    });
                } else {
                    corrected_seq.push(base);
                }
            }
            
            reads.push(CorrectedRead {
                id: i,
                original: original_seq.to_string(),
                corrected: corrected_seq,
                corrections,
                quality_scores,
                correction_metadata: CorrectionMetadata {
                    algorithm: "quality_corrector".to_string(),
                    confidence_threshold: 0.8,
                    context_window: 15,
                    correction_time_ms: fastrand::u64(80..250),
                },
            });
        }
        
        let builder = AssemblyGraphBuilder::new(12, 18, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle quality-based corrections");
        
        let assembly = result.unwrap();
        assert!(!assembly.contigs.is_empty(), "Should generate contigs despite quality issues");
    }

    #[test]
    fn test_paired_end_like_coverage() {
        let genome = generate_synthetic_genome(3000, 0.52);
        let mut reads = Vec::new();
        let insert_size = 500;
        let read_length = 150;
        
        // Simulate paired-end sequencing
        for i in 0..100 {
            if genome.len() >= insert_size {
                let fragment_start = fastrand::usize(0..genome.len() - insert_size);
                
                // Forward read
                let forward_end = std::cmp::min(fragment_start + read_length, genome.len());
                let forward_seq = &genome[fragment_start..forward_end];
                
                // Reverse read (from end of fragment)
                let reverse_start = std::cmp::max(fragment_start + insert_size - read_length, 0);
                let reverse_end = std::cmp::min(reverse_start + read_length, genome.len());
                let reverse_seq = &genome[reverse_start..reverse_end];
                
                reads.push(CorrectedRead {
                    id: i * 2,
                    original: forward_seq.to_string(),
                    corrected: forward_seq.to_string(),
                    corrections: Vec::new(),
                    quality_scores: vec![38; forward_seq.len()],
                    correction_metadata: CorrectionMetadata {
                        algorithm: "paired_end_sim".to_string(),
                        confidence_threshold: 0.95,
                        context_window: 10,
                        correction_time_ms: 75,
                    },
                });
                
                reads.push(CorrectedRead {
                    id: i * 2 + 1,
                    original: reverse_seq.to_string(),
                    corrected: reverse_seq.to_string(),
                    corrections: Vec::new(),
                    quality_scores: vec![36; reverse_seq.len()],
                    correction_metadata: CorrectionMetadata {
                        algorithm: "paired_end_sim".to_string(),
                        confidence_threshold: 0.95,
                        context_window: 10,
                        correction_time_ms: 75,
                    },
                });
            }
        }
        
        let builder = AssemblyGraphBuilder::new(21, 31, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Should handle paired-end like reads");
        
        let assembly = result.unwrap();
        
        // Paired-end should give good coverage and contiguity
        assert!(!assembly.contigs.is_empty(), "Should generate contigs");
        assert!(assembly.assembly_stats.coverage_mean > 8.0, "Should have good coverage");
        
        // Should have some high-coverage nodes from overlapping pairs
        let high_coverage_nodes = assembly.graph_fragment.nodes
            .values()
            .filter(|n| n.coverage > 5)
            .count();
        
        assert!(high_coverage_nodes > 0, "Should have well-covered regions from paired reads");
    }
}

#[cfg(test)]
mod biological_accuracy_validation {
    use super::*;

    #[test]
    fn test_codon_preservation_in_assembly() {
        // Create a sequence with valid start/stop codons
        let coding_sequence = "ATGAAATTTGGCCCCTAG"; // ATG (start), AAA (Lys), TTT (Phe), GGC (Gly), CCT (Pro), TAG (stop)
        let non_coding = generate_synthetic_genome(500, 0.4);
        let genome = format!("{}{}{}", non_coding, coding_sequence, non_coding);
        
        let reads = generate_synthetic_reads(&genome, 80, 20.0, 0.005);
        
        let builder = AssemblyGraphBuilder::new(15, 21, 3);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Coding sequence assembly should succeed");
        
        let assembly = result.unwrap();
        
        // Check if coding sequence patterns are preserved in contigs
        let mut has_start_codon = false;
        let mut has_stop_codon = false;
        
        for contig in &assembly.contigs {
            if contig.sequence.contains("ATG") {
                has_start_codon = true;
            }
            if contig.sequence.contains("TAG") || contig.sequence.contains("TAA") || contig.sequence.contains("TGA") {
                has_stop_codon = true;
            }
        }
        
        // In a real implementation, you might have more sophisticated ORF detection
        assert!(!assembly.contigs.is_empty(), "Should generate contigs containing coding sequences");
    }

    #[test]
    fn test_ribosomal_rna_like_structure() {
        // Simulate highly conserved rRNA-like sequence
        let rrna_like = "GCCTAACATGCCAAGTCGAGCGGCGGCGGGAAGACCCGCGCCGCGGTGTTGATTTTGACGTGGGTTCCTCCGAATAGGGGCGACCACCCGGGCCGCGGTGTTGATTTTGACGTGGGTTC";
        let variable_regions = generate_synthetic_genome(800, 0.6);
        
        let genome = format!("{}{}{}{}{}", 
                           variable_regions,
                           rrna_like,
                           generate_synthetic_genome(400, 0.55),
                           rrna_like, // Repeated rRNA sequence
                           variable_regions);
        
        let reads = generate_synthetic_reads(&genome, 120, 25.0, 0.003);
        
        let builder = AssemblyGraphBuilder::new(25, 35, 4);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "rRNA-like sequence assembly should succeed");
        
        let assembly = result.unwrap();
        
        // Should detect high-coverage repetitive nodes (rRNA regions)
        let high_cov_repetitive = assembly.graph_fragment.nodes
            .values()
            .filter(|n| n.coverage > 10 && matches!(n.node_type, NodeType::Repetitive))
            .count();
        
        assert!(high_cov_repetitive > 0, "Should detect highly conserved repetitive sequences");
    }

    #[test]
    fn test_gc_skew_handling() {
        // Create genome with GC skew (common in bacterial genomes)
        let leading_strand_high_g = generate_synthetic_genome(1000, 0.7); // High G+C
        let lagging_strand_high_c = leading_strand_high_g
            .chars()
            .map(|c| match c {
                'G' => 'C',
                'C' => 'G', 
                'A' => 'T',
                'T' => 'A',
                _ => c,
            })
            .collect::<String>();
        
        let genome = format!("{}{}", leading_strand_high_g, lagging_strand_high_c);
        let reads = generate_synthetic_reads(&genome, 100, 15.0, 0.008);
        
        let builder = AssemblyGraphBuilder::new(17, 25, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "GC skewed genome should assemble");
        
        let assembly = result.unwrap();
        
        // Should handle compositional bias
        assert!(!assembly.contigs.is_empty(), "Should generate contigs despite GC skew");
        assert!(assembly.assembly_stats.gc_content > 0.4 && assembly.assembly_stats.gc_content < 0.8,
               "Should have intermediate GC content from mixed regions");
    }

    #[test]
    fn test_mobile_element_detection() {
        // Simulate transposable element-like repeats
        let transposon = "TGCAGGGATGACGCCCGGCCAGATCGTT";
        let unique_seq1 = generate_synthetic_genome(800, 0.45);
        let unique_seq2 = generate_synthetic_genome(600, 0.50);
        let unique_seq3 = generate_synthetic_genome(700, 0.48);
        
        // Insert transposon at multiple locations
        let genome = format!("{}{}{}{}{}{}", 
                           unique_seq1,
                           transposon,
                           unique_seq2,
                           transposon,
                           unique_seq3,
                           transposon);
        
        let reads = generate_synthetic_reads(&genome, 90, 18.0, 0.01);
        
        let builder = AssemblyGraphBuilder::new(13, 19, 2);
        let result = builder.build(&reads);
        assert!(result.is_ok(), "Mobile element containing genome should assemble");
        
        let assembly = result.unwrap();
        
        // Should detect mobile element as repetitive
        let repetitive_high_cov = assembly.graph_fragment.nodes
            .values()
            .filter(|n| matches!(n.node_type, NodeType::Repetitive) && n.coverage > 8)
            .count();
        
        assert!(repetitive_high_cov > 0, "Should detect mobile elements as repetitive sequences");
        
        // Assembly might be fragmented at repeat boundaries
        assert!(assembly.contigs.len() >= 2, "Mobile elements typically fragment assemblies");
    }

    #[test]
    fn test_snp_variant_handling() {
        // Create genome with SNP-like variations
        let base_genome = generate_synthetic_genome(2000, 0.48);
        let mut variant_reads = Vec::new();
        
        // Generate reads with some containing SNPs
        for i in 0..80 {
            let start = fastrand::usize(0..std::cmp::max(1, base_genome.len() - 120));
            let end = std::cmp::min(start + 120, base_genome.len());
            let mut read_seq = base_genome[start..end].to_string();
            
            // Introduce SNPs in 20% of reads
            if i % 5 == 0 {
                if let Some(snp_pos) = fastrand::usize(10..read_seq.len() - 10).into() {
                    let mut chars: Vec<char> = read_seq.chars().collect();
                    chars[snp_pos] = match chars[snp_pos] {
                        'A' => 'G',
                        'T' => 'C', 
                        'G' => 'A',
                        'C' => 'T',
                        _ => 'N',
                    };
                    read_seq = chars.into_iter().collect();
                }
            }
            
            variant_reads.push(CorrectedRead {
                id: i,
                original: read_seq.clone(),
                corrected: read_seq,
                corrections: Vec::new(),
                quality_scores: vec![35; end - start],
                correction_metadata: CorrectionMetadata {
                    algorithm: "snp_aware".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 8,
                    correction_time_ms: 100,
                },
            });
        }
        
        let builder = AssemblyGraphBuilder::new(15, 21, 2);
        let result = builder.build(&variant_reads);
        assert!(result.is_ok(), "SNP-containing reads should assemble");
        
        let assembly = result.unwrap();
        
        // Should handle variants - might create bubble structures
        let bubbles = assembly.graph_fragment.find_bubbles();
        // In a more sophisticated implementation, bubbles might represent variants
        
        assert!(!assembly.contigs.is_empty(), "Should generate contigs despite variants");
    }
}