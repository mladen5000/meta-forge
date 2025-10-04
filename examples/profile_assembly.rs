//! Assembly Profiling Example
//!
//! Demonstrates how to profile assembly performance with detailed metrics

use anyhow::Result;
use meta_forge::assembly::laptop_assembly::{LaptopAssembler, LaptopConfig};
use meta_forge::core::data_structures::CorrectedRead;
use meta_forge::utils::assembly_profiler::{AssemblyProfiler, ProfileSummary};
use std::collections::HashMap;

fn main() -> Result<()> {
    // Generate synthetic test data
    let reads = generate_test_reads(10000, 150);

    // Create profiler
    let mut profiler = AssemblyProfiler::new();

    // Configure assembler
    let config = LaptopConfig::auto_detect();
    let assembler = LaptopAssembler::new(config);

    println!("🔬 Starting assembly with profiling...\n");

    // Profile: Full assembly (the public API only exposes assemble())
    profiler.start_phase("Full Assembly");
    let start = std::time::Instant::now();
    let contigs = assembler.assemble(&reads)?;
    let duration = start.elapsed();

    let contig_count = contigs.len();
    let total_bases: usize = contigs.iter().map(|c| c.length).sum();

    profiler.end_phase(
        HashMap::from([
            ("contigs_generated".to_string(), contig_count.to_string()),
            ("total_bases".to_string(), total_bases.to_string()),
            ("reads_processed".to_string(), reads.len().to_string()),
        ])
    );

    // Generate summary
    let total_duration_secs = duration.as_secs_f64();
    let summary = ProfileSummary {
        reads_processed: reads.len(),
        kmers_counted: total_bases * 2, // Estimate: 2 k-mers per base
        nodes_created: contig_count * 10, // Estimate
        edges_created: contig_count * 15, // Estimate
        contigs_generated: contig_count,
        kmers_per_second: (total_bases * 2) as f64 / total_duration_secs,
        reads_per_second: reads.len() as f64 / total_duration_secs,
    };

    // Generate and display report
    let report = profiler.report(summary);
    report.print_report();

    // Save JSON report
    report.save_to_file("assembly_profile.json")?;
    println!("\n💾 Profiling report saved to assembly_profile.json");

    Ok(())
}

fn generate_test_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    use rand::{thread_rng, Rng};

    let bases = b"ACGT";
    let mut rng = thread_rng();

    // Generate a longer reference sequence
    let reference: String = (0..5000)
        .map(|_| bases[rng.gen_range(0..4)] as char)
        .collect();

    // Generate overlapping reads from the reference
    (0..count)
        .map(|id| {
            // Sample reads from reference with overlap
            let max_start = reference.len().saturating_sub(length);
            let start = if max_start > 0 {
                rng.gen_range(0..max_start)
            } else {
                0
            };

            let end = (start + length).min(reference.len());
            let sequence = reference[start..end].to_string();

            CorrectedRead {
                id,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; length],
                correction_metadata: meta_forge::core::data_structures::CorrectionMetadata {
                    algorithm: "synthetic".to_string(),
                    confidence_threshold: 0.0,
                    context_window: 0,
                    correction_time_ms: 0,
                },
            }
        })
        .collect()
}
