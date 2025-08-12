use anyhow::Result;
use meta_forge::utils::intermediate_output::{
    IntermediateOutputManager, OutputConfig, PipelineSection,
};
use serde_json::json;
use std::collections::HashMap;

/// Demonstration of the intermediate output system
fn main() -> Result<()> {
    println!("🧬 MetaForge Intermediate Output System Demo");
    println!("===========================================\n");

    // Create output configuration
    let config = OutputConfig {
        enable_json: true,
        enable_binary: false, // Simplified for demo
        enable_fasta: true,
        enable_tsv: true,
        compress_files: false, // Disabled for easier inspection
        max_file_size_mb: 100,
    };

    // Initialize output manager - creates run_DDMMYY_HHMMSS directory
    let output_manager =
        IntermediateOutputManager::new(std::path::PathBuf::from("output"), config)?;

    println!(
        "📁 Created run directory: {}",
        output_manager.run_dir.display()
    );
    println!("🆔 Run ID: {}\n", output_manager.run_id);

    // Simulate pipeline data flow with intermediate saves

    // 1. Preprocessing Phase
    println!("📋 Phase 1: Preprocessing");
    let preprocessing_data = json!({
        "reads_processed": 15000,
        "corrections_made": 245,
        "quality_score": 0.94,
        "total_length": 75000000,
        "processing_time_ms": 12500
    });

    output_manager.save_intermediate(
        PipelineSection::Preprocessing,
        "corrected_reads",
        &preprocessing_data,
        json!({
            "algorithm": "statistical_consensus",
            "parameters": {"k": 15, "threshold": 0.8}
        }),
    )?;
    println!("   ✅ Saved preprocessing intermediate files");

    // 2. Assembly Phase
    println!("🧬 Phase 2: Assembly");
    let assembly_data = json!({
        "num_contigs": 1247,
        "total_length": 68500000,
        "n50": 45000,
        "max_contig_length": 125000,
        "coverage_mean": 25.5,
        "k_range": [21, 31]
    });

    output_manager.save_intermediate(
        PipelineSection::Assembly,
        "assembly_results",
        &assembly_data,
        json!({
            "assembler": "adaptive_k_mer",
            "graph_nodes": 45000,
            "graph_edges": 78000
        }),
    )?;

    // Save contigs as FASTA
    let contigs = vec![
        (
            "contig_1".to_string(),
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
        ),
        (
            "contig_2".to_string(),
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA".to_string(),
        ),
        (
            "contig_3".to_string(),
            "TTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC".to_string(),
        ),
    ];

    output_manager.save_sequences(
        PipelineSection::Assembly,
        "contigs",
        &contigs,
        json!({"num_contigs": contigs.len()}),
    )?;
    println!("   ✅ Saved assembly intermediate files (JSON + FASTA)");

    // 3. Classification Phase
    println!("🏷️  Phase 3: Classification");
    let classification_data = json!({
        "classified_contigs": 1180,
        "unclassified_contigs": 67,
        "unique_taxa": 15,
        "mean_confidence": 0.87,
        "top_species": ["E. coli", "B. subtilis", "S. aureus"]
    });

    output_manager.save_intermediate(
        PipelineSection::Classification,
        "taxonomic_results",
        &classification_data,
        json!({
            "classifier": "ML_ensemble",
            "database": "NCBI_RefSeq",
            "features_used": ["composition", "k_mers", "codon_usage"]
        }),
    )?;

    // Save classification details as TSV
    let headers = vec![
        "contig_id".to_string(),
        "taxonomy".to_string(),
        "confidence".to_string(),
        "lineage".to_string(),
    ];
    let rows = vec![
        vec![
            "contig_1".to_string(),
            "Escherichia coli".to_string(),
            "0.95".to_string(),
            "Bacteria;Proteobacteria;Gammaproteobacteria".to_string(),
        ],
        vec![
            "contig_2".to_string(),
            "Bacillus subtilis".to_string(),
            "0.89".to_string(),
            "Bacteria;Firmicutes;Bacilli".to_string(),
        ],
        vec![
            "contig_3".to_string(),
            "Staphylococcus aureus".to_string(),
            "0.92".to_string(),
            "Bacteria;Firmicutes;Bacilli".to_string(),
        ],
    ];

    output_manager.save_tsv(
        PipelineSection::Classification,
        "detailed_classifications",
        &headers,
        &rows,
        json!({"format": "NCBI_taxonomy"}),
    )?;
    println!("   ✅ Saved classification intermediate files (JSON + TSV)");

    // 4. Abundance Phase
    println!("📊 Phase 4: Abundance");
    let mut abundant_kmers = HashMap::new();
    abundant_kmers.insert(1u64, 150.5);
    abundant_kmers.insert(2u64, 95.2);
    abundant_kmers.insert(3u64, 78.9);

    let abundance_data = json!({
        "unique_kmers": 45000,
        "total_kmers": 750000,
        "abundant_kmers": abundant_kmers,
        "diversity_index": 3.2,
        "evenness": 0.75
    });

    output_manager.save_intermediate(
        PipelineSection::Abundance,
        "abundance_profile",
        &abundance_data,
        json!({
            "method": "HyperLogLog_L0_sampling",
            "k_size": 21,
            "sampling_rate": 0.01
        }),
    )?;
    println!("   ✅ Saved abundance intermediate files");

    // 5. Final Report Phase
    println!("📝 Phase 5: Final Report");
    let report_data = json!({
        "sample_name": "demo_sample",
        "timestamp": chrono::Utc::now().to_rfc3339(),
        "summary": {
            "total_contigs": 1247,
            "total_length": 68500000,
            "n50": 45000,
            "unique_species": 15,
            "diversity_index": 3.2
        },
        "quality_metrics": {
            "assembly_completeness": 0.92,
            "classification_confidence": 0.87,
            "coverage_uniformity": 0.78
        }
    });

    output_manager.save_intermediate(
        PipelineSection::Report,
        "final_report",
        &report_data,
        json!({
            "format_version": "1.0",
            "generator": "MetaForge_v0.4.0"
        }),
    )?;
    println!("   ✅ Saved final report intermediate files");

    // Generate and display run summary
    println!("\n📋 Run Summary");
    println!("==============");
    let summary = output_manager.generate_run_summary()?;
    println!("Run ID: {}", summary.run_id);
    println!("Run Directory: {}", summary.run_dir.display());
    println!(
        "Timestamp: {}",
        summary.timestamp.format("%Y-%m-%d %H:%M:%S UTC")
    );

    println!("\nFiles by Section:");
    for (section, files) in &summary.sections {
        println!("  📁 {} ({} files)", section.as_str(), files.len());
        for file in files {
            println!("    📄 {} ({} bytes)", file.name, file.size);
        }
    }

    // Demonstrate restart capability
    println!("\n🔄 Restart Capability Demo");
    println!("=========================");

    println!("Available checkpoints:");
    let sections = [
        PipelineSection::Assembly,
        PipelineSection::Features,
        PipelineSection::Classification,
        PipelineSection::Abundance,
        PipelineSection::Report,
    ];

    for section in &sections {
        let has_data = output_manager.has_intermediate(
            section.clone(),
            match section {
                PipelineSection::Assembly => "assembly_results",
                PipelineSection::Classification => "taxonomic_results",
                PipelineSection::Abundance => "abundance_profile",
                PipelineSection::Report => "final_report",
                _ => "data",
            },
        );
        let status = if has_data {
            "✅ Available"
        } else {
            "❌ Missing"
        };
        println!("  🔄 {} - {}", section.as_str(), status);
    }

    println!("\nTo restart from a checkpoint, use:");
    println!(
        "  meta-pipeline resume --run-id {} --section classification --sample-name demo_sample",
        summary.run_id
    );

    println!("\n🎉 Intermediate Output System Demo Complete!");
    println!("📁 All files saved to: {}", summary.run_dir.display());

    Ok(())
}
