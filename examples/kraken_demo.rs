use anyhow::Result;
use meta_forge::utils::kraken_reporter::{
    ClassificationStatus, ConfidenceMetrics, KrakenClassification, KrakenReporter,
    TaxonomicLineage, TaxonomicRank,
};
use std::path::Path;

/// Demonstration of the new Kraken-style reporting functionality
fn main() -> Result<()> {
    println!("🧬 Kraken-style Metagenomic Report Demo");
    println!("=====================================");

    // Create a new Kraken reporter
    let mut reporter = KrakenReporter::new("demo_sample".to_string());

    // Create example classifications with hierarchical taxonomy
    let classifications = create_demo_classifications();

    // Add classifications to reporter
    for classification in classifications {
        reporter.add_classification(classification);
    }

    // Create output directory
    std::fs::create_dir_all("output/demo")?;

    // Generate Kraken-style output files
    println!("\n📄 Generating Kraken-style reports...");

    // Standard Kraken output format
    let kraken_path = Path::new("output/demo/demo_sample.kraken");
    reporter.write_kraken_output(&kraken_path)?;
    println!("✅ Standard Kraken output: {}", kraken_path.display());

    // Kraken-style abundance report
    let report_path = Path::new("output/demo/demo_sample_kraken_report.txt");
    reporter.write_kraken_report(&report_path)?;
    println!("✅ Kraken abundance report: {}", report_path.display());

    // Enhanced JSON report with detailed information
    let json_path = Path::new("output/demo/demo_sample_enhanced.json");
    reporter.write_enhanced_json_report(&json_path)?;
    println!("✅ Enhanced JSON report: {}", json_path.display());

    // Display sample of the reports
    show_sample_outputs(&kraken_path, &report_path)?;

    println!("\n🎉 Demo completed successfully!");
    println!("📁 Output files are available in: output/demo/");

    Ok(())
}

fn create_demo_classifications() -> Vec<KrakenClassification> {
    vec![
        // E. coli classification
        KrakenClassification {
            classification_status: ClassificationStatus::Classified,
            sequence_id: "contig_001".to_string(),
            taxonomy_id: 511145,
            sequence_length: 1547,
            kmer_lca_mappings: vec![511145, 511145, 562, 511145, 561, 511145],
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: 0.92,
                classified_kmers: 147,
                total_kmers: 160,
                confidence_score: 0.92,
                classification_method: "exact_kmer_match".to_string(),
            },
            taxonomic_lineage: TaxonomicLineage {
                domain: Some(TaxonomicRank {
                    taxonomy_id: 2,
                    name: "Bacteria".to_string(),
                    rank: "superkingdom".to_string(),
                }),
                kingdom: None,
                phylum: Some(TaxonomicRank {
                    taxonomy_id: 1224,
                    name: "Proteobacteria".to_string(),
                    rank: "phylum".to_string(),
                }),
                class: Some(TaxonomicRank {
                    taxonomy_id: 1236,
                    name: "Gammaproteobacteria".to_string(),
                    rank: "class".to_string(),
                }),
                order: Some(TaxonomicRank {
                    taxonomy_id: 91347,
                    name: "Enterobacterales".to_string(),
                    rank: "order".to_string(),
                }),
                family: Some(TaxonomicRank {
                    taxonomy_id: 543,
                    name: "Enterobacteriaceae".to_string(),
                    rank: "family".to_string(),
                }),
                genus: Some(TaxonomicRank {
                    taxonomy_id: 561,
                    name: "Escherichia".to_string(),
                    rank: "genus".to_string(),
                }),
                species: Some(TaxonomicRank {
                    taxonomy_id: 511145,
                    name: "Escherichia coli str. K-12 substr. MG1655".to_string(),
                    rank: "strain".to_string(),
                }),
            },
        },
        // Salmonella classification
        KrakenClassification {
            classification_status: ClassificationStatus::Classified,
            sequence_id: "contig_002".to_string(),
            taxonomy_id: 28901,
            sequence_length: 892,
            kmer_lca_mappings: vec![28901, 590, 28901, 28901],
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: 0.78,
                classified_kmers: 69,
                total_kmers: 89,
                confidence_score: 0.78,
                classification_method: "lca_mapping".to_string(),
            },
            taxonomic_lineage: TaxonomicLineage {
                domain: Some(TaxonomicRank {
                    taxonomy_id: 2,
                    name: "Bacteria".to_string(),
                    rank: "superkingdom".to_string(),
                }),
                kingdom: None,
                phylum: Some(TaxonomicRank {
                    taxonomy_id: 1224,
                    name: "Proteobacteria".to_string(),
                    rank: "phylum".to_string(),
                }),
                class: Some(TaxonomicRank {
                    taxonomy_id: 1236,
                    name: "Gammaproteobacteria".to_string(),
                    rank: "class".to_string(),
                }),
                order: Some(TaxonomicRank {
                    taxonomy_id: 91347,
                    name: "Enterobacterales".to_string(),
                    rank: "order".to_string(),
                }),
                family: Some(TaxonomicRank {
                    taxonomy_id: 543,
                    name: "Enterobacteriaceae".to_string(),
                    rank: "family".to_string(),
                }),
                genus: Some(TaxonomicRank {
                    taxonomy_id: 590,
                    name: "Salmonella".to_string(),
                    rank: "genus".to_string(),
                }),
                species: Some(TaxonomicRank {
                    taxonomy_id: 28901,
                    name: "Salmonella enterica".to_string(),
                    rank: "species".to_string(),
                }),
            },
        },
        // Unclassified sequence
        KrakenClassification {
            classification_status: ClassificationStatus::Unclassified,
            sequence_id: "contig_003".to_string(),
            taxonomy_id: 0,
            sequence_length: 423,
            kmer_lca_mappings: vec![],
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: 0.0,
                classified_kmers: 0,
                total_kmers: 42,
                confidence_score: 0.0,
                classification_method: "none".to_string(),
            },
            taxonomic_lineage: TaxonomicLineage {
                domain: None,
                kingdom: None,
                phylum: None,
                class: None,
                order: None,
                family: None,
                genus: None,
                species: None,
            },
        },
        // Bacillus classification
        KrakenClassification {
            classification_status: ClassificationStatus::Classified,
            sequence_id: "contig_004".to_string(),
            taxonomy_id: 1396,
            sequence_length: 2134,
            kmer_lca_mappings: vec![1396, 1386, 1396, 1396, 186817, 1396],
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: 0.85,
                classified_kmers: 181,
                total_kmers: 213,
                confidence_score: 0.85,
                classification_method: "exact_kmer_match".to_string(),
            },
            taxonomic_lineage: TaxonomicLineage {
                domain: Some(TaxonomicRank {
                    taxonomy_id: 2,
                    name: "Bacteria".to_string(),
                    rank: "superkingdom".to_string(),
                }),
                kingdom: None,
                phylum: Some(TaxonomicRank {
                    taxonomy_id: 1239,
                    name: "Firmicutes".to_string(),
                    rank: "phylum".to_string(),
                }),
                class: Some(TaxonomicRank {
                    taxonomy_id: 91061,
                    name: "Bacilli".to_string(),
                    rank: "class".to_string(),
                }),
                order: Some(TaxonomicRank {
                    taxonomy_id: 1385,
                    name: "Bacillales".to_string(),
                    rank: "order".to_string(),
                }),
                family: Some(TaxonomicRank {
                    taxonomy_id: 186817,
                    name: "Bacillaceae".to_string(),
                    rank: "family".to_string(),
                }),
                genus: Some(TaxonomicRank {
                    taxonomy_id: 1386,
                    name: "Bacillus".to_string(),
                    rank: "genus".to_string(),
                }),
                species: Some(TaxonomicRank {
                    taxonomy_id: 1396,
                    name: "Bacillus subtilis".to_string(),
                    rank: "species".to_string(),
                }),
            },
        },
    ]
}

fn show_sample_outputs(kraken_path: &Path, report_path: &Path) -> Result<()> {
    println!("\n📄 Sample Kraken Output:");
    println!("========================");
    let kraken_content = std::fs::read_to_string(kraken_path)?;
    println!("{}", kraken_content);

    println!("\n📊 Sample Kraken Report (first 15 lines):");
    println!("==========================================");
    let report_content = std::fs::read_to_string(report_path)?;
    let lines: Vec<&str> = report_content.lines().take(15).collect();
    for line in lines {
        println!("{}", line);
    }
    if report_content.lines().count() > 15 {
        println!("... (truncated)");
    }

    Ok(())
}
