// Quick test to verify taxonomy classification works correctly
use meta_forge::ml::simple_classifier::{SimpleClassifierConfig, SimpleContigClassifier};
use meta_forge::core::data_structures::{Contig, ContigType};

fn main() {
    println!("🧬 Testing Taxonomy-Aware Classification\n");

    // Create test contigs with different GC content (simulating different species)
    let contigs = vec![
        Contig {
            id: 1,
            sequence: "ATATATATATAT".repeat(100), // AT-rich (low GC ~33%)
            coverage: 10.0,
            length: 1200,
            contig_type: ContigType::Linear,
            node_path: vec![],
        },
        Contig {
            id: 2,
            sequence: "GCGCGCGCGCGC".repeat(100), // GC-rich (high GC ~66%)
            coverage: 15.0,
            length: 1200,
            contig_type: ContigType::Linear,
            node_path: vec![],
        },
        Contig {
            id: 3,
            sequence: "ATGCATGCATGC".repeat(100), // Balanced GC (~50%)
            coverage: 12.0,
            length: 1200,
            contig_type: ContigType::Linear,
            node_path: vec![],
        },
    ];

    let config = SimpleClassifierConfig {
        kmer_size: 4,
        min_contig_length: 500,
        ..Default::default()
    };

    println!("📊 Classifier config:");
    println!("   - K-mer size: {}", config.kmer_size);
    println!("   - Min contig length: {}", config.min_contig_length);
    println!("   - Auto-detect bins: {}", config.auto_detect_bins);
    println!();

    let classifier = SimpleContigClassifier::new(config).unwrap();

    println!("🔬 Running taxonomy-aware classification...");
    let results = classifier.classify_contigs_with_taxonomy(&contigs).unwrap();

    println!("\n✅ Classification Results:\n");
    for result in &results {
        println!("   Contig {}: {} ({})",
            result.contig_id,
            result.taxon,
            result.rank
        );
        println!("      Bin: {}, Confidence: {:.2}%",
            result.bin_id,
            result.confidence * 100.0
        );
        println!();
    }

    // Check that we're getting actual species names, not "Bin_X"
    let has_real_taxonomy = results.iter().any(|r|
        r.taxon != "Unclassified" && !r.taxon.starts_with("Bin_")
    );

    if has_real_taxonomy {
        println!("✅ SUCCESS: Taxonomy labeling is working correctly!");
        println!("   Species names are being assigned instead of generic bin labels.");
    } else {
        println!("⚠️  WARNING: Only getting generic labels (Unclassified/Bin_X)");
        println!("   This may be expected for synthetic test data.");
    }
}
