//! Debug test to verify coverage filtering behavior
//! This will test the AdvancedAssemblyGraphBuilder directly with min_coverage=1

use meta_forge::assembly::graph_construction::AdvancedAssemblyGraphBuilder;
use meta_forge::core::data_structures::*;

#[test]
fn test_coverage_filtering_with_min_coverage_1() {
    println!("🧪 Testing AdvancedAssemblyGraphBuilder with min_coverage=1");

    // Create test reads
    let test_reads = vec![
        CorrectedRead {
            id: 0,
            original: "ATCGATCGATCGATCGATCG".to_string(),
            corrected: "ATCGATCGATCGATCGATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 20],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 1,
            original: "TCGATCGATCGATCGATCGA".to_string(),
            corrected: "TCGATCGATCGATCGATCGA".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 20],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
    ];

    println!("📊 Test data: {} reads", test_reads.len());

    // Test with min_coverage=1 (should keep all nodes)
    println!("\n🔧 Testing with min_coverage=1...");
    let builder1 = AdvancedAssemblyGraphBuilder::new(15, 31, 1, 8).unwrap(); // min_coverage=1

    match builder1.build_graph(&test_reads) {
        Ok(assembly_graph) => {
            println!("✅ Graph construction succeeded with min_coverage=1!");
            println!("📊 Results with min_coverage=1:");
            println!(
                "  Graph fragment nodes: {}",
                assembly_graph.graph_fragment.nodes.len()
            );
            println!(
                "  Graph fragment edges: {}",
                assembly_graph.graph_fragment.edges.len()
            );
            println!("  Contigs: {}", assembly_graph.contigs.len());

            if assembly_graph.graph_fragment.nodes.len() > 0 {
                println!("🎉 SUCCESS: min_coverage=1 preserves nodes!");

                // Print node details
                for (hash, node) in &assembly_graph.graph_fragment.nodes {
                    println!(
                        "  Node {}: k-mer='{}', coverage={}",
                        hash, node.kmer.sequence, node.coverage
                    );
                }
            } else {
                println!("🚨 ISSUE: Even with min_coverage=1, no nodes created");
            }
        }
        Err(e) => {
            println!("❌ Graph construction failed with min_coverage=1: {}", e);
        }
    }

    // Test with min_coverage=2 (should filter out low-coverage nodes)
    println!("\n🔧 Testing with min_coverage=2...");
    let builder2 = AdvancedAssemblyGraphBuilder::new(15, 31, 2, 8).unwrap(); // min_coverage=2

    match builder2.build_graph(&test_reads) {
        Ok(assembly_graph) => {
            println!("✅ Graph construction succeeded with min_coverage=2!");
            println!("📊 Results with min_coverage=2:");
            println!(
                "  Graph fragment nodes: {}",
                assembly_graph.graph_fragment.nodes.len()
            );
            println!(
                "  Graph fragment edges: {}",
                assembly_graph.graph_fragment.edges.len()
            );
            println!("  Contigs: {}", assembly_graph.contigs.len());

            if assembly_graph.graph_fragment.nodes.len() == 0 {
                println!("✅ EXPECTED: min_coverage=2 filters out single-coverage nodes");
            } else {
                println!("🚨 UNEXPECTED: min_coverage=2 still has nodes:");

                // Print node details
                for (hash, node) in &assembly_graph.graph_fragment.nodes {
                    println!(
                        "  Node {}: k-mer='{}', coverage={}",
                        hash, node.kmer.sequence, node.coverage
                    );
                }
            }
        }
        Err(e) => {
            println!("❌ Graph construction failed with min_coverage=2: {}", e);
        }
    }
}

#[test]
fn test_higher_coverage_dataset() {
    println!("🧪 Testing with higher coverage dataset");

    // Create reads with overlaps to generate k-mers with coverage > 1
    let test_reads = vec![
        // These reads have significant overlaps, so some k-mers should appear multiple times
        CorrectedRead {
            id: 0,
            original: "ATCGATCGATCGATCGATCGATCGATCGATCG".to_string(), // 32bp
            corrected: "ATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 32],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 1,
            original: "TCGATCGATCGATCGATCGATCGATCGATCGA".to_string(), // Overlaps by 31bp
            corrected: "TCGATCGATCGATCGATCGATCGATCGATCGA".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 32],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 2,
            original: "CGATCGATCGATCGATCGATCGATCGATCGAT".to_string(), // Overlaps by 30bp
            corrected: "CGATCGATCGATCGATCGATCGATCGATCGAT".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 32],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
        CorrectedRead {
            id: 3,
            original: "GATCGATCGATCGATCGATCGATCGATCGATC".to_string(), // Overlaps by 29bp
            corrected: "GATCGATCGATCGATCGATCGATCGATCGATC".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 32],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        },
    ];

    println!(
        "📊 High coverage test data: {} overlapping reads",
        test_reads.len()
    );

    // Test with min_coverage=2 (should work with overlapping reads)
    println!("\n🔧 Testing high coverage data with min_coverage=2...");
    let builder = AdvancedAssemblyGraphBuilder::new(15, 31, 2, 8).unwrap();

    match builder.build_graph(&test_reads) {
        Ok(assembly_graph) => {
            println!("✅ Graph construction succeeded!");
            println!("📊 High coverage results:");
            println!(
                "  Graph fragment nodes: {}",
                assembly_graph.graph_fragment.nodes.len()
            );
            println!(
                "  Graph fragment edges: {}",
                assembly_graph.graph_fragment.edges.len()
            );
            println!("  Contigs: {}", assembly_graph.contigs.len());

            if assembly_graph.contigs.len() > 0 {
                println!("🎉 SUCCESS: High coverage dataset produces contigs!");

                for (i, contig) in assembly_graph.contigs.iter().enumerate() {
                    println!(
                        "  Contig {}: {}bp, coverage={:.1}x, sequence='{}'",
                        i + 1,
                        contig.length,
                        contig.coverage,
                        if contig.sequence.len() > 50 {
                            format!("{}...", &contig.sequence[..47])
                        } else {
                            contig.sequence.clone()
                        }
                    );
                }
            } else {
                println!("🚨 ISSUE: Even with high coverage, no contigs generated");
            }
        }
        Err(e) => {
            println!("❌ Graph construction failed: {}", e);
        }
    }
}
