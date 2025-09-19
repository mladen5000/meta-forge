//! Debug test specifically for AssemblyChunk processing
//! This will help identify why k-mers are extracted but not converted to graph nodes.

use meta_forge::assembly::graph_construction::AssemblyChunk;
use meta_forge::core::data_structures::*;

#[test]
fn test_assembly_chunk_finalize() {
    println!("🧪 Testing AssemblyChunk finalize process");

    // Create test read
    let test_read = CorrectedRead {
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
    };

    let k = 15; // Use same k as pipeline
    println!(
        "📖 Input read: '{}' ({}bp)",
        test_read.corrected,
        test_read.corrected.len()
    );
    println!("🔢 K-mer size: {}", k);

    // Test k-mer extraction manually
    if test_read.corrected.len() >= k {
        let expected_kmers = test_read.corrected.len() - k + 1;
        println!("📊 Expected k-mers: {}", expected_kmers);

        // Manual k-mer extraction to verify
        for i in 0..=test_read.corrected.len() - k {
            let kmer_str = &test_read.corrected[i..i + k];
            match CanonicalKmer::new(kmer_str) {
                Ok(canonical) => {
                    println!(
                        "  K-mer {}: '{}' -> canonical '{}' (hash={})",
                        i + 1,
                        kmer_str,
                        canonical.sequence,
                        canonical.hash
                    );
                }
                Err(e) => {
                    println!("  ❌ Failed to create k-mer '{}': {}", kmer_str, e);
                }
            }
        }
    }

    // Create and test AssemblyChunk
    let mut chunk = AssemblyChunk::new(0, k);

    println!("\n🔧 Testing AssemblyChunk processing...");
    println!("Before adding read:");
    println!(
        "  GraphFragment nodes: {}",
        chunk.graph_fragment.nodes.len()
    );
    println!(
        "  GraphFragment edges: {}",
        chunk.graph_fragment.edges.len()
    );

    // Add read
    match chunk.add_read(test_read) {
        Ok(()) => println!("✅ Read added successfully"),
        Err(e) => {
            println!("❌ Failed to add read: {}", e);
            panic!("Read addition failed");
        }
    }

    println!("After adding read (before finalize):");
    println!(
        "  GraphFragment nodes: {}",
        chunk.graph_fragment.nodes.len()
    );
    println!(
        "  GraphFragment edges: {}",
        chunk.graph_fragment.edges.len()
    );

    // Finalize the chunk
    println!("\n🏁 Finalizing chunk...");
    chunk.finalize();

    println!("After finalize:");
    println!(
        "  GraphFragment nodes: {}",
        chunk.graph_fragment.nodes.len()
    );
    println!(
        "  GraphFragment edges: {}",
        chunk.graph_fragment.edges.len()
    );

    if chunk.graph_fragment.nodes.is_empty() {
        println!("🚨 BUG FOUND: AssemblyChunk.finalize() created no nodes!");
        println!("This explains why the pipeline produces zero contigs.");

        // Let's look at the reads in the chunk
        println!("\nDEBUG: Chunk contains {} reads:", chunk.reads.len());
        for (i, read) in chunk.reads.iter().enumerate() {
            println!(
                "  Read {}: '{}' ({}bp)",
                i,
                read.corrected,
                read.corrected.len()
            );
        }

        panic!("AssemblyChunk.finalize() is not working correctly!");
    } else {
        println!(
            "✅ AssemblyChunk.finalize() created {} nodes and {} edges",
            chunk.graph_fragment.nodes.len(),
            chunk.graph_fragment.edges.len()
        );

        // Print node details
        for (hash, node) in &chunk.graph_fragment.nodes {
            println!(
                "  Node {}: k-mer='{}', coverage={}",
                hash, node.kmer.sequence, node.coverage
            );
        }

        // Print edge details
        for (i, edge) in chunk.graph_fragment.edges.iter().enumerate() {
            println!(
                "  Edge {}: {} -> {} (weight={})",
                i + 1,
                edge.from_hash,
                edge.to_hash,
                edge.weight
            );
        }
    }
}
