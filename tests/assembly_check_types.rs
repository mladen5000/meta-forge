//! Quick type check test

use meta_forge::assembly::performance_optimizations::*;
use meta_forge::core::data_structures::*;

#[test]
fn test_contig_types() {
    println!("🧪 Testing contig types returned by ParallelContigGenerator");

    let mut graph = CacheOptimizedGraph::new(1);
    graph.add_node(123, 5);

    match ParallelContigGenerator::generate_contigs_parallel(&graph) {
        Ok(contigs) => {
            println!("✅ Successfully generated {} contigs", contigs.len());

            for (i, contig) in contigs.iter().enumerate() {
                // Check the actual fields available
                println!(
                    "Contig {}: id={}, length={}, coverage={:.1}, sequence='{}', type={:?}",
                    i + 1,
                    contig.id,
                    contig.length,
                    contig.coverage,
                    contig.sequence, // This should work if it's the right type
                    contig.contig_type
                );
            }
        }
        Err(e) => println!("❌ Failed: {}", e),
    }
}
