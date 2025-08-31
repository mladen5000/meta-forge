//! Test-Driven Development (TDD) Suite for Assembly Zero Contig Bug
//! ================================================================
//!
//! This comprehensive test suite systematically tests each component of the 
//! metagenomic assembly pipeline to identify why zero contigs are being generated.
//!
//! **TDD Approach:**
//! 1. Write failing tests that reveal the bug
//! 2. Fix the minimal code to make tests pass
//! 3. Refactor while keeping tests green
//! 4. Repeat until zero contig issue is resolved
//!
//! **Biological Context:**
//! - Metagenomic assembly should produce fewer contigs than input reads
//! - A 1:1 read:contig ratio indicates assembly failure
//! - Proper assembly connects overlapping k-mers into longer sequences
//! - Edge cases: short reads, low coverage, repetitive sequences, N bases

use crate::core::data_structures::*;
use crate::assembly::graph_construction::*;
use crate::assembly::performance_optimizations::*;
use anyhow::Result;
use std::collections::HashMap;

/// Test data generator for biologically realistic sequences
pub struct BiologicalTestData;

impl BiologicalTestData {
    /// Generate overlapping reads that should assembly into a single contig
    pub fn overlapping_reads() -> Vec<CorrectedRead> {
        vec![
            Self::create_read(0, "ATCGATCGATCGATCGATCG"), // 20bp read
            Self::create_read(1, "TCGATCGATCGATCGATCGA"), // Overlaps by 19bp
            Self::create_read(2, "CGATCGATCGATCGATCGAT"), // Overlaps by 19bp  
            Self::create_read(3, "GATCGATCGATCGATCGATC"), // Overlaps by 19bp
        ]
    }

    /// Generate reads with no overlaps (should produce multiple contigs)
    pub fn non_overlapping_reads() -> Vec<CorrectedRead> {
        vec![
            Self::create_read(0, "AAAAAAAAAAAAAAAAAAAA"), // 20bp A's
            Self::create_read(1, "CCCCCCCCCCCCCCCCCCCC"), // 20bp C's  
            Self::create_read(2, "GGGGGGGGGGGGGGGGGGGG"), // 20bp G's
            Self::create_read(3, "TTTTTTTTTTTTTTTTTTTT"), // 20bp T's
        ]
    }

    /// Generate reads with ambiguous bases (N's)
    pub fn ambiguous_base_reads() -> Vec<CorrectedRead> {
        vec![
            Self::create_read(0, "ATCGATCGANNNNATCGATC"), // Contains N's
            Self::create_read(1, "TCGATCGATCGATCGATCGA"), // Normal read
            Self::create_read(2, "NNNNNNNNNNNNNNNNNNN"), // All N's
        ]
    }

    /// Generate short reads (edge case)
    pub fn short_reads() -> Vec<CorrectedRead> {
        vec![
            Self::create_read(0, "ATG"),   // 3bp - too short for k=21
            Self::create_read(1, "CGA"),   // 3bp - too short for k=21  
            Self::create_read(2, "TCGA"),  // 4bp - too short for k=21
        ]
    }

    /// Generate low-coverage reads (single occurrence)
    pub fn low_coverage_reads() -> Vec<CorrectedRead> {
        vec![
            Self::create_read(0, "ATCGATCGATCGATCGATCGATCGATCGATCG"), // 32bp - single copy
        ]
    }

    /// Generate high-coverage reads (repetitive)
    pub fn high_coverage_reads() -> Vec<CorrectedRead> {
        let sequence = "ATCGATCGATCGATCGATCGATCGATCGATCG";
        (0..10).map(|i| Self::create_read(i, sequence)).collect()
    }

    /// Generate repetitive sequence reads
    pub fn repetitive_sequence_reads() -> Vec<CorrectedRead> {
        vec![
            Self::create_read(0, "ATATATATATATATATATATATATATATATA"), // AT repeat
            Self::create_read(1, "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"), // GC repeat
        ]
    }

    /// Create a CorrectedRead with realistic metadata
    fn create_read(id: usize, sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()], // High quality scores
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 3,
                correction_time_ms: 0,
            },
        }
    }
}

/// **TEST 1: K-mer Extraction and Counting**
/// This tests the fundamental step of breaking reads into k-mers
#[cfg(test)]
mod kmer_extraction_tests {
    use super::*;

    #[test]
    fn test_kmer_extraction_basic() {
        // FAILING TEST: Verify k-mer extraction works
        let read = BiologicalTestData::create_read(0, "ATCGATCG");
        let k = 4;
        
        // Extract k-mers manually to verify
        let expected_kmers = vec!["ATCG", "TCGA", "CGAT", "GATC", "ATCG"];
        
        // Test CanonicalKmer creation
        let mut canonical_kmers = Vec::new();
        for i in 0..=(read.corrected.len() - k) {
            let kmer_str = &read.corrected[i..i+k];
            match CanonicalKmer::new(kmer_str) {
                Ok(canonical) => canonical_kmers.push(canonical),
                Err(e) => panic!("Failed to create k-mer '{}': {}", kmer_str, e),
            }
        }
        
        assert_eq!(canonical_kmers.len(), expected_kmers.len(), 
                   "K-mer extraction should produce {} k-mers from {}-mer sequence with k={}",
                   expected_kmers.len(), read.corrected.len(), k);
        
        // Verify each k-mer is valid DNA
        for kmer in &canonical_kmers {
            assert!(kmer.sequence.chars().all(|c| matches!(c, 'A' | 'T' | 'C' | 'G')),
                   "All k-mers should contain only ATCG: {}", kmer.sequence);
        }
    }

    #[test]
    fn test_kmer_extraction_edge_cases() {
        // Test with short sequence
        let short_read = BiologicalTestData::create_read(0, "ATG");
        let k = 21; // Default k-mer size
        
        // Should handle short sequences gracefully
        let result = extract_kmers_from_read(&short_read, k);
        assert!(result.is_empty() || result.is_err(), 
               "Short reads ({}bp) should not produce k-mers with k={}", 
               short_read.corrected.len(), k);
    }

    #[test]
    fn test_kmer_extraction_ambiguous_bases() {
        let ambiguous_read = BiologicalTestData::create_read(0, "ATCGNNNATCG");
        let k = 4;
        
        // Should handle or skip k-mers with N's based on configuration
        let result = extract_kmers_from_read(&ambiguous_read, k);
        
        // This test should either skip N-containing k-mers or handle them appropriately
        if let Ok(kmers) = result {
            for kmer in kmers {
                // Each k-mer should either be valid DNA or contain N's if allowed
                assert!(kmer.sequence.chars().all(|c| matches!(c, 'A' | 'T' | 'C' | 'G' | 'N')),
                       "K-mer should contain only ATCGN: {}", kmer.sequence);
            }
        }
    }

    /// Extract k-mers from a read for testing
    fn extract_kmers_from_read(read: &CorrectedRead, k: usize) -> Result<Vec<CanonicalKmer>, String> {
        if read.corrected.len() < k {
            return Ok(Vec::new());
        }

        let mut kmers = Vec::new();
        for i in 0..=(read.corrected.len() - k) {
            let kmer_str = &read.corrected[i..i+k];
            match CanonicalKmer::new(kmer_str) {
                Ok(canonical) => kmers.push(canonical),
                Err(e) => return Err(format!("Failed to create k-mer '{}': {}", kmer_str, e)),
            }
        }
        Ok(kmers)
    }
}

/// **TEST 2: Graph Construction from K-mers**
/// This tests that k-mers are properly connected in the assembly graph
#[cfg(test)]
mod graph_construction_tests {
    use super::*;

    #[test]
    fn test_graph_construction_single_read() {
        // FAILING TEST: Single read should create connected nodes
        let read = BiologicalTestData::create_read(0, "ATCGATCG");
        let k = 4;
        
        let mut chunk = AssemblyChunk::new(0, k);
        chunk.add_read(read).unwrap();
        chunk.finalize();
        
        // Verify nodes were created
        assert!(!chunk.graph_fragment.nodes.is_empty(), 
               "Graph should contain nodes from k-mer extraction");
        
        // Verify edges were created between consecutive k-mers
        assert!(!chunk.graph_fragment.edges.is_empty(),
               "Graph should contain edges connecting consecutive k-mers");
        
        let node_count = chunk.graph_fragment.nodes.len();
        let edge_count = chunk.graph_fragment.edges.len();
        
        // For a single read, we expect (read_length - k + 1) nodes
        // and (read_length - k) edges connecting them
        let expected_nodes = 8 - k + 1; // = 5 nodes
        let expected_edges = expected_nodes - 1; // = 4 edges
        
        assert_eq!(node_count, expected_nodes, 
                   "Should have {} nodes for {}-bp read with k={}", 
                   expected_nodes, 8, k);
        assert_eq!(edge_count, expected_edges,
                   "Should have {} edges for {} nodes in linear chain",
                   expected_edges, node_count);
    }

    #[test] 
    fn test_graph_construction_overlapping_reads() {
        // FAILING TEST: Overlapping reads should merge nodes
        let reads = BiologicalTestData::overlapping_reads();
        let k = 4;
        
        let mut chunk = AssemblyChunk::new(0, k);
        for read in reads {
            chunk.add_read(read).unwrap();
        }
        chunk.finalize();
        
        // With overlapping reads, we expect fewer unique nodes than total k-mers
        let unique_nodes = chunk.graph_fragment.nodes.len();
        assert!(unique_nodes > 0, "Should have created nodes from overlapping reads");
        
        // Verify that overlapping k-mers increased coverage
        let mut high_coverage_nodes = 0;
        for node in chunk.graph_fragment.nodes.values() {
            if node.coverage > 1 {
                high_coverage_nodes += 1;
            }
        }
        
        assert!(high_coverage_nodes > 0, 
               "Overlapping reads should create high-coverage nodes");
    }

    #[test]
    fn test_graph_construction_non_overlapping_reads() {
        // FAILING TEST: Non-overlapping reads should create separate components
        let reads = BiologicalTestData::non_overlapping_reads();
        let k = 4;
        
        let mut chunk = AssemblyChunk::new(0, k);
        for read in reads {
            chunk.add_read(read).unwrap();
        }
        chunk.finalize();
        
        // Each non-overlapping read should contribute its own nodes
        let node_count = chunk.graph_fragment.nodes.len();
        assert!(node_count > 0, "Should have created nodes from non-overlapping reads");
        
        // Verify connectivity - should have multiple disconnected components
        let connected_components = find_connected_components(&chunk.graph_fragment);
        assert!(connected_components.len() > 1, 
               "Non-overlapping reads should create multiple connected components");
    }

    /// Find connected components in a graph fragment (helper function)
    fn find_connected_components(fragment: &GraphFragment) -> Vec<Vec<u64>> {
        let mut visited = std::collections::HashSet::new();
        let mut components = Vec::new();
        
        for &node_hash in fragment.nodes.keys() {
            if !visited.contains(&node_hash) {
                let mut component = Vec::new();
                dfs_component(node_hash, fragment, &mut visited, &mut component);
                if !component.is_empty() {
                    components.push(component);
                }
            }
        }
        
        components
    }
    
    fn dfs_component(node_hash: u64, fragment: &GraphFragment, 
                     visited: &mut std::collections::HashSet<u64>, 
                     component: &mut Vec<u64>) {
        if visited.contains(&node_hash) {
            return;
        }
        
        visited.insert(node_hash);
        component.push(node_hash);
        
        // Find neighbors through edges
        for edge in &fragment.edges {
            if edge.from_hash == node_hash && !visited.contains(&edge.to_hash) {
                dfs_component(edge.to_hash, fragment, visited, component);
            }
            if edge.to_hash == node_hash && !visited.contains(&edge.from_hash) {
                dfs_component(edge.from_hash, fragment, visited, component);
            }
        }
    }
}

/// **TEST 3: Path Traversal and Contig Generation**
/// This tests the critical step where graphs are converted to contigs
#[cfg(test)]
mod contig_generation_tests {
    use super::*;

    #[test]
    fn test_contig_generation_simple_path() {
        // FAILING TEST: Simple linear path should generate one contig
        let reads = BiologicalTestData::overlapping_reads();
        let k = 4;
        
        // Build graph
        let mut chunk = AssemblyChunk::new(0, k);
        for read in reads {
            chunk.add_read(read).unwrap();
        }
        chunk.finalize();
        
        // Test path traversal
        let paths = find_eulerian_paths(&chunk.graph_fragment);
        assert!(!paths.is_empty(), "Should find at least one path through the graph");
        
        // Convert paths to contigs
        let mut contigs = Vec::new();
        for (i, path) in paths.iter().enumerate() {
            match path_to_contig(i, path, &chunk.graph_fragment) {
                Ok(contig) => contigs.push(contig),
                Err(e) => panic!("Failed to convert path to contig: {}", e),
            }
        }
        
        assert!(!contigs.is_empty(), "Should generate at least one contig");
        
        // Verify contig properties
        for contig in &contigs {
            assert!(contig.length > 0, "Contig should have positive length");
            assert!(!contig.sequence.is_empty(), "Contig should have sequence");
            assert!(contig.coverage > 0.0, "Contig should have positive coverage");
        }
    }

    #[test]
    fn test_contig_generation_multiple_components() {
        // FAILING TEST: Multiple components should generate multiple contigs
        let reads = BiologicalTestData::non_overlapping_reads();
        let k = 4;
        
        // Build graph
        let mut chunk = AssemblyChunk::new(0, k);
        for read in reads {
            chunk.add_read(read).unwrap();
        }
        chunk.finalize();
        
        // Use ParallelContigGenerator to generate contigs
        let cache_graph = convert_to_cache_optimized_graph(&chunk.graph_fragment);
        
        match ParallelContigGenerator::generate_contigs_parallel(&cache_graph) {
            Ok(contigs) => {
                assert!(!contigs.is_empty(), "Should generate contigs from graph");
                
                // With non-overlapping reads, expect multiple contigs
                assert!(contigs.len() > 1, 
                       "Non-overlapping reads should produce multiple contigs");
                
                // Verify each contig has valid properties
                for contig in &contigs {
                    assert!(contig.length > 0, "Each contig should have positive length");
                    assert!(!contig.sequence.is_empty(), "Each contig should have sequence");
                }
            }
            Err(e) => panic!("Contig generation failed: {}", e),
        }
    }

    #[test]
    fn test_minimum_contig_length_filtering() {
        // FAILING TEST: Short contigs should be filtered out
        let reads = BiologicalTestData::short_reads();
        let k = 4;
        
        // Build graph from short reads
        let mut chunk = AssemblyChunk::new(0, k);
        for read in reads {
            if let Ok(()) = chunk.add_read(read) {
                // Successfully added read
            }
        }
        chunk.finalize();
        
        // Convert and generate contigs
        let cache_graph = convert_to_cache_optimized_graph(&chunk.graph_fragment);
        
        match ParallelContigGenerator::generate_contigs_parallel(&cache_graph) {
            Ok(contigs) => {
                // Filter by minimum length (e.g., 50 bp)
                let min_length = 50;
                let filtered_contigs: Vec<_> = contigs.into_iter()
                    .filter(|c| c.length >= min_length)
                    .collect();
                
                // Short reads should not produce long contigs
                assert_eq!(filtered_contigs.len(), 0,
                          "Short reads should not produce contigs >= {} bp", min_length);
            }
            Err(_) => {
                // Expected for very short reads
                println!("Short reads correctly failed to generate contigs");
            }
        }
    }

    /// Helper: Find Eulerian paths in a graph fragment
    fn find_eulerian_paths(fragment: &GraphFragment) -> Vec<Vec<u64>> {
        // Simplified path finding - just return all possible linear chains
        let mut paths = Vec::new();
        let mut visited_edges = std::collections::HashSet::new();
        
        // Find start nodes (nodes with in-degree = 0 or out-degree > in-degree)
        let adjacency = fragment.get_adjacency_list();
        let mut in_degrees = std::collections::HashMap::new();
        let mut out_degrees = std::collections::HashMap::new();
        
        // Calculate degrees
        for (&from, neighbors) in &adjacency {
            out_degrees.insert(from, neighbors.len());
            for &to in neighbors {
                *in_degrees.entry(to).or_insert(0) += 1;
            }
        }
        
        // Find potential start nodes
        for &node in fragment.nodes.keys() {
            let in_deg = in_degrees.get(&node).copied().unwrap_or(0);
            let out_deg = out_degrees.get(&node).copied().unwrap_or(0);
            
            if in_deg == 0 || out_deg > in_deg {
                // Try to find a path starting from this node
                let mut path = Vec::new();
                dfs_path(node, &adjacency, &mut visited_edges, &mut path);
                if !path.is_empty() {
                    paths.push(path);
                }
            }
        }
        
        paths
    }
    
    fn dfs_path(node: u64, adjacency: &std::collections::HashMap<u64, Vec<u64>>, 
                visited_edges: &mut std::collections::HashSet<(u64, u64)>,
                path: &mut Vec<u64>) {
        path.push(node);
        
        if let Some(neighbors) = adjacency.get(&node) {
            for &neighbor in neighbors {
                let edge = (node, neighbor);
                if !visited_edges.contains(&edge) {
                    visited_edges.insert(edge);
                    dfs_path(neighbor, adjacency, visited_edges, path);
                    return; // Follow only one path for simplicity
                }
            }
        }
    }
    
    /// Helper: Convert path to contig
    fn path_to_contig(id: usize, path: &[u64], fragment: &GraphFragment) -> Result<Contig> {
        if path.is_empty() {
            return Err(anyhow!("Empty path cannot be converted to contig"));
        }
        
        let sequence = fragment.reconstruct_sequence_from_path(path)?;
        let coverage = fragment.calculate_path_coverage_from_hashes(path);
        
        Ok(Contig {
            id,
            sequence: sequence.clone(),
            coverage,
            length: sequence.len(),
            node_path: path.to_vec(),
            contig_type: ContigType::Linear,
        })
    }
    
    /// Helper: Convert GraphFragment to CacheOptimizedGraph
    fn convert_to_cache_optimized_graph(fragment: &GraphFragment) -> CacheOptimizedGraph {
        let mut cache_graph = CacheOptimizedGraph::new(fragment.nodes.len());
        
        // Add all nodes
        for node in fragment.nodes.values() {
            cache_graph.add_node(node.kmer.hash, node.coverage);
        }
        
        // Add all edges
        for edge in &fragment.edges {
            let _ = cache_graph.add_edge(edge.from_hash, edge.to_hash);
        }
        
        cache_graph
    }
}

/// **TEST 4: Integration Tests - Full Pipeline**
/// These tests exercise the complete assembly pipeline
#[cfg(test)]
mod integration_tests {
    use super::*;

    #[test] 
    fn test_full_assembly_pipeline_overlapping() {
        // INTEGRATION TEST: Full pipeline with overlapping reads
        let reads = BiologicalTestData::overlapping_reads();
        let k = 4;
        
        println!("🧬 Testing full assembly pipeline with {} overlapping reads", reads.len());
        
        // Step 1: Build assembly chunks
        let mut chunk = AssemblyChunk::new(0, k);
        for read in &reads {
            chunk.add_read(read.clone()).unwrap();
        }
        chunk.finalize();
        
        println!("📊 Graph fragment: {} nodes, {} edges", 
                chunk.graph_fragment.nodes.len(),
                chunk.graph_fragment.edges.len());
        
        assert!(!chunk.graph_fragment.nodes.is_empty(), "Should create nodes");
        assert!(!chunk.graph_fragment.edges.is_empty(), "Should create edges");
        
        // Step 2: Convert to optimized graph
        let cache_graph = convert_to_cache_optimized_graph(&chunk.graph_fragment);
        let (nodes, edges, _, _) = cache_graph.get_statistics();
        
        println!("⚡ Cache graph: {} nodes, {} edges", nodes, edges);
        assert_eq!(nodes, chunk.graph_fragment.nodes.len());
        
        // Step 3: Generate contigs
        match ParallelContigGenerator::generate_contigs_parallel(&cache_graph) {
            Ok(contigs) => {
                println!("✅ Generated {} contigs", contigs.len());
                
                // Critical assertion: should produce fewer contigs than reads
                assert!(contigs.len() < reads.len(),
                       "Assembly should produce fewer contigs ({}) than reads ({})", 
                       contigs.len(), reads.len());
                
                // Verify contig properties
                for (i, contig) in contigs.iter().enumerate() {
                    println!("  Contig {}: {}bp, coverage={:.1}x", 
                            i, contig.length, contig.coverage);
                    
                    assert!(contig.length > 0, "Contig should have positive length");
                    assert!(!contig.sequence.is_empty(), "Contig should have sequence");
                    assert!(contig.coverage > 0.0, "Contig should have coverage");
                }
                
                // Check for proper assembly (longer contigs than reads)
                let max_contig_length = contigs.iter().map(|c| c.length).max().unwrap_or(0);
                let max_read_length = reads.iter().map(|r| r.corrected.len()).max().unwrap_or(0);
                
                assert!(max_contig_length >= max_read_length,
                       "At least one contig should be as long as the longest read");
            }
            Err(e) => panic!("Full pipeline failed at contig generation: {}", e),
        }
    }

    #[test]
    fn test_full_assembly_pipeline_edge_cases() {
        // INTEGRATION TEST: Edge cases that might cause zero contigs
        let test_cases = vec![
            ("Short reads", BiologicalTestData::short_reads()),
            ("Ambiguous bases", BiologicalTestData::ambiguous_base_reads()),
            ("Low coverage", BiologicalTestData::low_coverage_reads()),
        ];
        
        let k = 4; // Use smaller k for short reads
        
        for (name, reads) in test_cases {
            println!("🧪 Testing edge case: {}", name);
            
            let mut chunk = AssemblyChunk::new(0, k);
            for read in &reads {
                if let Ok(()) = chunk.add_read(read.clone()) {
                    // Successfully added
                }
            }
            chunk.finalize();
            
            println!("  Graph: {} nodes, {} edges", 
                    chunk.graph_fragment.nodes.len(),
                    chunk.graph_fragment.edges.len());
            
            if !chunk.graph_fragment.nodes.is_empty() {
                let cache_graph = convert_to_cache_optimized_graph(&chunk.graph_fragment);
                
                match ParallelContigGenerator::generate_contigs_parallel(&cache_graph) {
                    Ok(contigs) => {
                        println!("  ✅ Generated {} contigs", contigs.len());
                        
                        if contigs.is_empty() {
                            println!("  ⚠️  Zero contigs generated - this might be expected for {}", name);
                        }
                    }
                    Err(e) => {
                        println!("  ❌ Contig generation failed: {}", e);
                    }
                }
            } else {
                println!("  ⚠️  No graph nodes created - likely due to filtering");
            }
        }
    }

    /// Helper function reused from contig_generation_tests
    fn convert_to_cache_optimized_graph(fragment: &GraphFragment) -> CacheOptimizedGraph {
        let mut cache_graph = CacheOptimizedGraph::new(fragment.nodes.len());
        
        // Add all nodes
        for node in fragment.nodes.values() {
            cache_graph.add_node(node.kmer.hash, node.coverage);
        }
        
        // Add all edges  
        for edge in &fragment.edges {
            let _ = cache_graph.add_edge(edge.from_hash, edge.to_hash);
        }
        
        cache_graph
    }
}

/// **TEST 5: Property-Based Tests**
/// Tests that verify invariants across different data inputs
#[cfg(test)]
mod property_tests {
    use super::*;

    #[test]
    fn test_assembly_conservation_properties() {
        // PROPERTY TEST: Assembly should conserve biological properties
        let reads = BiologicalTestData::high_coverage_reads();
        let k = 4;
        
        // Calculate input properties
        let total_input_length: usize = reads.iter().map(|r| r.corrected.len()).sum();
        let input_gc_content = calculate_input_gc_content(&reads);
        
        // Run assembly
        let mut chunk = AssemblyChunk::new(0, k);
        for read in &reads {
            chunk.add_read(read.clone()).unwrap();
        }
        chunk.finalize();
        
        let cache_graph = convert_to_cache_optimized_graph(&chunk.graph_fragment);
        
        match ParallelContigGenerator::generate_contigs_parallel(&cache_graph) {
            Ok(contigs) => {
                if !contigs.is_empty() {
                    // Calculate output properties
                    let total_output_length: usize = contigs.iter().map(|c| c.length).sum();
                    let output_gc_content = calculate_output_gc_content(&contigs);
                    
                    // Property: GC content should be approximately conserved
                    let gc_diff = (input_gc_content - output_gc_content).abs();
                    assert!(gc_diff < 0.1, 
                           "GC content should be conserved: input={:.2}, output={:.2}", 
                           input_gc_content, output_gc_content);
                    
                    // Property: Output length should not exceed input (except for overlaps)
                    assert!(total_output_length <= total_input_length + (reads.len() * k),
                           "Output length should not greatly exceed input length");
                    
                    println!("✅ Conservation properties maintained:");
                    println!("  Input: {}bp, GC={:.2}%", total_input_length, input_gc_content * 100.0);
                    println!("  Output: {}bp, GC={:.2}%", total_output_length, output_gc_content * 100.0);
                } else {
                    panic!("Property test failed: zero contigs generated");
                }
            }
            Err(e) => panic!("Property test failed: {}", e),
        }
    }

    #[test]
    fn test_assembly_scalability_properties() {
        // PROPERTY TEST: Assembly should scale reasonably with input size
        let sizes = vec![10, 50, 100];
        let k = 4;
        
        for &size in &sizes {
            println!("🔍 Testing scalability with {} reads", size);
            
            // Generate test data
            let reads: Vec<_> = (0..size).map(|i| {
                BiologicalTestData::create_read(i, "ATCGATCGATCGATCGATCG")
            }).collect();
            
            let start = std::time::Instant::now();
            
            // Run assembly
            let mut chunk = AssemblyChunk::new(0, k);
            for read in &reads {
                chunk.add_read(read.clone()).unwrap();
            }
            chunk.finalize();
            
            let cache_graph = convert_to_cache_optimized_graph(&chunk.graph_fragment);
            
            match ParallelContigGenerator::generate_contigs_parallel(&cache_graph) {
                Ok(contigs) => {
                    let elapsed = start.elapsed();
                    let reads_per_sec = size as f64 / elapsed.as_secs_f64();
                    
                    println!("  ✅ {} reads → {} contigs in {:.3}s ({:.0} reads/sec)",
                            size, contigs.len(), elapsed.as_secs_f64(), reads_per_sec);
                    
                    // Property: Should maintain reasonable performance
                    assert!(elapsed.as_secs() < 10, 
                           "Assembly should complete within 10 seconds for {} reads", size);
                }
                Err(e) => println!("  ❌ Failed for size {}: {}", size, e),
            }
        }
    }

    /// Helper functions
    fn calculate_input_gc_content(reads: &[CorrectedRead]) -> f64 {
        let mut g_count = 0;
        let mut c_count = 0;
        let mut total = 0;
        
        for read in reads {
            for c in read.corrected.chars() {
                match c.to_ascii_uppercase() {
                    'G' => { g_count += 1; total += 1; }
                    'C' => { c_count += 1; total += 1; }
                    'A' | 'T' => { total += 1; }
                    _ => {} // Skip non-DNA characters
                }
            }
        }
        
        if total > 0 {
            (g_count + c_count) as f64 / total as f64
        } else {
            0.0
        }
    }
    
    fn calculate_output_gc_content(contigs: &[Contig]) -> f64 {
        let mut g_count = 0;
        let mut c_count = 0;
        let mut total = 0;
        
        for contig in contigs {
            for c in contig.sequence.chars() {
                match c.to_ascii_uppercase() {
                    'G' => { g_count += 1; total += 1; }
                    'C' => { c_count += 1; total += 1; }
                    'A' | 'T' => { total += 1; }
                    _ => {} // Skip non-DNA characters
                }
            }
        }
        
        if total > 0 {
            (g_count + c_count) as f64 / total as f64
        } else {
            0.0
        }
    }

    fn convert_to_cache_optimized_graph(fragment: &GraphFragment) -> CacheOptimizedGraph {
        let mut cache_graph = CacheOptimizedGraph::new(fragment.nodes.len());
        
        // Add all nodes
        for node in fragment.nodes.values() {
            cache_graph.add_node(node.kmer.hash, node.coverage);
        }
        
        // Add all edges
        for edge in &fragment.edges {
            let _ = cache_graph.add_edge(edge.from_hash, edge.to_hash);
        }
        
        cache_graph
    }
}

/// **DEBUGGING UTILITIES**
/// Functions to help diagnose the zero contig issue
pub mod debug_utils {
    use super::*;

    /// Print detailed information about a graph fragment
    pub fn debug_graph_fragment(fragment: &GraphFragment, name: &str) {
        println!("\n🔍 DEBUG: GraphFragment '{}' Analysis", name);
        println!("{'=':.50}", "");
        println!("Nodes: {}", fragment.nodes.len());
        println!("Edges: {}", fragment.edges.len());
        
        if fragment.nodes.is_empty() {
            println!("⚠️  NO NODES FOUND - This is likely the bug!");
            return;
        }
        
        // Analyze nodes
        println!("\n📊 Node Analysis:");
        let mut coverage_dist = HashMap::new();
        for (hash, node) in &fragment.nodes {
            *coverage_dist.entry(node.coverage).or_insert(0) += 1;
            if fragment.nodes.len() <= 10 {
                println!("  Node {}: coverage={}, type={:?}", 
                        hash, node.coverage, node.node_type);
            }
        }
        
        println!("Coverage distribution: {:?}", coverage_dist);
        
        // Analyze edges
        println!("\n🔗 Edge Analysis:");
        if fragment.edges.is_empty() {
            println!("⚠️  NO EDGES FOUND - Graph nodes are disconnected!");
        } else {
            let mut weight_dist = HashMap::new();
            for edge in &fragment.edges {
                *weight_dist.entry(edge.weight).or_insert(0) += 1;
                if fragment.edges.len() <= 10 {
                    println!("  Edge: {} → {} (weight={})", 
                            edge.from_hash, edge.to_hash, edge.weight);
                }
            }
            println!("Weight distribution: {:?}", weight_dist);
        }
        
        // Connectivity analysis
        let components = find_connected_components(fragment);
        println!("\n🔗 Connectivity Analysis:");
        println!("Connected components: {}", components.len());
        for (i, component) in components.iter().enumerate() {
            println!("  Component {}: {} nodes", i, component.len());
        }
        
        if components.len() == fragment.nodes.len() {
            println!("⚠️  ALL NODES ARE ISOLATED - This explains zero contigs!");
        }
    }
    
    /// Debug the ParallelContigGenerator process
    pub fn debug_contig_generation(cache_graph: &CacheOptimizedGraph, name: &str) {
        println!("\n🧬 DEBUG: ContigGeneration '{}' Analysis", name);
        println!("{'=':.50}", "");
        
        let (nodes, edges, memory, cache_rate) = cache_graph.get_statistics();
        println!("CacheOptimizedGraph: {} nodes, {} edges", nodes, edges);
        println!("Memory usage: {} bytes, Cache hit rate: {:.1}%", memory, cache_rate * 100.0);
        
        if nodes == 0 {
            println!("⚠️  NO NODES IN CACHE GRAPH - Cannot generate contigs!");
            return;
        }
        
        // Check individual nodes
        let node_ids = cache_graph.get_node_ids();
        println!("Node IDs: {:?}", &node_ids[..node_ids.len().min(10)]);
        
        // Check connectivity
        for &node_id in &node_ids[..node_ids.len().min(5)] {
            let neighbors = cache_graph.get_neighbors(node_id);
            println!("Node {} has {} neighbors: {:?}", 
                    node_id, neighbors.len(), &neighbors[..neighbors.len().min(3)]);
        }
    }
    
    fn find_connected_components(fragment: &GraphFragment) -> Vec<Vec<u64>> {
        let mut visited = std::collections::HashSet::new();
        let mut components = Vec::new();
        
        for &node_hash in fragment.nodes.keys() {
            if !visited.contains(&node_hash) {
                let mut component = Vec::new();
                dfs_component(node_hash, fragment, &mut visited, &mut component);
                if !component.is_empty() {
                    components.push(component);
                }
            }
        }
        
        components
    }
    
    fn dfs_component(node_hash: u64, fragment: &GraphFragment, 
                     visited: &mut std::collections::HashSet<u64>, 
                     component: &mut Vec<u64>) {
        if visited.contains(&node_hash) {
            return;
        }
        
        visited.insert(node_hash);
        component.push(node_hash);
        
        // Find neighbors through edges
        for edge in &fragment.edges {
            if edge.from_hash == node_hash && !visited.contains(&edge.to_hash) {
                dfs_component(edge.to_hash, fragment, visited, component);
            }
            if edge.to_hash == node_hash && !visited.contains(&edge.from_hash) {
                dfs_component(edge.from_hash, fragment, visited, component);
            }
        }
    }
}