# Additional Test Coverage Recommendations

## Biological Edge Cases Requiring Test Coverage

### 1. DNA Sequence Edge Cases

#### Ambiguous Base Handling
```rust
#[test]
fn test_ambiguous_base_iupac_codes() {
    // Test IUPAC ambiguous nucleotide codes
    let sequences = [
        "ATCGRYWSMKHVDBN",  // All IUPAC codes
        "ATCGATNNNGATCG",    // Multiple N's (common in real data)
        "NNNNNNATCGATCG",    // N's at start (5' end artifacts)  
        "ATCGATCGNNNNNNN",   // N's at end (3' end artifacts)
    ];
    
    for seq in sequences {
        let result = validator.validate_sequence(seq, None);
        // Should handle gracefully based on N content threshold
    }
}
```

#### Homopolymer Runs (Sequencing Artifacts)
```rust
#[test] 
fn test_homopolymer_edge_cases() {
    let test_cases = [
        ("AAAAAAAAAA", "poly-A tract"),           // 10 A's
        ("TTTTTTTTTTTTTT", "poly-T tract"),      // 14 T's  
        ("GGGGGGGGGGGGGGGG", "poly-G tract"),    // 16 G's (harder to sequence)
        ("CCCCCCCCCCCCCCCCCC", "poly-C tract"),  // 18 C's (harder to sequence)
    ];
    
    for (seq, desc) in test_cases {
        // Test both in isolation and embedded in longer sequences
        let embedded = format!("ATCGATCG{}GATCGATC", seq);
        // Should detect but handle appropriately
    }
}
```

#### Extreme GC Content
```rust
#[test]
fn test_extreme_gc_content_organisms() {
    let test_cases = [
        ("ATATATATATATATATATATATATATATATATATAT", 0.0, "AT-rich (Plasmodium)"),
        ("GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC", 1.0, "GC-rich (Streptomyces)"),
        ("AAAAAATCGATCGAAAAAATCGATCGAAAAAA", 0.25, "Low GC (some archaea)"),
        ("GCCCGGGCCCGGGCCCGGGCCCGGGCCCGGG", 0.8, "High GC (Mycobacterium)"),
    ];
    
    for (seq, expected_gc, organism) in test_cases {
        let result = validator.validate_sequence(seq, None);
        assert!((result.metrics.gc_content - expected_gc).abs() < 0.05);
    }
}
```

### 2. Assembly Graph Algorithm Edge Cases

#### Short Read Overlap Detection
```rust
#[test]
fn test_minimal_overlap_detection() {
    // Test minimum overlap requirements for assembly
    let reads = [
        "ATCGATCGATCGATCG",           // Read 1
        "GATCGATCGATCGAAG",           // Read 2 (12bp overlap)
        "TCGATCGATCGAAGCC",           // Read 3 (10bp overlap) 
        "GATCGATCGAAGCCTT",           // Read 4 (8bp overlap)
    ];
    
    // Test various minimum overlap thresholds
    for min_overlap in [6, 8, 10, 12] {
        let builder = AssemblyGraphBuilder::new(15, min_overlap, 1);
        let result = builder.build(&create_reads(&reads));
        // Verify graph connectivity vs overlap threshold
    }
}
```

#### Repetitive Element Handling
```rust
#[test]
fn test_repetitive_element_assembly() {
    // Simulate transposable elements, tandem repeats
    let repeat_unit = "ATCGATCG";
    let sequences = [
        format!("GCTAGCTA{}{}{}{}", repeat_unit, repeat_unit, repeat_unit, "CGATCGAT"),
        format!("AGTCAGTC{}{}{}{}", repeat_unit, repeat_unit, repeat_unit, "TAGCTAGC"), 
        format!("TTCGTTCG{}{}{}{}", repeat_unit, repeat_unit, repeat_unit, "AAGCAAGC"),
    ];
    
    // Assembly should detect but not get confused by repeats
    let result = build_assembly(&sequences);
    
    // Check that repetitive regions are properly marked
    let repetitive_nodes = result.graph_fragment.nodes
        .values()
        .filter(|n| matches!(n.node_type, NodeType::Repetitive))
        .count();
    
    assert!(repetitive_nodes > 0, "Should detect repetitive elements");
}
```

#### Coverage Distribution Edge Cases
```rust
#[test]
fn test_variable_coverage_assembly() {
    let mut reads = Vec::new();
    
    // Simulate realistic metagenomic coverage distribution
    // High abundance species (100x coverage)
    for i in 0..100 {
        reads.push(create_read(i, "ATCGATCGATCGATCGATCGATCGATCGATCG", 35));
    }
    
    // Medium abundance species (10x coverage) 
    for i in 100..110 {
        reads.push(create_read(i, "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA", 30));
    }
    
    // Low abundance species (2x coverage)
    for i in 110..112 {
        reads.push(create_read(i, "TTAACCGGTTAACCGGTTAACCGGTTAACCGG", 25));
    }
    
    let result = build_assembly(&reads);
    
    // Verify assembly handles coverage variations appropriately
    let coverage_stats = &result.graph_fragment.coverage_stats;
    assert!(coverage_stats.mean_coverage > 5.0);
    assert!(coverage_stats.median_coverage > 0);
    
    // Check that low-coverage sequences aren't completely filtered out
    let low_cov_nodes = result.graph_fragment.nodes
        .values()
        .filter(|n| n.coverage <= 5)
        .count();
    
    assert!(low_cov_nodes > 0, "Should retain some low-coverage sequences");
}
```

### 3. Machine Learning Feature Edge Cases

#### Graph Topology Extremes
```rust
#[test]
fn test_extreme_graph_topologies() {
    let test_cases = [
        ("linear", create_linear_graph(10)),           // No branches
        ("star", create_star_graph(10)),               // One central hub
        ("complete", create_complete_graph(5)),        // All nodes connected
        ("disconnected", create_disconnected_graph()), // Multiple components
    ];
    
    for (topology_name, graph) in test_cases {
        let features = extract_graph_features(&graph);
        
        // Features should capture topological properties
        match topology_name {
            "linear" => assert!(features.average_degree() < 2.5),
            "star" => assert!(features.max_degree() > features.average_degree() * 2.0),
            "complete" => assert!(features.clustering_coefficient() > 0.8),
            "disconnected" => assert!(features.connected_components() > 1),
            _ => {}
        }
    }
}
```

#### Feature Extraction Robustness
```rust
#[test]
fn test_feature_extraction_robustness() {
    let edge_cases = [
        create_empty_graph(),                    // No nodes
        create_single_node_graph(),              // Isolated node
        create_self_loop_graph(),                // Self-referencing edges
        create_highly_connected_graph(1000),    // Large graph
    ];
    
    for graph in edge_cases {
        let result = std::panic::catch_unwind(|| {
            extract_node_features(&graph)
        });
        
        // Should not panic on edge cases
        assert!(result.is_ok(), "Feature extraction should handle edge cases");
    }
}
```

### 4. Database Integration Edge Cases

#### Taxonomic Hierarchy Validation
```rust
#[test]
fn test_taxonomic_hierarchy_consistency() {
    let taxonomy_entries = [
        TaxonomyEntry { id: 1, name: "Root", parent_id: None, rank: "no rank" },
        TaxonomyEntry { id: 2, name: "Bacteria", parent_id: Some(1), rank: "superkingdom" },
        TaxonomyEntry { id: 3, name: "Proteobacteria", parent_id: Some(2), rank: "phylum" },
        TaxonomyEntry { id: 4, name: "Gammaproteobacteria", parent_id: Some(3), rank: "class" },
        TaxonomyEntry { id: 5, name: "Enterobacterales", parent_id: Some(4), rank: "order" },
    ];
    
    db.insert_taxonomy_entries(&taxonomy_entries).unwrap();
    
    // Test hierarchy traversal
    let lineage = db.get_full_lineage(5).unwrap();
    assert_eq!(lineage.len(), 5);
    
    // Test circular reference detection
    let circular_entry = TaxonomyEntry {
        id: 6, 
        name: "Circular", 
        parent_id: Some(6), // Self-reference
        rank: "species"
    };
    
    let result = db.insert_taxonomy_entries(&[circular_entry]);
    assert!(result.is_err(), "Should reject circular references");
}
```

#### Large Dataset Handling
```rust
#[test]
fn test_large_dataset_performance() {
    let db = create_test_database();
    
    // Insert large number of sequences (simulate real metagenome)
    let mut sequences = Vec::new();
    for i in 0..10000 {
        sequences.push(SequenceEntry {
            id: i,
            sequence_hash: format!("hash_{}", i),
            sequence_data: generate_random_sequence(150),
            length: 150,
            gc_content: 0.5,
            taxonomy_id: Some(1),
            source: "synthetic".to_string(),
            created_at: chrono::Utc::now(),
        });
    }
    
    let start_time = std::time::Instant::now();
    db.insert_sequences(&sequences).unwrap();
    let insert_time = start_time.elapsed();
    
    // Should complete within reasonable time
    assert!(insert_time.as_secs() < 30, "Large insert should complete in reasonable time");
    
    // Test query performance
    let start_time = std::time::Instant::now();
    let results = db.find_similar_sequences("ATCGATCGATCG", 8, 100).unwrap();
    let query_time = start_time.elapsed();
    
    assert!(query_time.as_millis() < 1000, "Query should be fast even with large dataset");
}
```

### 5. Quality Control Edge Cases

#### Paired-End Read Validation
```rust
#[test]
fn test_paired_end_validation() {
    let paired_reads = [
        ("ATCGATCGATCGATCG", "CGATCGATCGATCGAT"), // Normal pair
        ("ATCGATCGATCGATCG", "ATCGATCGATCGATCG"), // Identical (invalid)
        ("ATCGATCGATCGATCG", "NNNNNNNNNNNNNNNN"), // One read failed
        ("ATCGATCGATCGATCG", ""),                   // Missing mate
    ];
    
    for (read1, read2) in paired_reads {
        let result = validate_paired_reads(read1, read2, 150..500); // Expected insert size
        // Should detect various paired-end anomalies
    }
}
```

#### Contamination Detection
```rust
#[test]
fn test_contamination_detection() {
    let sequences = [
        "ATCGATCGATCGATCGATCGATCGATCG",           // Clean sequence
        "AGATCGGAAGAGCGTCGTGTAGGGAAA",            // Illumina adapter
        "CTGTCTCTTATACACATCTCCGAGCCC",            // Nextera adapter  
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAA",           // Poly-A tail
        "TTTTTTTTTTTTTTTTTTTTTTTTTTTT",           // Poly-T primer
    ];
    
    for seq in sequences {
        let result = detect_contamination(seq);
        // Should identify various contamination sources
    }
}
```

### 6. Integration Test Edge Cases

#### Memory Pressure Scenarios
```rust
#[test]
fn test_memory_pressure_handling() {
    // Test with limited memory settings
    let config = PipelineConfiguration {
        max_memory_gb: 1, // Very limited
        chunk_size: 100,
        ..Default::default()
    };
    
    // Process large dataset with memory constraints
    let large_dataset = generate_large_test_dataset(50000);
    let result = run_pipeline_with_config(&large_dataset, config);
    
    // Should complete without memory errors
    assert!(result.is_ok(), "Should handle memory pressure gracefully");
}
```

#### Multi-threaded Race Conditions
```rust
#[test]
fn test_concurrent_processing() {
    use std::sync::Arc;
    use std::thread;
    
    let db = Arc::new(create_test_database());
    let mut handles = Vec::new();
    
    // Spawn multiple threads doing database operations
    for i in 0..10 {
        let db_clone = Arc::clone(&db);
        let handle = thread::spawn(move || {
            let sequences = generate_test_sequences(100, i);
            db_clone.insert_sequences(&sequences).unwrap();
        });
        handles.push(handle);
    }
    
    // Wait for all threads to complete
    for handle in handles {
        handle.join().unwrap();
    }
    
    // Verify data integrity
    let stats = db.get_database_stats().unwrap();
    assert_eq!(stats.sequences_count, 1000); // All sequences should be inserted
}
```

## Implementation Priority

### Immediate (Critical for Scientific Validity)
1. Ambiguous base handling (N bases)
2. Extreme GC content validation  
3. Repetitive element assembly
4. Variable coverage handling

### High Priority (Important for Robustness)
1. Homopolymer detection
2. Graph topology extremes
3. Taxonomic hierarchy validation
4. Paired-end read validation

### Medium Priority (Enhanced Testing)
1. Contamination detection
2. Large dataset performance
3. Memory pressure handling
4. Multi-threaded safety

### Low Priority (Future Enhancement)
1. Advanced feature extraction edge cases
2. Complex graph algorithms
3. Specialized organism handling
4. Performance benchmarking

These additional tests will ensure the bioinformatics pipeline handles real-world biological data complexities while maintaining scientific accuracy and computational efficiency.