//! Enhanced Graph Construction Utilities
//! ====================================
//!
//! Advanced algorithms for high-accuracy DNA sequence graph construction optimized for laptop use.
//! Implements sophisticated assembly quality metrics and adaptive algorithms for improved accuracy.

use ahash::{AHashMap, AHashSet};
use anyhow::{Result, anyhow};
use std::time::Instant;
use rayon::prelude::*;

/// Calculate sequence complexity using Shannon entropy with enhanced accuracy
pub fn calculate_sequence_complexity(sequence: &str) -> f64 {
    let mut counts = [0usize; 4]; // A, C, G, T
    let mut total = 0;
    let mut ambiguous_count = 0;

    for c in sequence.chars() {
        match c.to_ascii_uppercase() {
            'A' => {
                counts[0] += 1;
                total += 1;
            }
            'C' => {
                counts[1] += 1;
                total += 1;
            }
            'G' => {
                counts[2] += 1;
                total += 1;
            }
            'T' => {
                counts[3] += 1;
                total += 1;
            }
            'N' | 'X' => ambiguous_count += 1, // Track ambiguous bases
            _ => {} // Ignore other invalid characters
        }
    }

    if total == 0 {
        return 0.0;
    }

    // Calculate Shannon entropy
    let mut entropy = 0.0;
    for count in counts {
        if count > 0 {
            let p = count as f64 / total as f64;
            entropy -= p * p.log2();
        }
    }

    // Penalize sequences with many ambiguous bases
    let ambiguous_penalty = (ambiguous_count as f64) / (total + ambiguous_count) as f64;
    entropy * (1.0 - ambiguous_penalty * 0.5)
}

/// Calculate GC content percentage
pub fn calculate_gc_content(sequence: &str) -> f64 {
    let mut gc_count = 0;
    let mut total_count = 0;

    for c in sequence.chars() {
        match c.to_ascii_uppercase() {
            'G' | 'C' => {
                gc_count += 1;
                total_count += 1;
            }
            'A' | 'T' => {
                total_count += 1;
            }
            _ => {} // Ignore invalid characters
        }
    }

    if total_count == 0 {
        0.0
    } else {
        (gc_count as f64 / total_count as f64) * 100.0
    }
}

/// Enhanced k-mer overlap calculation with quality scoring
pub fn calculate_kmer_overlap(kmer1: &str, kmer2: &str) -> usize {
    let k1_len = kmer1.len();
    let k2_len = kmer2.len();

    // Find the maximum overlap between suffix of kmer1 and prefix of kmer2
    let max_overlap = k1_len.min(k2_len);

    for overlap_len in (1..=max_overlap).rev() {
        let suffix = &kmer1[k1_len - overlap_len..];
        let prefix = &kmer2[..overlap_len];

        if suffix == prefix {
            return overlap_len;
        }
    }

    0
}

/// Calculate k-mer overlap quality based on complexity and GC content
pub fn calculate_overlap_quality(kmer1: &str, kmer2: &str, overlap_len: usize) -> f64 {
    if overlap_len == 0 {
        return 0.0;
    }

    let k1_len = kmer1.len();
    let overlap_region = &kmer1[k1_len - overlap_len..];

    // Quality factors
    let complexity = calculate_sequence_complexity(overlap_region);
    let gc_content = calculate_gc_content(overlap_region);

    // Ideal GC content is around 50%, complexity should be high
    let gc_score = 1.0 - ((gc_content - 50.0) / 50.0).abs();
    let complexity_score = (complexity / 2.0).min(1.0); // Normalize to [0,1]

    // Length bonus for longer overlaps
    let length_bonus = (overlap_len as f64 / kmer1.len() as f64).min(1.0);

    (gc_score + complexity_score + length_bonus) / 3.0
}

/// Enhanced edge validation with comprehensive graph analysis
pub fn validate_edge_connectivity(edges: &[(String, String)]) -> Result<GraphValidationResult> {
    if edges.is_empty() {
        return Ok(GraphValidationResult {
            is_valid: true,
            start_nodes: 0,
            end_nodes: 0,
            total_nodes: 0,
            connected_components: 0,
            issues: vec!["No edges found - graph consists of isolated nodes".to_string()],
        });
    }

    let mut node_degrees: AHashMap<String, (usize, usize)> = AHashMap::new();
    let mut all_nodes: AHashSet<String> = AHashSet::new();

    // Count in-degree and out-degree for each node
    for (from, to) in edges {
        all_nodes.insert(from.clone());
        all_nodes.insert(to.clone());

        let from_entry = node_degrees.entry(from.clone()).or_insert((0, 0));
        from_entry.1 += 1; // out-degree

        let to_entry = node_degrees.entry(to.clone()).or_insert((0, 0));
        to_entry.0 += 1; // in-degree
    }

    // Analyze connectivity patterns
    let mut start_nodes = 0;
    let mut end_nodes = 0;
    let mut branching_nodes = 0;
    let mut linear_nodes = 0;
    let mut issues = Vec::new();

    for (node, (in_deg, out_deg)) in &node_degrees {
        match (in_deg, out_deg) {
            (0, out_deg) if *out_deg > 0 => start_nodes += 1,
            (in_deg, 0) if *in_deg > 0 => end_nodes += 1,
            (1, 1) => linear_nodes += 1,
            (in_deg, out_deg) if *in_deg > 1 || *out_deg > 1 => {
                branching_nodes += 1;
                if *in_deg > 2 || *out_deg > 2 {
                    issues.push(format!("High degree node '{}': in={}, out={}", node, in_deg, out_deg));
                }
            },
            _ => {},
        }
    }

    // Estimate connected components (simplified)
    let connected_components = estimate_connected_components(edges);

    // Validate graph structure
    let is_valid = start_nodes > 0 && end_nodes > 0;

    if start_nodes == 0 {
        issues.push("No start nodes found (nodes with in-degree 0)".to_string());
    }
    if end_nodes == 0 {
        issues.push("No end nodes found (nodes with out-degree 0)".to_string());
    }
    if connected_components > 1 {
        issues.push(format!("Graph has {} disconnected components", connected_components));
    }

    Ok(GraphValidationResult {
        is_valid,
        start_nodes,
        end_nodes,
        total_nodes: all_nodes.len(),
        connected_components,
        issues,
    })
}

/// Result of graph validation analysis
#[derive(Debug, Clone)]
pub struct GraphValidationResult {
    pub is_valid: bool,
    pub start_nodes: usize,
    pub end_nodes: usize,
    pub total_nodes: usize,
    pub connected_components: usize,
    pub issues: Vec<String>,
}

/// Estimate number of connected components in the graph
fn estimate_connected_components(edges: &[(String, String)]) -> usize {
    let mut adjacency: AHashMap<String, Vec<String>> = AHashMap::new();
    let mut all_nodes: AHashSet<String> = AHashSet::new();

    // Build adjacency list (undirected for connectivity analysis)
    for (from, to) in edges {
        all_nodes.insert(from.clone());
        all_nodes.insert(to.clone());

        adjacency.entry(from.clone()).or_default().push(to.clone());
        adjacency.entry(to.clone()).or_default().push(from.clone());
    }

    let mut visited: AHashSet<String> = AHashSet::new();
    let mut components = 0;

    for node in all_nodes {
        if !visited.contains(&node) {
            // DFS to mark all nodes in this component
            let mut stack = vec![node.clone()];
            while let Some(current) = stack.pop() {
                if visited.insert(current.clone()) {
                    if let Some(neighbors) = adjacency.get(&current) {
                        for neighbor in neighbors {
                            if !visited.contains(neighbor) {
                                stack.push(neighbor.clone());
                            }
                        }
                    }
                }
            }
            components += 1;
        }
    }

    components
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sequence_complexity() {
        // Uniform sequence should have maximum entropy
        assert!(calculate_sequence_complexity("ATCGATCGATCG") > 1.5);

        // Homogeneous sequence should have zero entropy
        assert!(calculate_sequence_complexity("AAAAAAAAAA") < 0.1);
    }

    #[test]
    fn test_gc_content() {
        assert!((calculate_gc_content("ATCG") - 50.0).abs() < 0.1);
        assert!((calculate_gc_content("AAAA") - 0.0).abs() < 0.1);
        assert!((calculate_gc_content("GGCC") - 100.0).abs() < 0.1);
    }

    #[test]
    fn test_kmer_overlap() {
        assert_eq!(calculate_kmer_overlap("ATCG", "TCGA"), 3);
        assert_eq!(calculate_kmer_overlap("ATCG", "CGAA"), 2);
        assert_eq!(calculate_kmer_overlap("ATCG", "GAAT"), 1); // "G" overlaps
        assert_eq!(calculate_kmer_overlap("ATCG", "TTTX"), 0); // No overlap
    }

    #[test]
    fn test_edge_validation() {
        // Valid linear path: A -> B -> C
        let valid_edges = vec![
            ("A".to_string(), "B".to_string()),
            ("B".to_string(), "C".to_string()),
        ];
        let result = validate_edge_connectivity(&valid_edges).unwrap();
        assert!(result.is_valid);
        assert_eq!(result.connected_components, 1);

        // Valid branching: A -> B, C -> B (has start A,C and end B)
        let valid_branching = vec![
            ("A".to_string(), "B".to_string()),
            ("C".to_string(), "B".to_string()),
        ];
        let result = validate_edge_connectivity(&valid_branching).unwrap();
        assert!(result.is_valid);
        assert_eq!(result.start_nodes, 2); // A and C
        assert_eq!(result.end_nodes, 1);   // B

        // Invalid: no end nodes (circular reference)
        let invalid_circular = vec![
            ("A".to_string(), "B".to_string()),
            ("B".to_string(), "A".to_string()),
        ];
        let result = validate_edge_connectivity(&invalid_circular).unwrap();
        assert!(!result.is_valid);
        assert_eq!(result.start_nodes, 0);
        assert_eq!(result.end_nodes, 0);

        // Invalid: no start nodes (all nodes have incoming edges)
        let invalid_no_start = vec![
            ("A".to_string(), "B".to_string()),
            ("C".to_string(), "A".to_string()),
            ("B".to_string(), "C".to_string()),
        ];
        let result = validate_edge_connectivity(&invalid_no_start).unwrap();
        assert!(!result.is_valid);
        assert_eq!(result.start_nodes, 0);
    }

    #[test]
    fn test_laptop_assembly_scenarios() {
        // Test realistic k-mer overlaps for laptop assembly
        // Typical k-mer sizes: 21, 31, 41
        let kmer_21_1 = "ATCGATCGATCGATCGATCGA"; // 21-mer
        let kmer_21_2 = "TCGATCGATCGATCGATCGAT"; // 20-overlap with first
        assert_eq!(calculate_kmer_overlap(kmer_21_1, kmer_21_2), 20);

        // Test edge validation for typical contig assembly graph
        // This represents a simplified assembly graph with proper ratios
        let assembly_edges = vec![
            ("contig_1".to_string(), "contig_2".to_string()),
            ("contig_2".to_string(), "contig_3".to_string()),
            ("contig_4".to_string(), "contig_2".to_string()), // Branch merging
        ];
        let result = validate_edge_connectivity(&assembly_edges).unwrap();
        assert!(result.is_valid);
        assert_eq!(result.start_nodes, 2); // contig_1 and contig_4
        assert_eq!(result.end_nodes, 1);   // contig_3

        // Test sequence complexity for assembly quality
        let high_complexity = "ATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        let low_complexity = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

        assert!(calculate_sequence_complexity(high_complexity) > 1.0);
        assert!(calculate_sequence_complexity(low_complexity) < 0.1);

        // Test GC content for bacterial genomes (typical range 30-70%)
        let typical_bacterial = "ATCGATCGCCGGATCGATCGATCG"; // Mixed content
        let gc_content = calculate_gc_content(typical_bacterial);
        assert!(gc_content >= 30.0 && gc_content <= 70.0);
    }
}