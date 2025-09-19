//! Basic Graph Construction Utilities
//! ==================================
//!
//! Simple utilities for DNA sequence graph construction optimized for laptop use.
//! Focused on clarity and efficiency over complex parallelization.

use ahash::AHashMap;
use anyhow::Result;

/// Calculate sequence complexity using Shannon entropy
pub fn calculate_sequence_complexity(sequence: &str) -> f64 {
    let mut counts = [0usize; 4]; // A, C, G, T
    let mut total = 0;

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
            _ => {} // Ignore invalid characters
        }
    }

    if total == 0 {
        return 0.0;
    }

    let mut entropy = 0.0;
    for count in counts {
        if count > 0 {
            let p = count as f64 / total as f64;
            entropy -= p * p.log2();
        }
    }

    entropy
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

/// Simple k-mer overlap calculation
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

/// Simple edge validation for basic connectivity
pub fn validate_edge_connectivity(edges: &[(String, String)]) -> Result<()> {
    let mut node_degrees: AHashMap<String, (usize, usize)> = AHashMap::new();

    // Count in-degree and out-degree for each node
    for (from, to) in edges {
        let from_entry = node_degrees.entry(from.clone()).or_insert((0, 0));
        from_entry.1 += 1; // out-degree

        let to_entry = node_degrees.entry(to.clone()).or_insert((0, 0));
        to_entry.0 += 1; // in-degree
    }

    // Check for basic connectivity patterns
    let mut start_nodes = 0;
    let mut end_nodes = 0;

    for (in_deg, out_deg) in node_degrees.values() {
        if *in_deg == 0 && *out_deg > 0 {
            start_nodes += 1;
        } else if *out_deg == 0 && *in_deg > 0 {
            end_nodes += 1;
        }
    }

    if start_nodes == 0 || end_nodes == 0 {
        return Err(anyhow::anyhow!("Graph lacks proper start or end nodes"));
    }

    Ok(())
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
        assert_eq!(calculate_kmer_overlap("ATCG", "GAAT"), 0);
    }

    #[test]
    fn test_edge_validation() {
        let valid_edges = vec![
            ("A".to_string(), "B".to_string()),
            ("B".to_string(), "C".to_string()),
        ];
        assert!(validate_edge_connectivity(&valid_edges).is_ok());

        let invalid_edges = vec![
            ("A".to_string(), "B".to_string()),
            ("C".to_string(), "B".to_string()),
        ];
        // This should fail because there's no clear path
        assert!(validate_edge_connectivity(&invalid_edges).is_err());
    }
}