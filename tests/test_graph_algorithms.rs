//! Graph algorithm correctness tests
//! Tests assembly graph construction, repeat resolution, tip removal, bubble detection, and path finding

use anyhow::Result;
use meta_forge::core::data_structures::*;
use std::collections::HashSet;

#[cfg(test)]
mod graph_construction_tests {
    use super::*;

    fn create_test_fragment_with_nodes(node_sequences: Vec<&str>) -> Result<GraphFragment> {
        let mut fragment = GraphFragment::new(0);

        for (i, seq) in node_sequences.iter().enumerate() {
            let kmer = CanonicalKmer::new(seq)?;
            let mut node = GraphNode::new(kmer, seq.len());
            node.coverage = (i + 1) as u32; // Give different coverage values
            fragment.add_node(node);
        }

        Ok(fragment)
    }

    #[test]
    fn test_linear_path_construction() {
        let mut fragment = GraphFragment::new(0);

        // Create a linear path: A -> B -> C -> D
        let sequences = vec!["ATCG", "TCGA", "CGAT", "GATC"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // Add edges to create linear path
        for i in 0..nodes.len() - 1 {
            let edge = GraphEdge::new(nodes[i].hash, nodes[i + 1].hash, 3);
            fragment.add_edge(edge);
        }

        assert_eq!(fragment.nodes.len(), 4);
        assert_eq!(fragment.edges.len(), 3);

        // Check adjacency list
        let adj = fragment.get_adjacency_list();
        assert_eq!(adj.len(), 3); // Three nodes have outgoing edges

        // Each internal node should have degree 2 (except ends)
        for i in 0..nodes.len() - 1 {
            let neighbors = adj.get(&nodes[i].hash).unwrap();
            assert_eq!(neighbors.len(), 1);
            assert_eq!(neighbors[0], nodes[i + 1].hash);
        }
    }

    #[test]
    fn test_branching_structure_construction() {
        let mut fragment = GraphFragment::new(0);

        // Create branching structure: A -> B -> C
        //                                  \-> D
        let sequences = vec!["ATCG", "TCGA", "CGAT", "CGAA"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // A -> B
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        // B -> C
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[2].hash, 3));
        // B -> D (branch)
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[3].hash, 3));

        assert_eq!(fragment.nodes.len(), 4);
        assert_eq!(fragment.edges.len(), 3);

        // Check branching node (B) has two outgoing edges
        let adj = fragment.get_adjacency_list();
        let b_neighbors = adj.get(&nodes[1].hash).unwrap();
        assert_eq!(b_neighbors.len(), 2);
        assert!(b_neighbors.contains(&nodes[2].hash));
        assert!(b_neighbors.contains(&nodes[3].hash));
    }

    #[test]
    fn test_circular_structure_construction() {
        let mut fragment = GraphFragment::new(0);

        // Create circular structure: A -> B -> C -> A
        let sequences = vec!["ATCG", "TCGA", "CGAT"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // Create cycle
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[2].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[2].hash, nodes[0].hash, 3)); // Close the cycle

        assert_eq!(fragment.nodes.len(), 3);
        assert_eq!(fragment.edges.len(), 3);

        // In a cycle, all nodes should have both in-degree and out-degree of 1
        let adj = fragment.get_adjacency_list();
        for node in &nodes {
            let neighbors = adj.get(&node.hash).unwrap();
            assert_eq!(neighbors.len(), 1);
        }

        // Check that tips detection doesn't find any tips in a perfect cycle
        let tips = fragment.find_tips();
        assert!(tips.is_empty(), "Perfect cycle should have no tips");
    }

    #[test]
    fn test_complex_graph_with_multiple_components() {
        let mut fragment = GraphFragment::new(0);

        // Component 1: A -> B -> C
        let comp1_seqs = vec!["AAAA", "TTTT", "GGGG"];
        let mut comp1_nodes = Vec::new();

        for seq in &comp1_seqs {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            comp1_nodes.push(kmer);
            fragment.add_node(node);
        }

        fragment.add_edge(GraphEdge::new(comp1_nodes[0].hash, comp1_nodes[1].hash, 3));
        fragment.add_edge(GraphEdge::new(comp1_nodes[1].hash, comp1_nodes[2].hash, 3));

        // Component 2: D -> E (separate component)
        let comp2_seqs = vec!["CCCC", "ATGC"];
        let mut comp2_nodes = Vec::new();

        for seq in &comp2_seqs {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            comp2_nodes.push(kmer);
            fragment.add_node(node);
        }

        fragment.add_edge(GraphEdge::new(comp2_nodes[0].hash, comp2_nodes[1].hash, 3));

        assert_eq!(fragment.nodes.len(), 5);
        assert_eq!(fragment.edges.len(), 3);

        // Check that adjacency list reflects disconnected components
        let adj = fragment.get_adjacency_list();
        assert_eq!(adj.len(), 3); // Three nodes have outgoing edges

        // No edges should exist between components
        let comp1_neighbors: HashSet<u64> = comp1_nodes
            .iter()
            .filter_map(|n| adj.get(&n.hash))
            .flat_map(|neighbors| neighbors.iter())
            .cloned()
            .collect();

        let comp2_hashes: HashSet<u64> = comp2_nodes.iter().map(|n| n.hash).collect();
        assert!(comp1_neighbors.is_disjoint(&comp2_hashes));
    }

    #[test]
    fn test_self_loop_handling() {
        let mut fragment = GraphFragment::new(0);

        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let node = GraphNode::new(kmer.clone(), 4);
        fragment.add_node(node);

        // Create self-loop
        let self_edge = GraphEdge::new(kmer.hash, kmer.hash, 4);
        fragment.add_edge(self_edge);

        assert_eq!(fragment.nodes.len(), 1);
        assert_eq!(fragment.edges.len(), 1);

        let adj = fragment.get_adjacency_list();
        let neighbors = adj.get(&kmer.hash).unwrap();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(neighbors[0], kmer.hash); // Points to itself
    }
}

#[cfg(test)]
mod tip_detection_tests {
    use super::*;

    #[test]
    fn test_simple_tip_detection() {
        let mut fragment = GraphFragment::new(0);

        // Linear path with tips: A -> B -> C -> D
        let sequences = vec!["AAAA", "TTTT", "GGGG", "CCCC"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // Create linear chain
        for i in 0..nodes.len() - 1 {
            fragment.add_edge(GraphEdge::new(nodes[i].hash, nodes[i + 1].hash, 3));
        }

        let tips = fragment.find_tips();
        assert_eq!(tips.len(), 2); // Start and end are tips
        assert!(tips.contains(&nodes[0].hash)); // First node (no incoming)
        assert!(tips.contains(&nodes[3].hash)); // Last node (no outgoing)
    }

    #[test]
    fn test_branching_tips() {
        let mut fragment = GraphFragment::new(0);

        // Structure:   B -> D (tip)
        //             /
        //        A ->
        //             \
        //              C -> E (tip)

        let sequences = vec!["AAAA", "TTTT", "GGGG", "CCCC", "ATGC"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // A -> B, A -> C
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[2].hash, 3));
        // B -> D, C -> E
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[3].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[2].hash, nodes[4].hash, 3));

        let tips = fragment.find_tips();
        assert_eq!(tips.len(), 3); // A (no incoming), D, E (no outgoing)
        assert!(tips.contains(&nodes[0].hash)); // A
        assert!(tips.contains(&nodes[3].hash)); // D
        assert!(tips.contains(&nodes[4].hash)); // E
    }

    #[test]
    fn test_no_tips_in_cycle() {
        let mut fragment = GraphFragment::new(0);

        // Perfect cycle: A -> B -> C -> A
        let sequences = vec!["AAAA", "TTTT", "GGGG"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[2].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[2].hash, nodes[0].hash, 3));

        let tips = fragment.find_tips();
        assert!(tips.is_empty(), "Cycle should have no tips");
    }

    #[test]
    fn test_isolated_node_as_tip() {
        let mut fragment = GraphFragment::new(0);

        // Add isolated nodes (no edges)
        let sequences = vec!["AAAA", "TTTT", "GGGG"];

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            fragment.add_node(node);
        }

        // No edges added - all nodes are isolated
        let tips = fragment.find_tips();
        assert_eq!(tips.len(), 3); // All isolated nodes are tips
    }
}

#[cfg(test)]
mod bubble_detection_tests {
    use super::*;

    #[test]
    fn test_simple_bubble_detection() {
        let mut fragment = GraphFragment::new(0);

        // Simple bubble:     B
        //                   / \
        //              A ->     -> D
        //                   \ /
        //                    C

        let sequences = vec!["AAAA", "TTTT", "GGGG", "CCCC"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // A -> B, A -> C
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[2].hash, 3));
        // B -> D, C -> D
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[3].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[2].hash, nodes[3].hash, 3));

        let bubbles = fragment.find_bubbles();
        assert!(!bubbles.is_empty(), "Should detect bubble structure");

        if let Some(bubble) = bubbles.first() {
            assert_eq!(bubble.start_node, nodes[0].hash);
            assert_eq!(bubble.end_node, nodes[3].hash);
            assert_eq!(bubble.alternative_paths.len(), 2);
            assert_eq!(bubble.bubble_type, BubbleType::Simple);
        }
    }

    #[test]
    fn test_complex_bubble_with_multiple_paths() {
        let mut fragment = GraphFragment::new(0);

        // Complex bubble:    B -> E
        //                   /      \
        //              A ->  C -> F  -> H
        //                   \      /
        //                    D -> G

        let sequences = vec![
            "AAAA", "TTTT", "GGGG", "CCCC", "ATGC", "CGTA", "TACG", "GCTA",
        ];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // A -> B, C, D
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[2].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[3].hash, 3));

        // B -> E, C -> F, D -> G
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[4].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[2].hash, nodes[5].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[3].hash, nodes[6].hash, 3));

        // E, F, G -> H
        fragment.add_edge(GraphEdge::new(nodes[4].hash, nodes[7].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[5].hash, nodes[7].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[6].hash, nodes[7].hash, 3));

        let bubbles = fragment.find_bubbles();
        assert!(!bubbles.is_empty(), "Should detect complex bubble");

        if let Some(bubble) = bubbles.first() {
            assert_eq!(bubble.start_node, nodes[0].hash);
            assert_eq!(bubble.end_node, nodes[7].hash);
            assert_eq!(bubble.alternative_paths.len(), 3); // Three alternative paths
            assert_eq!(bubble.bubble_type, BubbleType::Complex);
        }
    }

    #[test]
    fn test_no_bubble_in_linear_path() {
        let mut fragment = GraphFragment::new(0);

        // Simple linear path: A -> B -> C -> D
        let sequences = vec!["AAAA", "TTTT", "GGGG", "CCCC"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        for i in 0..nodes.len() - 1 {
            fragment.add_edge(GraphEdge::new(nodes[i].hash, nodes[i + 1].hash, 3));
        }

        let bubbles = fragment.find_bubbles();
        assert!(bubbles.is_empty(), "Linear path should have no bubbles");
    }

    #[test]
    fn test_branching_without_reconvergence() {
        let mut fragment = GraphFragment::new(0);

        // Branching without reconvergence: A -> B -> C
        //                                     \-> D -> E

        let sequences = vec!["AAAA", "TTTT", "GGGG", "CCCC", "ATGC"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // A -> B
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        // B -> C, B -> D
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[2].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[3].hash, 3));
        // D -> E
        fragment.add_edge(GraphEdge::new(nodes[3].hash, nodes[4].hash, 3));

        let bubbles = fragment.find_bubbles();
        assert!(
            bubbles.is_empty(),
            "Branching without reconvergence should have no bubbles"
        );
    }
}

#[cfg(test)]
mod path_reconstruction_tests {
    use super::*;

    #[test]
    fn test_simple_path_reconstruction() {
        let mut fragment = GraphFragment::new(0);

        // Create path with known k-mers
        let sequences = vec!["ATCG", "TCGA", "CGAT"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        let path: Vec<u64> = nodes.iter().map(|n| n.hash).collect();
        let reconstructed = fragment.reconstruct_sequence_from_path(&path).unwrap();

        assert!(!reconstructed.is_empty());
        assert!(reconstructed.len() >= 4); // At least one k-mer length
    }

    #[test]
    fn test_empty_path_reconstruction() {
        let fragment = GraphFragment::new(0);
        let empty_path = vec![];

        let result = fragment
            .reconstruct_sequence_from_path(&empty_path)
            .unwrap();
        assert!(result.is_empty());
    }

    #[test]
    fn test_single_node_path() {
        let mut fragment = GraphFragment::new(0);

        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let node = GraphNode::new(kmer.clone(), 4);
        fragment.add_node(node);

        let path = vec![kmer.hash];
        let reconstructed = fragment.reconstruct_sequence_from_path(&path).unwrap();

        assert_eq!(reconstructed, "ATCG");
    }

    #[test]
    fn test_path_with_missing_nodes() {
        let fragment = GraphFragment::new(0);

        // Path with non-existent node hashes
        let path = vec![12345, 67890];
        let reconstructed = fragment.reconstruct_sequence_from_path(&path).unwrap();

        // Should handle missing nodes gracefully
        assert!(!reconstructed.is_empty());
    }

    #[test]
    fn test_coverage_calculation_for_path() {
        let mut fragment = GraphFragment::new(0);

        let sequences = ["ATCG", "TCGA"];
        let mut nodes = Vec::new();

        for (i, seq) in sequences.iter().enumerate() {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let mut node = GraphNode::new(kmer.clone(), 4);
            node.coverage = (i + 1) as u32 * 10; // Coverage 10, 20
            nodes.push(kmer);
            fragment.add_node(node);
        }

        let path: Vec<u64> = nodes.iter().map(|n| n.hash).collect();
        let avg_coverage = fragment.calculate_path_coverage_from_hashes(&path);

        assert_eq!(avg_coverage, 15.0); // (10 + 20) / 2
    }

    #[test]
    fn test_coverage_calculation_empty_path() {
        let fragment = GraphFragment::new(0);
        let empty_path = vec![];

        let coverage = fragment.calculate_path_coverage_from_hashes(&empty_path);
        assert_eq!(coverage, 0.0);
    }
}

#[cfg(test)]
mod graph_statistics_tests {
    use super::*;

    #[test]
    fn test_adjacency_list_generation() {
        let mut fragment = GraphFragment::new(0);

        let sequences = vec!["AAAA", "TTTT", "GGGG"];
        let mut nodes = Vec::new();

        for seq in &sequences {
            let kmer = CanonicalKmer::new(seq).unwrap();
            let node = GraphNode::new(kmer.clone(), 4);
            nodes.push(kmer);
            fragment.add_node(node);
        }

        // A -> B -> C
        fragment.add_edge(GraphEdge::new(nodes[0].hash, nodes[1].hash, 3));
        fragment.add_edge(GraphEdge::new(nodes[1].hash, nodes[2].hash, 3));

        let adj = fragment.get_adjacency_list();

        assert_eq!(adj.len(), 2); // Two nodes have outgoing edges
        assert!(adj.contains_key(&nodes[0].hash));
        assert!(adj.contains_key(&nodes[1].hash));
        assert!(!adj.contains_key(&nodes[2].hash)); // Leaf node has no outgoing edges

        let a_neighbors = adj.get(&nodes[0].hash).unwrap();
        assert_eq!(a_neighbors.len(), 1);
        assert_eq!(a_neighbors[0], nodes[1].hash);

        let b_neighbors = adj.get(&nodes[1].hash).unwrap();
        assert_eq!(b_neighbors.len(), 1);
        assert_eq!(b_neighbors[0], nodes[2].hash);
    }

    #[test]
    fn test_coverage_statistics_update() {
        let mut fragment = GraphFragment::new(0);

        // Add nodes with different coverages
        let coverages = [5, 10, 15, 20];
        for (i, &cov) in coverages.iter().enumerate() {
            let seq = format!("A{i:03}"); // A000, A001, etc.
            let kmer = CanonicalKmer::new(&seq).unwrap();
            let mut node = GraphNode::new(kmer, 4);
            node.coverage = cov;
            fragment.add_node(node);
        }

        // Coverage stats should be automatically updated
        assert_eq!(fragment.coverage_stats.total_nodes, 4);
        assert_eq!(fragment.coverage_stats.mean_coverage, 12.5); // (5+10+15+20)/4
        assert_eq!(fragment.coverage_stats.median_coverage, 15); // Median of sorted [5,10,15,20]

        // Check distribution
        assert_eq!(
            fragment.coverage_stats.coverage_distribution.get(&5),
            Some(&1)
        );
        assert_eq!(
            fragment.coverage_stats.coverage_distribution.get(&10),
            Some(&1)
        );
        assert_eq!(
            fragment.coverage_stats.coverage_distribution.get(&15),
            Some(&1)
        );
        assert_eq!(
            fragment.coverage_stats.coverage_distribution.get(&20),
            Some(&1)
        );
    }

    #[test]
    fn test_graph_merging_statistics() {
        let mut fragment1 = GraphFragment::new(0);
        let mut fragment2 = GraphFragment::new(1);

        // Add nodes to each fragment
        let kmer1 = CanonicalKmer::new("AAAA").unwrap();
        let kmer2 = CanonicalKmer::new("TTTT").unwrap();
        let kmer_common = CanonicalKmer::new("GGGG").unwrap();

        let mut node1 = GraphNode::new(kmer1, 4);
        node1.coverage = 10;
        fragment1.add_node(node1);

        let mut node_common1 = GraphNode::new(kmer_common.clone(), 4);
        node_common1.coverage = 5;
        fragment1.add_node(node_common1);

        let mut node2 = GraphNode::new(kmer2, 4);
        node2.coverage = 15;
        fragment2.add_node(node2);

        let mut node_common2 = GraphNode::new(kmer_common.clone(), 4);
        node_common2.coverage = 7;
        fragment2.add_node(node_common2);

        let original_stats1 = fragment1.coverage_stats.clone();

        fragment1.merge_with(fragment2).unwrap();

        // Check merged statistics
        assert_eq!(fragment1.coverage_stats.total_nodes, 3); // 3 unique nodes
        assert!(fragment1.coverage_stats.mean_coverage > original_stats1.mean_coverage);

        // Common node should have combined coverage
        let merged_common = fragment1.nodes.get(&kmer_common.hash).unwrap();
        assert_eq!(merged_common.coverage, 12); // 5 + 7
    }
}
