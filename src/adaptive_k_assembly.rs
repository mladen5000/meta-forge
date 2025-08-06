use std::collections::{HashMap, HashSet, VecDeque};
use anyhow::Result;

/// Variable-order de Bruijn graph that can adapt k-mer size per edge
#[derive(Default)]
pub struct AdaptiveGraph {
    /// Base k-mer size (minimum)
    base_k: usize,
    /// Max k-mer size allowed
    max_k: usize,
    /// Node data: minimiser hash -> (sequence context, coverage)
    nodes: HashMap<u64, NodeData>,
    /// Edges with variable k: (from, to) -> EdgeData
    edges: HashMap<(u64, u64), EdgeData>,
    /// Coordinate mapping: contig_pos -> (read_id, read_pos)
    coord_map: HashMap<ContigCoord, Vec<ReadCoord>>,
}

#[derive(Clone)]
struct NodeData {
    /// Actual k-mer sequence (for reconstruction)
    sequence: String,
    /// Coverage depth
    coverage: u32,
    /// Local complexity score (higher = more repetitive)
    complexity: f64,
}

#[derive(Clone)]
struct EdgeData {
    /// Effective k-mer size for this edge
    k_size: usize,
    /// Support count
    weight: u32,
    /// Overlap sequence between nodes
    overlap: String,
}

#[derive(Hash, PartialEq, Eq, Clone)]
struct ContigCoord {
    contig_id: usize,
    position: usize,
}

#[derive(Clone)]
struct ReadCoord {
    read_id: usize,
    read_pos: usize,
}

impl AdaptiveGraph {
    pub fn new(base_k: usize, max_k: usize) -> Self {
        Self {
            base_k,
            max_k,
            nodes: HashMap::new(),
            edges: HashMap::new(),
            coord_map: HashMap::new(),
        }
    }

    /// Add a read path, dynamically adjusting k in repetitive regions
    pub fn add_read_with_adaptive_k(&mut self, read_id: usize, sequence: &str) -> Result<()> {
        let minimisers = self.extract_minimisers_with_positions(sequence)?;
        
        // Calculate local complexity for each position
        let complexity_scores = self.calculate_local_complexity(sequence);
        
        // Build path with variable k
        for window in minimisers.windows(2) {
            let (pos1, min1) = window[0];
            let (pos2, min2) = window[1];
            
            // Determine k-size based on local complexity
            let avg_complexity = (complexity_scores[pos1] + complexity_scores[pos2]) / 2.0;
            let k_size = self.adaptive_k_size(avg_complexity);
            
            // Extract actual k-mer sequences at chosen k-size
            let kmer1 = self.extract_kmer_at_k(sequence, pos1, k_size)?;
            let kmer2 = self.extract_kmer_at_k(sequence, pos2, k_size)?;
            
            // Update nodes
            self.update_node(min1, &kmer1, avg_complexity);
            self.update_node(min2, &kmer2, avg_complexity);
            
            // Update edge with overlap information
            let overlap = self.calculate_overlap(&kmer1, &kmer2, k_size);
            self.update_edge((min1, min2), k_size, overlap);
            
            // Store coordinate mapping
            self.add_coordinate_mapping(read_id, pos1, pos2);
        }
        
        Ok(())
    }

    /// Calculate local sequence complexity using entropy
    fn calculate_local_complexity(&self, sequence: &str) -> Vec<f64> {
        let window_size = 50;
        let seq_bytes = sequence.as_bytes();
        let mut complexity = vec![0.0; sequence.len()];
        
        for i in 0..sequence.len() {
            let start = i.saturating_sub(window_size / 2);
            let end = (i + window_size / 2).min(sequence.len());
            let window = &seq_bytes[start..end];
            
            // Calculate Shannon entropy
            let mut counts = [0u32; 4]; // A, C, G, T
            for &byte in window {
                match byte {
                    b'A' | b'a' => counts[0] += 1,
                    b'C' | b'c' => counts[1] += 1,
                    b'G' | b'g' => counts[2] += 1,
                    b'T' | b't' => counts[3] += 1,
                    _ => {} // Skip ambiguous bases
                }
            }
            
            let total = counts.iter().sum::<u32>() as f64;
            if total > 0.0 {
                let entropy = counts
                    .iter()
                    .filter(|&&c| c > 0)
                    .map(|&c| {
                        let p = c as f64 / total;
                        -p * p.log2()
                    })
                    .sum::<f64>();
                
                // Low entropy = high repetitiveness = high complexity score
                complexity[i] = 2.0 - entropy; // Max entropy is 2.0 for DNA
            }
        }
        
        complexity
    }

    /// Determine k-size based on complexity score
    fn adaptive_k_size(&self, complexity: f64) -> usize {
        // Higher complexity (more repetitive) -> larger k
        let normalized = (complexity / 2.0).clamp(0.0, 1.0);
        let k_range = self.max_k - self.base_k;
        self.base_k + (normalized * k_range as f64) as usize
    }

    fn extract_kmer_at_k(&self, sequence: &str, pos: usize, k: usize) -> Result<String> {
        if pos + k <= sequence.len() {
            Ok(sequence[pos..pos + k].to_string())
        } else {
            // Handle edge case - use available sequence
            Ok(sequence[pos..].to_string())
        }
    }

    fn calculate_overlap(&self, kmer1: &str, kmer2: &str, k: usize) -> String {
        // Find maximum overlap between k-mers
        let max_overlap = (k - 1).min(kmer1.len()).min(kmer2.len());
        
        for i in (1..=max_overlap).rev() {
            if kmer1.ends_with(&kmer2[..i]) {
                return kmer2[i..].to_string();
            }
        }
        
        // No overlap found - full second k-mer
        kmer2.to_string()
    }

    fn update_node(&mut self, hash: u64, sequence: &str, complexity: f64) {
        self.nodes
            .entry(hash)
            .and_modify(|node| {
                node.coverage += 1;
                node.complexity = (node.complexity + complexity) / 2.0; // Running average
            })
            .or_insert(NodeData {
                sequence: sequence.to_string(),
                coverage: 1,
                complexity,
            });
    }

    fn update_edge(&mut self, edge_key: (u64, u64), k_size: usize, overlap: String) {
        self.edges
            .entry(edge_key)
            .and_modify(|edge| {
                edge.weight += 1;
                // Keep the k_size that has most support
                if edge.weight == 1 {
                    edge.k_size = k_size;
                    edge.overlap = overlap.clone();
                }
            })
            .or_insert(EdgeData {
                k_size,
                weight: 1,
                overlap,
            });
    }

    fn add_coordinate_mapping(&mut self, read_id: usize, pos1: usize, pos2: usize) {
        // This is simplified - in practice you'd track more detailed mappings
        // for each assembled position back to source reads
    }

    /// Assembly with variable-k consideration
    pub fn assemble_adaptive_contigs(&self) -> Vec<AdaptiveContig> {
        let mut visited = HashSet::new();
        let mut contigs = Vec::new();
        
        // Find high-coverage starting nodes
        let mut start_candidates: Vec<_> = self.nodes
            .iter()
            .filter(|(_, data)| data.coverage >= 3) // Minimum coverage threshold
            .map(|(&hash, _)| hash)
            .collect();
        
        start_candidates.sort_by_key(|&hash| std::cmp::Reverse(self.nodes[&hash].coverage));
        
        for start_node in start_candidates {
            if visited.contains(&start_node) {
                continue;
            }
            
            if let Some(contig) = self.extend_contig_adaptive(start_node, &mut visited) {
                contigs.push(contig);
            }
        }
        
        contigs
    }

    fn extend_contig_adaptive(&self, start: u64, visited: &mut HashSet<u64>) -> Option<AdaptiveContig> {
        let mut path = vec![start];
        let mut current = start;
        let mut total_sequence = self.nodes[&start].sequence.clone();
        
        // Extend forward
        while let Some(next) = self.find_best_next_node(current, visited) {
            if visited.contains(&next) {
                break;
            }
            
            path.push(next);
            
            // Get edge data for overlap calculation
            if let Some(edge_data) = self.edges.get(&(current, next)) {
                total_sequence.push_str(&edge_data.overlap);
            } else {
                // Fallback - just append the k-mer
                total_sequence.push_str(&self.nodes[&next].sequence);
            }
            
            current = next;
        }
        
        // Mark all nodes in path as visited
        for &node in &path {
            visited.insert(node);
        }
        
        if path.len() >= 2 {
            Some(AdaptiveContig {
                nodes: path,
                sequence: total_sequence,
                avg_coverage: self.calculate_path_coverage(&path),
            })
        } else {
            None
        }
    }

    fn find_best_next_node(&self, current: u64, visited: &HashSet<u64>) -> Option<u64> {
        self.edges
            .iter()
            .filter(|((from, to), _)| *from == current && !visited.contains(to))
            .max_by_key(|(_, edge_data)| edge_data.weight)
            .map(|((_, to), _)| *to)
    }

    fn calculate_path_coverage(&self, path: &[u64]) -> f64 {
        path.iter()
            .map(|&node| self.nodes[&node].coverage as f64)
            .sum::<f64>() / path.len() as f64
    }

    fn extract_minimisers_with_positions(&self, sequence: &str) -> Result<Vec<(usize, u64)>> {
        // Reuse your existing minimiser extraction logic
        // This should return (position, minimiser_hash) pairs
        Ok(vec![]) // Placeholder
    }
}

#[derive(Clone)]
pub struct AdaptiveContig {
    pub nodes: Vec<u64>,
    pub sequence: String,
    pub avg_coverage: f64,
}

/// Back-port contig coordinates to original read positions
pub fn map_contig_to_reads(
    contig: &AdaptiveContig,
    graph: &AdaptiveGraph,
) -> HashMap<usize, Vec<ReadCoord>> {
    let mut mapping = HashMap::new();
    
    for (contig_pos, &node_hash) in contig.nodes.iter().enumerate() {
        // Look up all read coordinates that contributed to this node
        let contig_coord = ContigCoord {
            contig_id: 0, // Would need to track contig IDs properly
            position: contig_pos,
        };
        
        if let Some(read_coords) = graph.coord_map.get(&contig_coord) {
            mapping.insert(contig_pos, read_coords.clone());
        }
    }
    
    mapping
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_adaptive_k_selection() {
        let graph = AdaptiveGraph::new(15, 31);
        
        // Low complexity (uniform) -> smaller k
        assert_eq!(graph.adaptive_k_size(0.5), 17);
        
        // High complexity (repetitive) -> larger k  
        assert_eq!(graph.adaptive_k_size(1.5), 27);
    }
    
    #[test]
    fn test_complexity_calculation() {
        let graph = AdaptiveGraph::new(15, 31);
        let uniform_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let complex_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTAC";
        
        let uniform_complexity = graph.calculate_local_complexity(uniform_seq);
        let complex_complexity = graph.calculate_local_complexity(complex_seq);
        
        // Uniform sequence should have higher complexity score (lower entropy)
        assert!(uniform_complexity[15] > complex_complexity[15]);
    }
}