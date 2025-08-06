//! Variable-order (adaptive‐k) de Bruijn graph for metagenomic assembly.
//!
//! * `base_k` is the lowest k-mer length allowed.
//! * `max_k` is the highest k-mer length allowed.
//!
//! Each **edge** can negotiate its own effective `k_size` depending on local
//! sequence complexity: low-entropy regions (repeats) get bigger k to disambiguate,
//! high-entropy regions keep smaller k for sensitivity.

use anyhow::{Result, bail};
use std::collections::{HashMap, HashSet};

/// Graph container.
#[derive(Default, Debug)]
pub struct AdaptiveGraph {
    /// Minimum k-mer length.
    base_k: usize,
    /// Maximum k-mer length.
    max_k: usize,

    /// Map `minimiser_hash → NodeData`.
    nodes: HashMap<u64, NodeData>,

    /// Map `(from_hash, to_hash) → EdgeData`.
    edges: HashMap<(u64, u64), EdgeData>,

    /// Map contig coordinates back to the reads that produced them.
    coord_map: HashMap<ContigCoord, Vec<ReadCoord>>,
}

/// Metadata stored per node (k-mer minimiser).
#[derive(Clone, Debug)]
struct NodeData {
    sequence: String,
    coverage: u32,
    complexity: f64,
}

/// Metadata stored per edge in the variable-k graph.
#[derive(Clone, Debug)]
struct EdgeData {
    k_size: usize,
    weight: u32,
    overlap: String,
}

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
struct ContigCoord {
    contig_id: usize,
    position: usize,
}

#[derive(Clone, Debug)]
struct ReadCoord {
    read_id: usize,
    read_pos: usize,
}

impl AdaptiveGraph {
    // --------------------------------------------------------------------- //
    //  Construction
    // --------------------------------------------------------------------- //

    pub fn new(base_k: usize, max_k: usize) -> Self {
        assert!(base_k > 0 && max_k >= base_k);
        Self {
            base_k,
            max_k,
            ..Default::default()
        }
    }

    /// Add a read and adapt k per window depending on local complexity.
    pub fn add_read_with_adaptive_k(&mut self, read_id: usize, sequence: &str) -> Result<()> {
        if sequence.len() < self.base_k {
            bail!("Sequence shorter than base_k");
        }
        let minimisers = self.extract_minimisers_with_positions(sequence)?;

        let complexity_scores = self.calculate_local_complexity(sequence);

        for window in minimisers.windows(2) {
            let (pos1, min1) = window[0];
            let (pos2, min2) = window[1];
            let avg_complexity = (complexity_scores[pos1] + complexity_scores[pos2]) / 2.0;
            let k_size = self.adaptive_k_size(avg_complexity);

            let kmer1 = self.extract_kmer_at_k(sequence, pos1, k_size)?;
            let kmer2 = self.extract_kmer_at_k(sequence, pos2, k_size)?;

            self.update_node(min1, &kmer1, avg_complexity);
            self.update_node(min2, &kmer2, avg_complexity);

            let overlap = self.calculate_overlap(&kmer1, &kmer2, k_size);
            self.update_edge((min1, min2), k_size, overlap);

            self.add_coordinate_mapping(read_id, pos1, pos2);
        }
        Ok(())
    }

    // --------------------------------------------------------------------- //
    //  Assembly
    // --------------------------------------------------------------------- //

    /// Greedy walk to produce contigs, respecting variable-k overlaps.
    pub fn assemble_adaptive_contigs(&self) -> Vec<AdaptiveContig> {
        let mut visited = HashSet::new();
        let mut contigs = Vec::new();

        // Start from high-coverage nodes.
        let mut start_candidates: Vec<_> = self
            .nodes
            .iter()
            .filter(|(_, d)| d.coverage >= 3)
            .map(|(&h, _)| h)
            .collect();
        start_candidates.sort_by_key(|&h| std::cmp::Reverse(self.nodes[&h].coverage));

        for start in start_candidates {
            if visited.contains(&start) {
                continue;
            }
            if let Some(contig) = self.extend_contig_adaptive(start, &mut visited) {
                contigs.push(contig);
            }
        }
        contigs
    }

    // --------------------------------------------------------------------- //
    //  Internal helpers
    // --------------------------------------------------------------------- //

    fn extend_contig_adaptive(
        &self,
        start: u64,
        visited: &mut HashSet<u64>,
    ) -> Option<AdaptiveContig> {
        let mut path = vec![start];
        let mut seq = self.nodes[&start].sequence.clone();
        let mut cur = start;

        while let Some(next) = self.find_best_next_node(cur, visited) {
            if visited.contains(&next) {
                break;
            }
            path.push(next);
            if let Some(edge) = self.edges.get(&(cur, next)) {
                seq.push_str(&edge.overlap);
            } else {
                seq.push_str(&self.nodes[&next].sequence);
            }
            cur = next;
        }
        for &n in &path {
            visited.insert(n);
        }
        if path.len() > 1 {
            Some(AdaptiveContig {
                nodes: path.clone(),
                sequence: seq,
                avg_coverage: self.calculate_path_coverage(&path),
            })
        } else {
            None
        }
    }

    fn find_best_next_node(&self, cur: u64, visited: &HashSet<u64>) -> Option<u64> {
        self.edges
            .iter()
            .filter(|((from, to), _)| *from == cur && !visited.contains(to))
            .max_by_key(|(_, e)| e.weight)
            .map(|((_, to), _)| *to)
    }

    fn calculate_path_coverage(&self, path: &[u64]) -> f64 {
        path.iter()
            .map(|n| self.nodes[n].coverage as f64)
            .sum::<f64>()
            / path.len() as f64
    }

    // ========== Complexity / adaptive k ================================= //

    fn calculate_local_complexity(&self, seq: &str) -> Vec<f64> {
        let win = 50;
        let mut out = vec![0.0; seq.len()];
        let bytes = seq.as_bytes();

        for i in 0..seq.len() {
            let s = i.saturating_sub(win / 2);
            let e = (i + win / 2).min(seq.len());
            let window = &bytes[s..e];

            let mut count = [0u32; 4];
            for &b in window {
                match b {
                    b'A' | b'a' => count[0] += 1,
                    b'C' | b'c' => count[1] += 1,
                    b'G' | b'g' => count[2] += 1,
                    b'T' | b't' => count[3] += 1,
                    _ => {}
                }
            }
            let total = count.iter().sum::<u32>() as f64;
            if total == 0.0 {
                continue;
            }
            let entropy: f64 = count
                .iter()
                .filter(|&&c| c > 0)
                .map(|&c| {
                    let p = c as f64 / total;
                    -p * p.log2()
                })
                .sum();

            // Reversed: low entropy => high "complexity score".
            out[i] = 2.0 - entropy;
        }
        out
    }

    fn adaptive_k_size(&self, complexity: f64) -> usize {
        let norm = (complexity / 2.0).clamp(0.0, 1.0);
        self.base_k + ((self.max_k - self.base_k) as f64 * norm).round() as usize
    }

    // ========== k-mer helpers =========================================== //

    fn extract_kmer_at_k(&self, seq: &str, pos: usize, k: usize) -> Result<String> {
        if pos + k <= seq.len() {
            Ok(seq[pos..pos + k].to_owned())
        } else {
            bail!("k-mer range out of bounds");
        }
    }

    fn calculate_overlap(&self, k1: &str, k2: &str, k: usize) -> String {
        let max = (k - 1).min(k1.len()).min(k2.len());
        for i in (1..=max).rev() {
            if k1.ends_with(&k2[..i]) {
                return k2[i..].to_owned();
            }
        }
        k2.to_owned()
    }

    fn update_node(&mut self, hash: u64, seq: &str, complexity: f64) {
        self.nodes
            .entry(hash)
            .and_modify(|n| {
                n.coverage += 1;
                n.complexity = (n.complexity + complexity) / 2.0;
            })
            .or_insert(NodeData {
                sequence: seq.to_owned(),
                coverage: 1,
                complexity,
            });
    }

    fn update_edge(&mut self, key: (u64, u64), k: usize, overlap: String) {
        self.edges
            .entry(key)
            .and_modify(|e| e.weight += 1)
            .or_insert(EdgeData {
                k_size: k,
                weight: 1,
                overlap,
            });
    }

    fn add_coordinate_mapping(&mut self, _read_id: usize, _p1: usize, _p2: usize) {
        // Implement as needed
    }

    // ========== Minimiser stub (replace with real hash) ================= //

    /// Very small demonstrator: returns every `base_k`-step minimiser hash.
    fn extract_minimisers_with_positions(&self, seq: &str) -> Result<Vec<(usize, u64)>> {
        let mut v = Vec::new();
        for i in (0..=seq.len() - self.base_k).step_by(self.base_k) {
            let kmer = &seq[i..i + self.base_k];
            let hash = fxhash::hash64(kmer.as_bytes()); // fast 64-bit hash
            v.push((i, hash));
        }
        if v.len() < 2 {
            bail!("Need at least 2 minimisers to add a read");
        }
        Ok(v)
    }
}

// ------------------------------------------------------------------------- //
//  Public Contig struct
// ------------------------------------------------------------------------- //

#[derive(Clone, Debug)]
pub struct AdaptiveContig {
    pub nodes: Vec<u64>,
    pub sequence: String,
    pub avg_coverage: f64,
}

// ------------------------------------------------------------------------- //
//  Unit-tests
// ------------------------------------------------------------------------- //

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adaptive_k_bounds() {
        let g = AdaptiveGraph::new(15, 31);
        assert_eq!(g.adaptive_k_size(0.0), 15);
        assert_eq!(g.adaptive_k_size(2.0), 31);
    }

    #[test]
    fn complexity_ordering() {
        let g = AdaptiveGraph::new(15, 31);
        let low_entropy = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let high_entropy = "ACGTACGTACGTACGTACGTACGTACGTACGTACG";
        let c1 = g.calculate_local_complexity(low_entropy);
        let c2 = g.calculate_local_complexity(high_entropy);
        assert!(c1[10] > c2[10]); // low entropy => higher score
    }

    #[test]
    fn dummy_read_add() {
        let mut g = AdaptiveGraph::new(15, 31);
        let read = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        assert!(g.add_read_with_adaptive_k(0, read).is_ok());
    }
}
