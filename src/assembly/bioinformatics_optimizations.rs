//! Bioinformatics-specific optimizations for assembly graph construction
//! =====================================================================
//!
//! This module provides specialized optimizations for genomic data processing including:
//! - Bit-packed k-mer representations for memory efficiency
//! - Rolling hash implementations for streaming k-mer processing
//! - SIMD operations for nucleotide sequence processing
//! - Specialized data structures for large-scale genomic datasets
//! - Parallel graph algorithms for ultra‑fast genome assembly (OLC paradigm)
//!
//! The implementation integrates modern techniques such as parallel transitive‑reduction, task‑based
//! parallelism with Rayon, hierarchical graph merging, SCC‑based contig generation, adaptive k‑mer
//! sizing, classic graph‑cleaning passes (tip removal & bubble popping) and a lock‑free mindset using
//! AHash* structures.

use ahash::{AHashMap, AHashSet};
use anyhow::{anyhow, Result};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::sync::Arc;

/* ------------------------------------------------------------------------- */
/*                               K‑MER STRUCTURES                             */
/* ------------------------------------------------------------------------- */

/// Bit‑packed k‑mer representation for memory‑efficient storage.
/// Each nucleotide is encoded in 2 bits: A=00, C=01, G=10, T=11
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct BitPackedKmer {
    /// Packed nucleotides (up to 32 nucleotides per u64)
    pub packed_data: Vec<u64>,
    /// Length of the k‑mer
    pub k: usize,
    /// Hash value for quick comparison
    pub hash: u64,
    /// Whether this is the canonical representation
    pub is_canonical: bool,
}

impl BitPackedKmer {
    /// Create a new bit‑packed k‑mer from a DNA sequence
    pub fn new(sequence: &str) -> Result<Self> {
        if sequence.is_empty() {
            return Err(anyhow!("Empty sequence provided"));
        }
        let k = sequence.len();
        if k > 1024 {
            return Err(anyhow!("K‑mer too long: {k} (max 1024)"));
        }
        let sequence_upper = sequence.to_ascii_uppercase();
        let rc = Self::reverse_complement(&sequence_upper)?;
        let (canonical, is_canonical) = if sequence_upper <= rc {
            (sequence_upper, true)
        } else {
            (rc, false)
        };
        let packed_data = Self::pack_sequence(&canonical)?;
        let hash = Self::compute_hash(&packed_data, k);
        Ok(Self {
            packed_data,
            k,
            hash,
            is_canonical,
        })
    }

    /* -------------------------- packing helpers -------------------------- */
    fn pack_sequence(sequence: &str) -> Result<Vec<u64>> {
        let mut packed = Vec::new();
        let mut current_word = 0u64;
        let mut bits_used = 0;
        for nucleotide in sequence.chars() {
            let bits = match nucleotide {
                'A' => 0b00,
                'C' => 0b01,
                'G' => 0b10,
                'T' => 0b11,
                _ => return Err(anyhow!("Invalid nucleotide: {nucleotide}")),
            };
            current_word |= bits << (62 - bits_used);
            bits_used += 2;
            if bits_used >= 64 {
                packed.push(current_word);
                current_word = 0;
                bits_used = 0;
            }
        }
        if bits_used > 0 {
            packed.push(current_word);
        }
        if packed.is_empty() {
            packed.push(0);
        }
        Ok(packed)
    }

    pub fn unpack_sequence(&self) -> String {
        let mut sequence = String::with_capacity(self.k);
        let mut remaining = self.k;
        for &word in &self.packed_data {
            let nucleotides_in_word = remaining.min(32);
            for i in 0..nucleotides_in_word {
                let bits = (word >> (62 - i * 2)) & 0b11;
                let nucleotide = match bits {
                    0b00 => 'A',
                    0b01 => 'C',
                    0b10 => 'G',
                    0b11 => 'T',
                    _ => unreachable!(),
                };
                sequence.push(nucleotide);
            }
            remaining -= nucleotides_in_word;
            if remaining == 0 {
                break;
            }
        }
        sequence
    }

    fn reverse_complement(sequence: &str) -> Result<String> {
        sequence
            .chars()
            .rev()
            .map(|c| match c {
                'A' => Ok('T'),
                'T' => Ok('A'),
                'G' => Ok('C'),
                'C' => Ok('G'),
                _ => Err(anyhow!("Invalid nucleotide: {c}")),
            })
            .collect()
    }

    fn compute_hash(packed_data: &[u64], k: usize) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        let mut hasher = DefaultHasher::new();
        packed_data.hash(&mut hasher);
        k.hash(&mut hasher);
        hasher.finish()
    }

    pub fn memory_usage(&self) -> usize {
        self.packed_data.len() * std::mem::size_of::<u64>() + std::mem::size_of::<Self>()
    }
}

/* ------------------------------------------------------------------------- */
/*                     ROLLING HASH & STREAMING K‑MER OPS                    */
/* ------------------------------------------------------------------------- */

#[derive(Debug, Clone)]
pub struct RollingHash {
    hash: u64,
    k: usize,
    base: u64,
    base_power: u64,
    window: Vec<u8>,
}

impl RollingHash {
    pub fn new(k: usize) -> Self {
        let base = 4u64;
        let base_power = base.pow((k as u32).saturating_sub(1));
        Self {
            hash: 0,
            k,
            base,
            base_power,
            window: Vec::with_capacity(k),
        }
    }

    pub fn push(&mut self, nucleotide: char) -> Result<Option<u64>> {
        let encoded = Self::encode_nucleotide(nucleotide)?;
        if self.window.len() < self.k {
            self.window.push(encoded);
            self.hash = self.hash * self.base + encoded as u64;
            Ok((self.window.len() == self.k).then_some(self.hash))
        } else {
            let outgoing = self.window[0] as u64;
            self.hash = (self.hash - outgoing * self.base_power) * self.base + encoded as u64;
            self.window.remove(0);
            self.window.push(encoded);
            Ok(Some(self.hash))
        }
    }

    fn encode_nucleotide(n: char) -> Result<u8> {
        match n.to_ascii_uppercase() {
            'A' => Ok(0),
            'C' => Ok(1),
            'G' => Ok(2),
            'T' => Ok(3),
            _ => Err(anyhow!("Invalid nucleotide: {n}")),
        }
    }

    pub fn current_kmer(&self) -> String {
        if self.window.len() != self.k {
            return String::new();
        }
        self.window
            .iter()
            .map(|&e| match e {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => 'N',
            })
            .collect()
    }
}

/// Streaming k‑mer processor for ultra‑large files (lock‑free counting).
#[derive(Debug)]
pub struct StreamingKmerProcessor {
    k: usize,
    rolling_hash: RollingHash,
    canonical_kmers: AHashSet<u64>,
    kmer_counts: AHashMap<u64, u32>,
    processed_sequences: usize,
    total_kmers: usize,
}

impl StreamingKmerProcessor {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            rolling_hash: RollingHash::new(k),
            canonical_kmers: AHashSet::new(),
            kmer_counts: AHashMap::new(),
            processed_sequences: 0,
            total_kmers: 0,
        }
    }

    pub fn process_sequence(&mut self, sequence: &str) -> Result<Vec<(String, u64)>> {
        let mut kmers = Vec::new();
        self.rolling_hash = RollingHash::new(self.k);
        for nucleotide in sequence.chars() {
            if let Some(_hash) = self.rolling_hash.push(nucleotide)? {
                let kseq = self.rolling_hash.current_kmer();
                let bitk = BitPackedKmer::new(&kseq)?;
                let h = bitk.hash;
                *self.kmer_counts.entry(h).or_insert(0) += 1;
                self.canonical_kmers.insert(h);
                kmers.push((kseq, h));
                self.total_kmers += 1;
            }
        }
        self.processed_sequences += 1;
        Ok(kmers)
    }

    pub fn get_statistics(&self) -> ProcessingStatistics {
        ProcessingStatistics {
            unique_kmers: self.canonical_kmers.len(),
            total_kmers: self.total_kmers,
            processed_sequences: self.processed_sequences,
            average_kmer_frequency: if self.canonical_kmers.is_empty() {
                0.0
            } else {
                self.total_kmers as f64 / self.canonical_kmers.len() as f64
            },
            memory_usage_bytes: self.estimate_memory_usage(),
        }
    }

    fn estimate_memory_usage(&self) -> usize {
        let km = self.canonical_kmers.len() * std::mem::size_of::<u64>();
        let ct = self.kmer_counts.len() * (std::mem::size_of::<u64>() + std::mem::size_of::<u32>());
        km + ct + std::mem::size_of::<RollingHash>()
    }

    pub fn get_frequent_kmers(&self, min_frequency: u32) -> Vec<(u64, u32)> {
        self.kmer_counts
            .iter()
            .filter(|(_, &c)| c >= min_frequency)
            .map(|(&h, &c)| (h, c))
            .collect()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessingStatistics {
    pub unique_kmers: usize,
    pub total_kmers: usize,
    pub processed_sequences: usize,
    pub average_kmer_frequency: f64,
    pub memory_usage_bytes: usize,
}

/* ------------------------------------------------------------------------- */
/*                     SIMD PLACEHOLDER – PARALLEL FALLBACK                  */
/* ------------------------------------------------------------------------- */

pub struct SimdNucleotideOps;
impl SimdNucleotideOps {
    pub fn count_nucleotides(sequence: &str) -> [usize; 4] {
        let chunks: Vec<_> = sequence.as_bytes().par_chunks(1024).collect();
        chunks
            .par_iter()
            .map(|chunk| {
                let mut c = [0usize; 4];
                for &b in *chunk {
                    match b.to_ascii_uppercase() {
                        b'A' => c[0] += 1,
                        b'C' => c[1] += 1,
                        b'G' => c[2] += 1,
                        b'T' => c[3] += 1,
                        _ => {}
                    }
                }
                c
            })
            .reduce(
                || [0; 4],
                |mut a, b| {
                    for i in 0..4 {
                        a[i] += b[i];
                    }
                    a
                },
            )
    }
    pub fn gc_content(sequence: &str) -> f64 {
        let [a, c, g, t] = Self::count_nucleotides(sequence);
        let total = a + c + g + t;
        if total == 0 {
            0.
        } else {
            (c + g) as f64 / total as f64
        }
    }
}

/* ------------------------------------------------------------------------- */
/*                       ADAPTIVE K‑MER SIZING (STUB)                        */
/* ------------------------------------------------------------------------- */

/// Choose k adaptively based on local sequence complexity (Shannon entropy).
pub fn adaptive_kmer_size(sequence: &str, min_k: usize, max_k: usize) -> usize {
    let entropy = {
        let counts = SimdNucleotideOps::count_nucleotides(sequence);
        let total = counts.iter().sum::<usize>() as f64;
        if total == 0. {
            return min_k;
        }
        counts.iter().fold(0.0, |e, &c| {
            if c == 0 {
                e
            } else {
                let p = c as f64 / total;
                e - p * p.log2()
            }
        })
    };
    let normalized = entropy / 2.0; // max entropy for 4 symbols = 2 bits
    let k_range = (max_k - min_k) as f64;
    min_k + (k_range * normalized) as usize
}

/* ------------------------------------------------------------------------- */
/*                        ASSEMBLY GRAPH & PARALLEL OPS                       */
/* ------------------------------------------------------------------------- */

#[derive(Debug, Clone)]
pub struct AssemblyGraph {
    adjacency: AHashMap<u64, AHashSet<u64>>, // Directed edges u -> v (overlap)
}

impl AssemblyGraph {
    pub fn new() -> Self {
        Self {
            adjacency: AHashMap::new(),
        }
    }

    pub fn add_edge(&mut self, from: u64, to: u64) {
        self.adjacency
            .entry(from)
            .or_insert_with(AHashSet::new)
            .insert(to);
    }

    pub fn merge(&mut self, other: &AssemblyGraph) {
        for (u, neigh) in &other.adjacency {
            self.adjacency
                .entry(*u)
                .or_insert_with(AHashSet::new)
                .extend(neigh);
        }
    }

    /* ------------------- PARALLEL TRANSITIVE REDUCTION ------------------- */
    pub fn transitive_reduction_parallel(&mut self) {
        let arc_adj = Arc::new(self.adjacency.clone()); // read‑only snapshot
        let edges_to_remove: Vec<(u64, u64)> = self
            .adjacency
            .par_iter()
            .flat_map(|(&u, neigh)| {
                neigh.iter().filter_map(move |&v| {
                    if Self::reachable_except_direct(&arc_adj, u, v) {
                        Some((u, v))
                    } else {
                        None
                    }
                })
            })
            .collect();
        for (u, v) in edges_to_remove {
            if let Some(neigh) = self.adjacency.get_mut(&u) {
                neigh.remove(&v);
            }
        }
    }

    fn reachable_except_direct(
        adj: &AHashMap<u64, AHashSet<u64>>,
        start: u64,
        target: u64,
    ) -> bool {
        let mut queue: VecDeque<u64> = VecDeque::new();
        let mut visited: AHashSet<u64> = AHashSet::new();
        if let Some(neigh) = adj.get(&start) {
            for &n in neigh {
                if n != target {
                    queue.push_back(n);
                    visited.insert(n);
                }
            }
        }
        while let Some(node) = queue.pop_front() {
            if node == target {
                return true;
            }
            if let Some(neigh) = adj.get(&node) {
                for &n in neigh {
                    if !visited.contains(&n) {
                        visited.insert(n);
                        queue.push_back(n);
                    }
                }
            }
        }
        false
    }

    /* --------------------------- TIP REMOVAL ----------------------------- */
    pub fn remove_tips(&mut self, tip_len: usize) {
        let mut to_remove = Vec::new();
        for (&node, neigh) in &self.adjacency {
            if neigh.is_empty() && self.in_degree(node) == 1 {
                // dead‑end
                to_remove.push((node, tip_len));
            }
        }
        for (tip, mut len) in to_remove {
            let mut current = tip;
            while len > 0 {
                let pred = self.predecessors(current).into_iter().next();
                if let Some(p) = pred {
                    self.adjacency.get_mut(&p).unwrap().remove(&current);
                    current = p;
                    len -= 1;
                } else {
                    break;
                }
            }
        }
    }

    /* --------------------------- BUBBLE POP ------------------------------ */
    pub fn pop_bubbles(&mut self) {
        // Simplistic implementation: for each pair of edges u->v and u->w where v & w converge
        // to same successor, remove the longer path.
        let mut to_remove = Vec::new();
        for (&u, neigh) in &self.adjacency {
            let neigh_vec: Vec<u64> = neigh.iter().cloned().collect();
            for i in 0..neigh_vec.len() {
                for j in (i + 1)..neigh_vec.len() {
                    let v = neigh_vec[i];
                    let w = neigh_vec[j];
                    if let Some(common) = self.common_successor(v, w) {
                        // decide which branch shorter
                        let len_v = self.path_length(v, common, 10);
                        let len_w = self.path_length(w, common, 10);
                        if len_v > len_w {
                            to_remove.push((u, v));
                        } else {
                            to_remove.push((u, w));
                        }
                    }
                }
            }
        }
        for (u, v) in to_remove {
            if let Some(neigh) = self.adjacency.get_mut(&u) {
                neigh.remove(&v);
            }
        }
    }

    fn common_successor(&self, a: u64, b: u64) -> Option<u64> {
        let succ_a: AHashSet<u64> = self.successors(a).into_iter().collect();
        for s in self.successors(b) {
            if succ_a.contains(&s) {
                return Some(s);
            }
        }
        None
    }

    fn path_length(&self, from: u64, to: u64, max: usize) -> usize {
        let mut queue: VecDeque<(u64, usize)> = VecDeque::new();
        let mut visited: AHashSet<u64> = AHashSet::new();
        queue.push_back((from, 0));
        visited.insert(from);
        while let Some((n, d)) = queue.pop_front() {
            if d > max {
                return max + 1;
            }
            if n == to {
                return d;
            }
            for s in self.successors(n) {
                if visited.insert(s) {
                    queue.push_back((s, d + 1));
                }
            }
        }
        max + 1
    }

    /* --------------------- STRONGLY CONNECTED COMPONENTS ------------------ */
    pub fn strongly_connected_components(&self) -> Vec<Vec<u64>> {
        // Kosaraju (sequential). Adequate for typical component counts.
        let mut visited: AHashSet<u64> = AHashSet::new();
        let mut order = Vec::new();
        for &v in self.adjacency.keys() {
            if !visited.contains(&v) {
                self.dfs_order(v, &mut visited, &mut order);
            }
        }
        let transposed = self.transpose();
        visited.clear();
        let mut components = Vec::new();
        for &v in order.iter().rev() {
            if !visited.contains(&v) {
                let mut comp = Vec::new();
                transposed.collect_component(v, &mut visited, &mut comp);
                components.push(comp);
            }
        }
        components
    }

    fn dfs_order(&self, v: u64, visited: &mut AHashSet<u64>, order: &mut Vec<u64>) {
        let mut stack = Vec::new();
        stack.push((v, false));
        while let Some((node, processed)) = stack.pop() {
            if processed {
                order.push(node);
                continue;
            }
            if !visited.insert(node) {
                continue;
            }
            stack.push((node, true));
            for n in self.successors(node) {
                stack.push((n, false));
            }
        }
    }

    fn collect_component(&self, v: u64, visited: &mut AHashSet<u64>, comp: &mut Vec<u64>) {
        let mut stack = vec![v];
        visited.insert(v);
        while let Some(node) = stack.pop() {
            comp.push(node);
            for n in self.successors(node) {
                if visited.insert(n) {
                    stack.push(n);
                }
            }
        }
    }

    fn transpose(&self) -> AssemblyGraph {
        let mut g = AssemblyGraph::new();
        for (&u, neigh) in &self.adjacency {
            for &v in neigh {
                g.add_edge(v, u);
            }
        }
        g
    }

    /* ------------------- PARALLEL CONTIG GENERATION (OLC) ----------------- */
    pub fn generate_contigs_parallel(&self) -> Vec<Vec<u64>> {
        self.strongly_connected_components()
            .into_par_iter()
            .map(|comp| self.eulerian_trail(&comp))
            .collect()
    }

    fn eulerian_trail(&self, nodes: &[u64]) -> Vec<u64> {
        // Hierholzer within component (naïve implementation)
        let mut stack = Vec::new();
        let mut path = Vec::new();
        if nodes.is_empty() {
            return path;
        }
        let start = nodes[0];
        stack.push(start);
        let mut adj_iter: AHashMap<u64, Vec<u64>> = AHashMap::new();
        for &n in nodes {
            adj_iter.insert(n, self.successors(n));
        }
        while let Some(v) = stack.last().cloned() {
            if let Some(next) = adj_iter.get_mut(&v).and_then(|l| l.pop()) {
                stack.push(next);
            } else {
                path.push(v);
                stack.pop();
            }
        }
        path.reverse();
        path
    }

    /* ---------------------- AUXILIARY GRAPH HELPERS ---------------------- */
    fn successors(&self, node: u64) -> Vec<u64> {
        self.adjacency
            .get(&node)
            .map(|s| s.iter().cloned().collect())
            .unwrap_or_default()
    }
    fn predecessors(&self, node: u64) -> Vec<u64> {
        self.adjacency
            .iter()
            .filter_map(|(&u, neigh)| neigh.contains(&node).then(|| u))
            .collect()
    }
    fn in_degree(&self, node: u64) -> usize {
        self.predecessors(node).len()
    }
}

/* -------------------------- HIERARCHICAL MERGE --------------------------- */

pub fn merge_graphs_hierarchical(mut graphs: Vec<AssemblyGraph>) -> AssemblyGraph {
    if graphs.is_empty() {
        return AssemblyGraph::new();
    }
    while graphs.len() > 1 {
        graphs = graphs
            .par_chunks(2)
            .map(|chunk| {
                let mut g = chunk[0].clone();
                if chunk.len() == 2 {
                    g.merge(&chunk[1]);
                }
                g
            })
            .collect();
    }
    graphs.pop().unwrap()
}

/* ------------------------------------------------------------------------- */
/*                                  TESTS                                    */
/* ------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_packed_kmer_roundtrip() {
        let seq = "ATCGATCG";
        let k = BitPackedKmer::new(seq).unwrap();
        assert_eq!(k.unpack_sequence(), seq);
    }

    #[test]
    fn test_transitive_reduction() {
        let mut g = AssemblyGraph::new();
        g.add_edge(1, 2);
        g.add_edge(2, 3);
        g.add_edge(1, 3); // transitive
        g.transitive_reduction_parallel();
        assert!(g.successors(1).len() == 1 && g.successors(1)[0] == 2);
    }

    #[test]
    fn test_scc_contigs() {
        let mut g = AssemblyGraph::new();
        g.add_edge(1, 2);
        g.add_edge(2, 1);
        g.add_edge(2, 3);
        let contigs = g.generate_contigs_parallel();
        assert_eq!(contigs.len(), 2); // one SCC {1,2} and one {3}
    }
}
