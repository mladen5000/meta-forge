//! Assembly Graph Construction and Simplification
//! ===========================================
//! A fully‑featured assembly graph builder with adaptive *k*-mer sizing,
//! graph‑level simplification, contig generation, and rich statistics.
//!
//! **Highlights**
//! * Adaptive *k*-mer size chosen per read‑chunk (`complexity_threshold`
//!   heuristics).
//! * Multi‑threaded extraction / assembly with `rayon` thread‑pool.
//! * Graph simplification passes: tip clipping, bubble popping, low‑confidence
//!   edge trimming, and linear‑chain compression.
//! * Contig extraction (linear / circular / scaffold) with assembly stats &
//!   FASTA + GFA exporters.
//!
//! **Note**  
//! Everything compiles on stable *Rust 1.77* assuming that the
//! `core_data_structures` module provides the requisite
//! `GraphNode`, `GraphEdge`, `GraphFragment`, `CanonicalKmer`, `CorrectedRead`
//! and helpers such as `calculate_sequence_complexity`.
//!
//! ```text
//! ├── AssemblyGraphBuilder
//! │   ├── create_assembly_chunks()
//! │   ├── build_fragment()
//! │   ├── merge_fragments()
//! │   └── simplify_graph()
//! └── AssemblyGraph
//!     ├── generate_contigs()
//!     ├── write_contigs_fasta()
//!     └── write_graph_gfa()
//! ```
//!
//! All public functions and data‑types are documented; run
//! `cargo doc --open` for an HTML view.

use ahash::{AHashMap, AHashSet};
use anyhow::{Result, anyhow};
use petgraph::algo::tarjan_scc;
use petgraph::{Directed, Graph, graph::NodeIndex};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};

use crate::core_data_structures::*;

/* ------------------------------------------------------------------------- */
/*                           ASSEMBLY GRAPH BUILDER                          */
/* ------------------------------------------------------------------------- */

/// Builds an [`AssemblyGraph`] from a slice of [`CorrectedRead`]s.
///
/// The builder operates in four high‑level phases:
///
/// 1. **Chunking / adaptive `k`** – reads are binned into 1 000‑read chunks and
///    an optimal *k* is chosen per chunk by measuring sequence complexity.
/// 2. **Per‑chunk graph construction** – each chunk is converted into an
///    independent [`GraphFragment`] in parallel.
/// 3. **Fragment merge** – all fragments are merged into a single graph and
///    converted into a `petgraph::Graph` for algorithmic convenience.
/// 4. **Graph simplification** – iterative clean‑up passes remove artefacts
///    (tips, bubbles, low‑confidence edges, linear chains).
#[derive(Debug)]
pub struct AssemblyGraphBuilder {
    /// Base *k*-mer size (smallest allowed).
    base_k: usize,
    /// Maximum *k*-mer size when complexity → 1.0.
    max_k: usize,
    /// Minimum per‑node coverage to keep a vertex.
    min_coverage: u32,
    /// Internal heuristic constant: when average complexity exceeds
    /// `complexity_threshold` the maximum *k* is selected.
    complexity_threshold: f64,
    /// Shared thread‑pool for all parallel work.
    thread_pool: rayon::ThreadPool,
}

impl AssemblyGraphBuilder {
    /// Construct a new builder.
    ///
    /// * `base_k` – smallest *k*-mer size.
    /// * `max_k` – largest *k*-mer size (inclusive).
    /// * `min_coverage` – nodes with coverage \< `min_coverage` are dropped.
    /// * `num_threads` – rayon thread‑pool size.
    pub fn new(base_k: usize, max_k: usize, min_coverage: u32, num_threads: usize) -> Result<Self> {
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .thread_name(|i| format!("assembly-{i}"))
            .build()?;

        Ok(Self {
            base_k,
            max_k,
            min_coverage,
            complexity_threshold: 0.7,
            thread_pool,
        })
    }

    /// Build an [`AssemblyGraph`] from corrected reads.
    pub fn build_graph(&self, reads: &[CorrectedRead]) -> Result<AssemblyGraph> {
        println!("🧬 Building assembly graph from {} reads", reads.len());

        // 1. chunking + adaptive k
        let chunks = self
            .thread_pool
            .install(|| self.create_assembly_chunks(reads))?;

        // 2. per‑chunk fragment construction
        let fragments: Result<Vec<_>> = chunks
            .par_iter()
            .map(|chunk| self.build_fragment(chunk))
            .collect();
        let fragments = fragments?;

        // 3. merge
        let merged = self.merge_fragments(fragments)?;

        // 4. simplify
        self.simplify_graph(merged)
    }

    /* ------------------- phase 1: chunking & adaptive k ------------------- */

    /// Split `reads` into fixed‑size chunks and calculate the optimal *k* per
    /// chunk based on average sequence complexity (Shannon entropy proxy).
    fn create_assembly_chunks(&self, reads: &[CorrectedRead]) -> Result<Vec<AssemblyChunk>> {
        const CHUNK_SIZE: usize = 1_000;
        let mut chunks = Vec::new();

        for (chunk_id, batch) in reads.chunks(CHUNK_SIZE).enumerate() {
            let complexity = self.analyze_batch_complexity(batch);
            let k = self.select_optimal_k(complexity);

            let mut chunk = AssemblyChunk::new(chunk_id, k);
            for read in batch {
                chunk.add_read(read.clone())?;
            }
            chunk.finalize();
            chunks.push(chunk);
        }
        Ok(chunks)
    }

    /// Compute mean complexity for a slice of reads.
    fn analyze_batch_complexity(&self, reads: &[CorrectedRead]) -> f64 {
        let sum: f64 = reads
            .par_iter()
            .map(|r| calculate_sequence_complexity(&r.corrected))
            .sum();
        sum / reads.len().max(1) as f64
    }

    /// Map a complexity score (0–1) to a *k*-mer size in `[base_k, max_k]`.
    #[inline]
    fn select_optimal_k(&self, complexity: f64) -> usize {
        let span = self.max_k - self.base_k;
        let factor = complexity.clamp(0.0, 1.0);
        self.base_k + ((factor * span as f64).round() as usize)
    }

    /* ---------------- phase 2: fragment construction per chunk ----------- */

    /// Build a [`GraphFragment`] from a single [`AssemblyChunk`], running the
    /// local clean‑up passes (low‑coverage filter, tip/bubble detection, edge
    /// weighting).
    fn build_fragment(&self, chunk: &AssemblyChunk) -> Result<GraphFragment> {
        let mut fragment = chunk.graph_fragment.clone();
        self.filter_low_coverage_nodes(&mut fragment)?;
        self.detect_structural_features(&mut fragment)?;
        self.calculate_edge_weights(&mut fragment)?;
        Ok(fragment)
    }

    /// Remove vertices whose `coverage` \< `min_coverage`.
    fn filter_low_coverage_nodes(&self, fragment: &mut GraphFragment) -> Result<()> {
        let initial = fragment.nodes.len();
        let to_remove: Vec<u64> = fragment
            .nodes
            .iter()
            .filter(|(_, n)| n.coverage < self.min_coverage)
            .map(|(&h, _)| h)
            .collect();
        for h in &to_remove {
            fragment.nodes.remove(h);
        }
        fragment
            .edges
            .retain(|e| !to_remove.contains(&e.from_hash) && !to_remove.contains(&e.to_hash));

        println!("   Filtered {} low‑coverage nodes", to_remove.len());
        Ok(())
    }

    /// Mark tips & bubbles (`node_type` metadata) for downstream heuristics.
    fn detect_structural_features(&self, fragment: &mut GraphFragment) -> Result<()> {
        let tips = fragment.find_tips();
        for h in &tips {
            if let Some(node) = fragment.nodes.get_mut(h) {
                node.node_type = NodeType::Tip;
            }
        }
        let bubbles = fragment.find_bubbles();
        for b in &bubbles {
            for path in &b.alternative_paths {
                for h in path {
                    if let Some(node) = fragment.nodes.get_mut(h) {
                        node.node_type = NodeType::Bubble;
                    }
                }
            }
        }
        println!(
            "   Detected {} tips · {} bubbles",
            tips.len(),
            bubbles.len()
        );
        Ok(())
    }

    /// Compute `edge.confidence` as √(support × min_cov) / 100 ∈ [0.1, 1].
    fn calculate_edge_weights(&self, fragment: &mut GraphFragment) -> Result<()> {
        for e in &mut fragment.edges {
            let from_cov = fragment
                .nodes
                .get(&e.from_hash)
                .map(|n| n.coverage)
                .unwrap_or(1);
            let to_cov = fragment
                .nodes
                .get(&e.to_hash)
                .map(|n| n.coverage)
                .unwrap_or(1);
            let min_cov = from_cov.min(to_cov);
            let support = e.supporting_reads.len() as f64;
            e.confidence = (support * min_cov as f64).sqrt() / 100.0;
            e.confidence = e.confidence.clamp(0.1, 1.0);
        }
        Ok(())
    }

    /* ------------------------- phase 3: merge ---------------------------- */

    /// Merge per‑chunk fragments into a single assembly graph.
    fn merge_fragments(&self, fragments: Vec<GraphFragment>) -> Result<AssemblyGraph> {
        println!("🔗 Merging {} fragments", fragments.len());
        let mut merged = GraphFragment::new(0);
        for f in fragments {
            merged.merge_with(f)?;
        }
        let pet = self.convert_to_petgraph(&merged)?;
        Ok(AssemblyGraph {
            graph_fragment: merged,
            petgraph: pet,
            contigs: Vec::new(),
            assembly_stats: AssemblyStats::default(),
        })
    }

    /// Translate a [`GraphFragment`] into `petgraph::Graph` (**O(V + E)**).
    fn convert_to_petgraph(
        &self,
        fragment: &GraphFragment,
    ) -> Result<Graph<u64, EdgeWeight, Directed>> {
        let mut g = Graph::<u64, EdgeWeight, Directed>::new();
        let mut idx = AHashMap::<u64, NodeIndex>::new();

        for &h in fragment.nodes.keys() {
            idx.insert(h, g.add_node(h));
        }
        for e in &fragment.edges {
            if let (Some(&from), Some(&to)) = (idx.get(&e.from_hash), idx.get(&e.to_hash)) {
                g.add_edge(
                    from,
                    to,
                    EdgeWeight {
                        weight: e.weight,
                        confidence: e.confidence,
                    },
                );
            }
        }
        println!(
            "   Petgraph: {} nodes · {} edges",
            g.node_count(),
            g.edge_count()
        );
        Ok(g)
    }

    /* ---------------------- phase 4: simplification ---------------------- */

    /// Run all clean‑up passes on `assembly_graph` and return the pruned graph.
    fn simplify_graph(&self, mut g: AssemblyGraph) -> Result<AssemblyGraph> {
        println!("✂️  Simplifying graph…");
        self.remove_tips(&mut g)?;
        self.pop_bubbles(&mut g)?;
        self.merge_linear_chains(&mut g)?;
        self.remove_low_confidence_edges(&mut g)?;
        Ok(g)
    }

    /* --------------- tip removal (dead‑end pruning) ---------------------- */

    fn remove_tips(&self, g: &mut AssemblyGraph) -> Result<()> {
        const MAX_TIP_LEN: usize = 100;
        const MIN_RATIO: f64 = 0.1;
        let tips = g.graph_fragment.find_tips();
        let mut removed = 0;

        for h in &tips {
            if let Some(node) = g.graph_fragment.nodes.get(h) {
                if node.kmer.len() <= MAX_TIP_LEN {
                    let neighbor_cov = self.get_neighbor_coverage(&g.graph_fragment, *h);
                    if neighbor_cov > 0.0 && (node.coverage as f64 / neighbor_cov) < MIN_RATIO {
                        g.graph_fragment.nodes.remove(h);
                        g.graph_fragment
                            .edges
                            .retain(|e| e.from_hash != *h && e.to_hash != *h);
                        removed += 1;
                    }
                }
            }
        }
        println!("   Removed {removed} tips");
        Ok(())
    }

    fn get_neighbor_coverage(&self, frag: &GraphFragment, h: u64) -> f64 {
        let mut covs = Vec::new();
        for e in &frag.edges {
            if e.from_hash == h {
                if let Some(n) = frag.nodes.get(&e.to_hash) {
                    covs.push(n.coverage as f64);
                }
            } else if e.to_hash == h {
                if let Some(n) = frag.nodes.get(&e.from_hash) {
                    covs.push(n.coverage as f64);
                }
            }
        }
        covs.iter().copied().sum::<f64>() / covs.len().max(1) as f64
    }

    /* ----------------- bubble popping (simple 2‑path bubbles) ------------ */

    fn pop_bubbles(&self, g: &mut AssemblyGraph) -> Result<()> {
        let bubbles = g.graph_fragment.find_bubbles();
        let mut popped = 0;

        for b in bubbles {
            if b.bubble_type == BubbleType::Simple && b.alternative_paths.len() == 2 {
                let p1_cov =
                    self.calculate_path_coverage(&g.graph_fragment, &b.alternative_paths[0]);
                let p2_cov =
                    self.calculate_path_coverage(&g.graph_fragment, &b.alternative_paths[1]);
                let remove = if p1_cov > p2_cov { 1 } else { 0 };
                let doomed = &b.alternative_paths[remove];

                for h in doomed {
                    g.graph_fragment.nodes.remove(h);
                }
                g.graph_fragment
                    .edges
                    .retain(|e| !doomed.contains(&e.from_hash) && !doomed.contains(&e.to_hash));
                popped += 1;
            }
        }
        println!("   Popped {popped} bubbles");
        Ok(())
    }

    fn calculate_path_coverage(&self, frag: &GraphFragment, path: &[u64]) -> f64 {
        if path.is_empty() {
            return 0.0;
        }
        let tot: u32 = path
            .iter()
            .filter_map(|h| frag.nodes.get(h))
            .map(|n| n.coverage)
            .sum();
        tot as f64 / path.len() as f64
    }

    /* ------------ linear‑chain compression (unitig collapse) ------------ */

    fn merge_linear_chains(&self, g: &mut AssemblyGraph) -> Result<()> {
        let chains = self.find_linear_chains(&g.graph_fragment);
        let mut merged = 0;

        for chain in chains {
            if chain.len() <= 1 {
                continue;
            }
            let new_node = self.merge_chain_nodes(&g.graph_fragment, &chain)?;
            for h in &chain[1..] {
                g.graph_fragment.nodes.remove(h);
            }
            if let Some(first) = g.graph_fragment.nodes.get_mut(&chain[0]) {
                *first = new_node;
            }
            g.graph_fragment
                .edges
                .retain(|e| !chain[1..].contains(&e.from_hash) && !chain[1..].contains(&e.to_hash));
            merged += 1;
        }
        println!("   Merged {merged} linear chains");
        Ok(())
    }

    fn find_linear_chains(&self, frag: &GraphFragment) -> Vec<Vec<u64>> {
        let adj = frag.get_adjacency_list();
        let mut indeg = AHashMap::<u64, usize>::new();
        for neigh in adj.values() {
            for n in neigh {
                *indeg.entry(*n).or_default() += 1;
            }
        }
        let mut chains = Vec::new();
        let mut visited = AHashSet::new();

        for (&n, neigh) in &adj {
            if visited.contains(&n) {
                continue;
            }
            let in_ = *indeg.get(&n).unwrap_or(&0);
            if in_ != 1 || neigh.len() != 1 {
                let chain = self.trace_linear_chain(n, &adj, &indeg, &mut visited);
                if chain.len() > 1 {
                    chains.push(chain);
                }
            }
        }
        chains
    }

    fn trace_linear_chain(
        &self,
        start: u64,
        adj: &AHashMap<u64, Vec<u64>>,
        indeg: &AHashMap<u64, usize>,
        visited: &mut AHashSet<u64>,
    ) -> Vec<u64> {
        let mut chain = vec![start];
        visited.insert(start);
        let mut cur = start;

        while let Some(neigh) = adj.get(&cur) {
            if neigh.len() == 1 {
                let next = neigh[0];
                if *indeg.get(&next).unwrap_or(&0) == 1 && !visited.contains(&next) {
                    chain.push(next);
                    visited.insert(next);
                    cur = next;
                    continue;
                }
            }
            break;
        }
        chain
    }

    fn merge_chain_nodes(&self, frag: &GraphFragment, chain: &[u64]) -> Result<GraphNode> {
        let first = frag
            .nodes
            .get(&chain[0])
            .ok_or_else(|| anyhow!("first node missing"))?;
        let mut seq = first.kmer.sequence.clone();
        let mut cov = first.coverage;
        let mut pos = first.read_positions.clone();
        let k = first.kmer_size;

        for h in &chain[1..] {
            if let Some(n) = frag.nodes.get(h) {
                let overlap = k - 1;
                if n.kmer.sequence.len() > overlap {
                    seq.push_str(&n.kmer.sequence[overlap..]);
                }
                cov += n.coverage;
                pos.extend(n.read_positions.iter().cloned());
            }
        }
        let merged_kmer = CanonicalKmer::new(&seq)?;
        let mut node = GraphNode::new(merged_kmer, seq.len());
        node.coverage = cov;
        node.read_positions = pos;
        node.complexity_score = calculate_sequence_complexity(&seq);
        node.update_node_type();
        Ok(node)
    }

    /* ----------- low‑confidence edge removal (global pass) -------------- */

    fn remove_low_confidence_edges(&self, g: &mut AssemblyGraph) -> Result<()> {
        const THRESH: f64 = 0.3;
        let before = g.graph_fragment.edges.len();
        g.graph_fragment.edges.retain(|e| e.confidence >= THRESH);
        println!(
            "   Dropped {} low‑confidence edges",
            before - g.graph_fragment.edges.len()
        );
        Ok(())
    }
}

/* ------------------------------------------------------------------------- */
/*                             ASSEMBLY GRAPH                                */
/* ------------------------------------------------------------------------- */

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyGraph {
    pub graph_fragment: GraphFragment,
    #[serde(skip)]
    pub petgraph: Graph<u64, EdgeWeight, Directed>,
    pub contigs: Vec<Contig>,
    pub assembly_stats: AssemblyStats,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EdgeWeight {
    pub weight: u32,
    pub confidence: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Contig {
    pub id: usize,
    pub sequence: String,
    pub coverage: f64,
    pub length: usize,
    pub node_path: Vec<u64>,
    pub contig_type: ContigType,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ContigType {
    Linear,
    Circular,
    Scaffold,
    Fragment,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AssemblyStats {
    pub total_contigs: usize,
    pub total_length: usize,
    pub longest_contig: usize,
    pub n50: usize,
    pub n90: usize,
    pub mean_coverage: f64,
    pub gc_content: f64,
    pub gaps: usize,
}

/* ................................ impl .................................. */

impl AssemblyGraph {
    /// Derive contigs from the simplified graph, update `assembly_stats`.
    pub fn generate_contigs(&mut self) -> Result<()> {
        println!("📝 Generating contigs…");
        let sccs = tarjan_scc(&self.petgraph);
        let mut contigs = Vec::new();
        let mut cid = 0;

        for comp in sccs {
            if comp.len() == 1 {
                if let Some(c) = self.build_linear_contig(cid, &comp)? {
                    contigs.push(c);
                    cid += 1;
                }
            } else if let Some(c) = self.build_complex_contig(cid, &comp)? {
                contigs.push(c);
                cid += 1;
            }
        }
        contigs.sort_by(|a, b| b.length.cmp(&a.length));
        for (i, c) in contigs.iter_mut().enumerate() {
            c.id = i;
        }
        self.contigs = contigs;
        self.calculate_assembly_stats();
        println!("   {} contigs generated", self.contigs.len());
        Ok(())
    }

    /* --- helpers for contig construction (linear vs complex) ------------ */

    fn build_linear_contig(&self, id: usize, comp: &[NodeIndex]) -> Result<Option<Contig>> {
        if let Some(&hash) = self.petgraph.node_weight(comp[0]) {
            if let Some(node) = self.graph_fragment.nodes.get(&hash) {
                return Ok(Some(Contig {
                    id,
                    sequence: node.kmer.sequence.clone(),
                    coverage: node.coverage as f64,
                    length: node.kmer.sequence.len(),
                    node_path: vec![hash],
                    contig_type: ContigType::Linear,
                }));
            }
        }
        Ok(None)
    }

    fn build_complex_contig(&self, id: usize, comp: &[NodeIndex]) -> Result<Option<Contig>> {
        let path = self.find_eulerian_path(comp)?;
        if let Some(hash_path) = path {
            let sequence = self.reconstruct_sequence_from_path(&hash_path)?;
            let coverage = self.calculate_path_coverage_from_hashes(&hash_path);
            let ctype = if self.is_circular_path(&hash_path) {
                ContigType::Circular
            } else {
                ContigType::Scaffold
            };
            let k = self
                .graph_fragment
                .nodes
                .values()
                .next()
                .map(|n| n.kmer_size)
                .unwrap_or(0);
            return Ok(Some(Contig {
                id,
                sequence,
                coverage,
                length: hash_path.len() * k,
                node_path: hash_path,
                contig_type: ctype,
            }));
        }
        Ok(None)
    }

    /* --- Eulerian path (placeholder simple walk) ------------------------ */

    fn find_eulerian_path(&self, comp: &[NodeIndex]) -> Result<Option<Vec<u64>>> {
        if comp.is_empty() {
            return Ok(None);
        }
        let mut path = Vec::with_capacity(comp.len());
        for n in comp {
            if let Some(&h) = self.petgraph.node_weight(*n) {
                path.push(h);
            }
        }
        Ok(Some(path))
    }

    fn reconstruct_sequence_from_path(&self, path: &[u64]) -> Result<String> {
        if path.is_empty() {
            return Ok(String::new());
        }
        let first = self
            .graph_fragment
            .nodes
            .get(&path[0])
            .ok_or_else(|| anyhow!("node missing"))?;
        let mut seq = first.kmer.sequence.clone();
        let k = first.kmer_size;

        for h in &path[1..] {
            if let Some(n) = self.graph_fragment.nodes.get(h) {
                let overlap = k - 1;
                seq.push_str(&n.kmer.sequence[overlap..]);
            }
        }
        Ok(seq)
    }

    fn calculate_path_coverage_from_hashes(&self, path: &[u64]) -> f64 {
        if path.is_empty() {
            return 0.0;
        }
        let tot: u32 = path
            .iter()
            .filter_map(|h| self.graph_fragment.nodes.get(h))
            .map(|n| n.coverage)
            .sum();
        tot as f64 / path.len() as f64
    }

    fn is_circular_path(&self, path: &[u64]) -> bool {
        path.len() >= 3
            && self
                .graph_fragment
                .edges
                .iter()
                .any(|e| e.from_hash == path[path.len() - 1] && e.to_hash == path[0])
    }

    /* --- assembly stats -------------------------------------------------- */

    fn calculate_assembly_stats(&mut self) {
        if self.contigs.is_empty() {
            return;
        }
        self.assembly_stats.total_contigs = self.contigs.len();
        self.assembly_stats.total_length = self.contigs.iter().map(|c| c.length).sum();
        self.assembly_stats.longest_contig =
            self.contigs.iter().map(|c| c.length).max().unwrap_or(0);

        let mut lens: Vec<_> = self.contigs.iter().map(|c| c.length).collect();
        lens.sort_unstable_by(|a, b| b.cmp(a));

        let total = self.assembly_stats.total_length;
        let mut cum = 0;
        for &l in &lens {
            cum += l;
            if self.assembly_stats.n50 == 0 && cum >= total / 2 {
                self.assembly_stats.n50 = l;
            }
            if self.assembly_stats.n90 == 0 && cum >= total * 9 / 10 {
                self.assembly_stats.n90 = l;
                break;
            }
        }
        let cov_sum: f64 = self
            .contigs
            .iter()
            .map(|c| c.coverage * c.length as f64)
            .sum();
        self.assembly_stats.mean_coverage = cov_sum / total as f64;

        let gc: usize = self
            .contigs
            .iter()
            .map(|c| {
                c.sequence
                    .chars()
                    .filter(|&ch| ch == 'G' || ch == 'C')
                    .count()
            })
            .sum();
        self.assembly_stats.gc_content = gc as f64 / total as f64;
    }

    /* --- exporters ------------------------------------------------------- */

    /// Write contigs in FASTA. 80 bp line‑wrap.
    pub fn write_contigs_fasta(&self, path: &std::path::Path) -> Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(path)?;
        for c in &self.contigs {
            writeln!(f, ">contig_{} len={} cov={:.2}", c.id, c.length, c.coverage)?;
            for chunk in c.sequence.as_bytes().chunks(80) {
                writeln!(f, "{}", std::str::from_utf8(chunk)?)?;
            }
        }
        println!("✅ FASTA written → {}", path.display());
        Ok(())
    }

    /// Export the graph in GFA v1.
    pub fn write_graph_gfa(&self, path: &std::path::Path) -> Result<()> {
        use std::io::Write;
        let mut f = std::fs::File::create(path)?;
        writeln!(f, "H\tVN:Z:1.0")?;

        for (h, n) in &self.graph_fragment.nodes {
            writeln!(
                f,
                "S\t{}\t{}\tLN:i:{}\tRC:i:{}",
                h,
                n.kmer.sequence,
                n.kmer.sequence.len(),
                n.coverage
            )?;
        }
        for e in &self.graph_fragment.edges {
            writeln!(
                f,
                "L\t{}\t+\t{}\t+\t{}M\tRC:i:{}",
                e.from_hash, e.to_hash, e.overlap_length, e.weight
            )?;
        }
        println!("✅ GFA written → {}", path.display());
        Ok(())
    }
}

/* ------------------------------------------------------------------------- */
/*                                   TESTS                                   */
/* ------------------------------------------------------------------------- */

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core_data_structures::*;

    #[test]
    fn builder_smoke_test() {
        let builder = AssemblyGraphBuilder::new(15, 25, 2, 2).unwrap();
        let reads = vec![
            CorrectedRead {
                id: 0,
                original: "ATCGATCGATCG".into(),
                corrected: "ATCGATCGATCG".into(),
                corrections: Vec::new(),
                quality_scores: vec![30; 12],
            },
            CorrectedRead {
                id: 1,
                original: "TCGATCGATCGA".into(),
                corrected: "TCGATCGATCGA".into(),
                corrections: Vec::new(),
                quality_scores: vec![30; 12],
            },
        ];
        let g = builder.build_graph(&reads).unwrap();
        assert!(!g.graph_fragment.nodes.is_empty());
        assert!(!g.graph_fragment.edges.is_empty());
    }

    #[test]
    fn stats_and_contigs() {
        let mut asm = mocked_assembly_graph();
        asm.generate_contigs().unwrap();
        assert!(asm.assembly_stats.total_contigs > 0);
        assert!(asm.assembly_stats.n50 > 0);
    }

    #[test]
    fn simplification_reduces_complexity() {
        let builder = AssemblyGraphBuilder::new(15, 25, 1, 1).unwrap();
        let mut asm = mocked_assembly_graph();
        let pre_nodes = asm.graph_fragment.nodes.len();
        let pre_edges = asm.graph_fragment.edges.len();
        let after = builder.simplify_graph(asm).unwrap();
        assert!(after.graph_fragment.nodes.len() <= pre_nodes);
        assert!(after.graph_fragment.edges.len() <= pre_edges);
    }

    /* ----------------------------- helpers ------------------------------ */

    fn mocked_assembly_graph() -> AssemblyGraph {
        let mut frag = GraphFragment::new(0);
        for i in 0..5 {
            let seq = format!("ATCG{}", "ATCG".repeat(i));
            let kmer = CanonicalKmer::new(&seq).unwrap();
            let node = GraphNode::new(kmer, seq.len());
            frag.add_node(node);
        }
        let hashes: Vec<u64> = frag.nodes.keys().copied().collect();
        for i in 0..hashes.len() - 1 {
            let e = GraphEdge::new(hashes[i], hashes[i + 1], 3);
            frag.add_edge(e);
        }
        AssemblyGraph {
            graph_fragment: frag.clone(),
            petgraph: Graph::new(),
            contigs: Vec::new(),
            assembly_stats: AssemblyStats::default(),
        }
    }
}
