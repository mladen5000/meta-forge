use crate::utils::configuration::{AmbiguousBaseConfig, AmbiguousBaseStrategy};
use ahash::{AHashMap, RandomState};
use anyhow::{anyhow, Result};
use petgraph::{Directed, Graph};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Core data structures for the enhanced metagenomics pipeline
/// Implements proper k-mer handling, minimizer extraction, and graph structures

/// Canonical k-mer representation with efficient hashing
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CanonicalKmer {
    pub sequence: String,
    pub hash: u64,
    pub is_canonical: bool,
}

impl CanonicalKmer {
    pub fn new(kmer: &str) -> Result<Self> {
        Self::new_with_config(kmer, &AmbiguousBaseConfig::default())
    }

    pub fn new_with_config(kmer: &str, config: &AmbiguousBaseConfig) -> Result<Self> {
        // Check for valid DNA characters, allowing N for ambiguous bases
        if !kmer
            .chars()
            .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T' | 'N' | 'a' | 'c' | 'g' | 't' | 'n'))
        {
            return Err(anyhow!("Invalid DNA sequence: {}", kmer));
        }

        // Handle ambiguous bases based on configuration
        let processed_kmer = Self::handle_ambiguous_bases(kmer, config)?;

        let kmer_upper = processed_kmer.to_uppercase();
        let rc = Self::reverse_complement(&kmer_upper)?;

        let (canonical, is_canonical) = if kmer_upper <= rc {
            (kmer_upper, true)
        } else {
            (rc, false)
        };

        let hash = Self::hash_sequence(&canonical);

        Ok(Self {
            sequence: canonical,
            hash,
            is_canonical,
        })
    }

    fn handle_ambiguous_bases(kmer: &str, config: &AmbiguousBaseConfig) -> Result<String> {
        let n_count = kmer.chars().filter(|&c| matches!(c, 'N' | 'n')).count();

        match config.strategy {
            AmbiguousBaseStrategy::Skip => {
                if n_count > 0 {
                    return Err(anyhow!("K-mer contains ambiguous bases (N): {}", kmer));
                }
                Ok(kmer.to_string())
            }

            AmbiguousBaseStrategy::Allow => {
                if n_count > config.max_n_count {
                    return Err(anyhow!(
                        "K-mer contains too many ambiguous bases ({} > {}): {}",
                        n_count,
                        config.max_n_count,
                        kmer
                    ));
                }
                Ok(kmer.to_string())
            }

            AmbiguousBaseStrategy::Replace => Ok(kmer
                .chars()
                .map(|c| {
                    if matches!(c, 'N' | 'n') {
                        config.replacement_base
                    } else {
                        c
                    }
                })
                .collect()),

            AmbiguousBaseStrategy::RandomReplace => {
                let mut result = String::new();
                let probs = config
                    .random_probabilities
                    .unwrap_or([0.25, 0.25, 0.25, 0.25]);

                for c in kmer.chars() {
                    if matches!(c, 'N' | 'n') {
                        // Simple random selection based on probabilities
                        let rand_val: f64 = fastrand::f64();
                        let replacement = if rand_val < probs[0] {
                            'A'
                        } else if rand_val < probs[0] + probs[1] {
                            'C'
                        } else if rand_val < probs[0] + probs[1] + probs[2] {
                            'G'
                        } else {
                            'T'
                        };
                        result.push(replacement);
                    } else {
                        result.push(c);
                    }
                }
                Ok(result)
            }

            AmbiguousBaseStrategy::ContextReplace => {
                // Simple context replacement: use most common base in k-mer
                let mut counts = [0; 4]; // A, C, G, T
                for c in kmer.chars() {
                    match c.to_ascii_uppercase() {
                        'A' => counts[0] += 1,
                        'C' => counts[1] += 1,
                        'G' => counts[2] += 1,
                        'T' => counts[3] += 1,
                        _ => {}
                    }
                }

                let max_idx = counts
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, &count)| count)
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);

                let replacement = ['A', 'C', 'G', 'T'][max_idx];

                Ok(kmer
                    .chars()
                    .map(|c| {
                        if matches!(c, 'N' | 'n') {
                            replacement
                        } else {
                            c
                        }
                    })
                    .collect())
            }
        }
    }

    fn reverse_complement(seq: &str) -> Result<String> {
        // Manual reverse complement to handle all DNA characters including N
        Ok(seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                'N' => 'N',
                _ => c, // Should not happen with validated input
            })
            .collect())
    }

    fn hash_sequence(seq: &str) -> u64 {
        // Use a deterministic hash function for consistent results
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let mut hasher = DefaultHasher::new();
        seq.hash(&mut hasher);
        hasher.finish()
    }

    pub fn to_string(&self) -> &str {
        &self.sequence
    }

    pub fn len(&self) -> usize {
        self.sequence.len()
    }
}

/// Minimizer with position information
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Minimizer {
    pub kmer: CanonicalKmer,
    pub position: usize,
    pub window_start: usize,
    pub window_end: usize,
}

/// Rolling hash minimizer extractor
pub struct MinimizerExtractor {
    k: usize,
    w: usize,
    hasher: RandomState,
    ambiguous_config: AmbiguousBaseConfig,
}

impl MinimizerExtractor {
    pub fn new(k: usize, w: usize) -> Self {
        Self::new_with_config(k, w, AmbiguousBaseConfig::default())
    }

    pub fn new_with_config(k: usize, w: usize, ambiguous_config: AmbiguousBaseConfig) -> Self {
        Self {
            k,
            w,
            hasher: RandomState::new(),
            ambiguous_config,
        }
    }

    pub fn extract_minimizers(&self, sequence: &str) -> Result<Vec<Minimizer>> {
        if sequence.len() < self.k {
            return Ok(Vec::new());
        }

        let seq_bytes = sequence.as_bytes();
        let mut minimizers = Vec::new();

        // Extract all k-mers with their hashes
        let mut kmer_hashes = Vec::new();
        for (i, window) in seq_bytes.windows(self.k).enumerate() {
            let kmer_str = std::str::from_utf8(window)?;
            match CanonicalKmer::new_with_config(kmer_str, &self.ambiguous_config) {
                Ok(canonical_kmer) => kmer_hashes.push((i, canonical_kmer)),
                Err(_) => {
                    // Skip k-mers that couldn't be processed (e.g., too many N's)
                    continue;
                }
            }
        }

        // Find minimizers using sliding window
        let window_size = self.w.min(kmer_hashes.len());
        for i in 0..=(kmer_hashes.len().saturating_sub(window_size)) {
            let window_end = (i + window_size).min(kmer_hashes.len());
            let window = &kmer_hashes[i..window_end];

            // Find minimum hash in window
            if let Some((min_idx, min_kmer)) = window
                .iter()
                .min_by_key(|(_, kmer)| kmer.hash)
                .map(|(idx, kmer)| (*idx, kmer.clone()))
            {
                minimizers.push(Minimizer {
                    kmer: min_kmer,
                    position: min_idx,
                    window_start: i,
                    window_end: window_end - 1,
                });
            }
        }

        // Remove duplicate consecutive minimizers
        minimizers.dedup_by(|a, b| a.kmer.hash == b.kmer.hash);

        Ok(minimizers)
    }
}

/// Assembly graph node representing a k-mer
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphNode {
    pub kmer: CanonicalKmer,
    pub coverage: u32,
    pub complexity_score: f64,
    pub kmer_size: usize,
    pub read_positions: Vec<ReadPosition>,
    pub node_type: NodeType,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadPosition {
    pub read_id: usize,
    pub position: usize,
    pub strand: Strand,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum NodeType {
    Unique,
    Repetitive,
    LowCoverage,
    HighCoverage,
    Tip,
    Bubble,
}

impl GraphNode {
    pub fn new(kmer: CanonicalKmer, kmer_size: usize) -> Self {
        Self {
            kmer,
            coverage: 1,
            complexity_score: 0.0,
            kmer_size,
            read_positions: Vec::new(),
            node_type: NodeType::Unique,
        }
    }

    pub fn add_read_position(&mut self, read_id: usize, position: usize, strand: Strand) {
        self.read_positions.push(ReadPosition {
            read_id,
            position,
            strand,
        });
        self.coverage += 1;
    }

    pub fn update_node_type(&mut self) {
        self.node_type = match self.coverage {
            0..=2 => NodeType::LowCoverage,
            3..=10 => NodeType::Unique,
            11..=100 => NodeType::HighCoverage,
            _ => NodeType::Repetitive,
        };
    }

    pub fn calculate_complexity(&mut self, sequence_context: &str) {
        self.complexity_score = calculate_sequence_complexity(sequence_context);
    }
}

/// Assembly graph edge representing k-mer overlaps
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphEdge {
    pub from_hash: u64,
    pub to_hash: u64,
    pub weight: u32,
    pub overlap_length: usize,
    pub confidence: f64,
    pub edge_type: EdgeType,
    pub supporting_reads: HashSet<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeType {
    Simple,
    Complex,
    Repeat,
    LowConfidence,
}

impl GraphEdge {
    pub fn new(from_hash: u64, to_hash: u64, overlap_length: usize) -> Self {
        Self {
            from_hash,
            to_hash,
            weight: 1,
            overlap_length,
            confidence: 1.0,
            edge_type: EdgeType::Simple,
            supporting_reads: HashSet::new(),
        }
    }

    pub fn add_support(&mut self, read_id: usize) {
        self.supporting_reads.insert(read_id);
        self.weight = self.supporting_reads.len() as u32;
        self.update_confidence();
    }

    fn update_confidence(&mut self) {
        // Confidence based on number of supporting reads
        self.confidence = (self.weight as f64).ln() / 10.0;
        self.confidence = self.confidence.min(1.0).max(0.1);

        // Update edge type based on weight
        self.edge_type = match self.weight {
            1 => EdgeType::LowConfidence,
            2..=5 => EdgeType::Simple,
            6..=20 => EdgeType::Complex,
            _ => EdgeType::Repeat,
        };
    }
}

/// Complete graph fragment representing a portion of the assembly graph
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GraphFragment {
    pub nodes: AHashMap<u64, GraphNode>,
    pub edges: Vec<GraphEdge>,
    pub fragment_id: usize,
    pub coverage_stats: CoverageStats,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CoverageStats {
    pub mean_coverage: f64,
    pub median_coverage: u32,
    pub coverage_distribution: HashMap<u32, usize>,
    pub total_nodes: usize,
    pub total_edges: usize,
}

impl GraphFragment {
    /// Reconstruct sequence from a path of node hashes
    pub fn reconstruct_sequence_from_path(&self, path: &[u64]) -> Result<String> {
        if path.is_empty() {
            return Ok(String::new());
        }

        // For now, create a placeholder sequence based on path length
        // In a real implementation, this would reconstruct from k-mer overlaps
        let mut sequence = String::new();
        for (i, &hash) in path.iter().enumerate() {
            if let Some(node) = self.nodes.get(&hash) {
                if i == 0 {
                    // Add full sequence for first node
                    sequence.push_str(&node.kmer.sequence);
                } else {
                    // Add only the last character for overlap (simple k-mer model)
                    if let Some(last_char) = node.kmer.sequence.chars().last() {
                        sequence.push(last_char);
                    }
                }
            } else {
                // Fallback: use hash as sequence component
                sequence.push_str(&format!("N{i}"));
            }
        }
        Ok(sequence)
    }

    /// Calculate coverage for a path based on node hashes
    pub fn calculate_path_coverage_from_hashes(&self, path: &[u64]) -> f64 {
        if path.is_empty() {
            return 0.0;
        }

        let total_coverage: u32 = path
            .iter()
            .filter_map(|&hash| self.nodes.get(&hash))
            .map(|node| node.coverage)
            .sum();

        total_coverage as f64 / path.len() as f64
    }

    pub fn new(fragment_id: usize) -> Self {
        Self {
            nodes: AHashMap::new(),
            edges: Vec::new(),
            fragment_id,
            coverage_stats: CoverageStats::default(),
        }
    }

    pub fn add_node(&mut self, node: GraphNode) {
        let hash = node.kmer.hash;
        self.nodes.insert(hash, node);
        self.update_coverage_stats();
    }

    pub fn add_edge(&mut self, edge: GraphEdge) -> Result<()> {
        // Validate that both nodes exist before adding edge
        if !self.nodes.contains_key(&edge.from_hash) {
            return Err(anyhow!("Source node with hash {} not found in fragment", edge.from_hash));
        }
        if !self.nodes.contains_key(&edge.to_hash) {
            return Err(anyhow!("Target node with hash {} not found in fragment", edge.to_hash));
        }
        
        // Check for duplicate edges
        let edge_exists = self.edges.iter().any(|e| 
            e.from_hash == edge.from_hash && e.to_hash == edge.to_hash
        );
        
        if edge_exists {
            // Update weight instead of creating duplicate
            if let Some(existing_edge) = self.edges.iter_mut().find(|e| 
                e.from_hash == edge.from_hash && e.to_hash == edge.to_hash
            ) {
                existing_edge.weight += edge.weight;
            }
        } else {
            self.edges.push(edge);
        }
        
        self.update_coverage_stats();
        Ok(())
    }

    pub fn merge_with(&mut self, other: GraphFragment) -> Result<()> {
        // Merge nodes
        for (hash, other_node) in other.nodes {
            if let Some(existing_node) = self.nodes.get_mut(&hash) {
                // Merge existing node
                existing_node.coverage += other_node.coverage;
                existing_node
                    .read_positions
                    .extend(other_node.read_positions);
                existing_node.update_node_type();
            } else {
                // Add new node
                self.nodes.insert(hash, other_node);
            }
        }

        // Merge edges
        for other_edge in other.edges {
            // Check if edge already exists
            if let Some(existing_edge) = self
                .edges
                .iter_mut()
                .find(|e| e.from_hash == other_edge.from_hash && e.to_hash == other_edge.to_hash)
            {
                existing_edge.weight += other_edge.weight;
                existing_edge
                    .supporting_reads
                    .extend(other_edge.supporting_reads);
                existing_edge.update_confidence();
            } else {
                self.edges.push(other_edge);
            }
        }

        self.update_coverage_stats();
        Ok(())
    }

    // Duplicate function definitions removed
    fn update_coverage_stats(&mut self) {
        let coverages: Vec<u32> = self.nodes.values().map(|n| n.coverage).collect();

        if !coverages.is_empty() {
            self.coverage_stats.mean_coverage =
                coverages.iter().sum::<u32>() as f64 / coverages.len() as f64;

            let mut sorted_coverages = coverages.clone();
            sorted_coverages.sort_unstable();
            self.coverage_stats.median_coverage = sorted_coverages[sorted_coverages.len() / 2];

            // Build coverage distribution
            let mut distribution = HashMap::new();
            for &coverage in &coverages {
                *distribution.entry(coverage).or_insert(0) += 1;
            }
            self.coverage_stats.coverage_distribution = distribution;
        }

        self.coverage_stats.total_nodes = self.nodes.len();
        self.coverage_stats.total_edges = self.edges.len();
    }

    pub fn get_adjacency_list(&self) -> AHashMap<u64, Vec<u64>> {
        let mut adjacency = AHashMap::new();

        for edge in &self.edges {
            adjacency
                .entry(edge.from_hash)
                .or_insert_with(Vec::new)
                .push(edge.to_hash);
        }

        adjacency
    }

    pub fn find_tips(&self) -> Vec<u64> {
        let adjacency = self.get_adjacency_list();
        let mut in_degree = AHashMap::new();
        let mut out_degree = AHashMap::new();

        // Calculate degrees
        for (&from, neighbors) in &adjacency {
            out_degree.insert(from, neighbors.len());
            for &to in neighbors {
                *in_degree.entry(to).or_insert(0) += 1;
            }
        }

        // Find nodes with in-degree or out-degree of 0
        let mut tips = Vec::new();
        for &node_hash in self.nodes.keys() {
            let in_deg = in_degree.get(&node_hash).copied().unwrap_or(0);
            let out_deg = out_degree.get(&node_hash).copied().unwrap_or(0);

            if in_deg == 0 || out_deg == 0 {
                tips.push(node_hash);
            }
        }

        tips
    }

    pub fn find_bubbles(&self) -> Vec<BubbleStructure> {
        let adjacency = self.get_adjacency_list();
        let mut bubbles = Vec::new();

        // Simple bubble detection: nodes with multiple outgoing edges that reconverge
        for (&start_node, neighbors) in &adjacency {
            if neighbors.len() > 1 {
                // Check if paths reconverge
                let mut paths = Vec::new();
                for &neighbor in neighbors {
                    if let Some(path) = self.trace_path(neighbor, &adjacency, 5) {
                        paths.push(path);
                    }
                }

                // Check for reconvergence
                if let Some(bubble) = self.detect_reconvergence(&paths, start_node) {
                    bubbles.push(bubble);
                }
            }
        }

        bubbles
    }

    fn trace_path(
        &self,
        start: u64,
        adjacency: &AHashMap<u64, Vec<u64>>,
        max_length: usize,
    ) -> Option<Vec<u64>> {
        let mut path = vec![start];
        let mut current = start;

        for _ in 0..max_length {
            if let Some(neighbors) = adjacency.get(&current) {
                if neighbors.len() == 1 {
                    current = neighbors[0];
                    path.push(current);
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        if path.len() > 1 {
            Some(path)
        } else {
            None
        }
    }

    fn detect_reconvergence(&self, paths: &[Vec<u64>], start_node: u64) -> Option<BubbleStructure> {
        if paths.len() < 2 {
            return None;
        }

        // Find common end nodes
        let end_nodes: HashSet<u64> = paths
            .iter()
            .filter_map(|path| path.last().copied())
            .collect();

        if end_nodes.len() == 1 {
            let end_node = *end_nodes.iter().next().unwrap();
            Some(BubbleStructure {
                start_node,
                end_node,
                alternative_paths: paths.to_vec(),
                bubble_type: if paths.len() == 2 {
                    BubbleType::Simple
                } else {
                    BubbleType::Complex
                },
            })
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BubbleStructure {
    pub start_node: u64,
    pub end_node: u64,
    pub alternative_paths: Vec<Vec<u64>>,
    pub bubble_type: BubbleType,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum BubbleType {
    Simple,
    Complex,
    Nested,
}

/// Types for assembly graph integration
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ContigType {
    Linear,
    Circular,
    Scaffold,
    Unitig,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum EulerianPathType {
    Circuit,
    Path,
    Multiple,
    None,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeWeight {
    pub weight: u32,
    pub confidence: f64,
}

// Duplicate ContigType removed - using the existing one above

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Contig {
    pub id: usize,
    pub sequence: String,
    pub coverage: f64,
    pub length: usize,
    pub node_path: Vec<u64>,
    pub contig_type: ContigType,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct AssemblyStats {
    pub total_length: usize,
    pub num_contigs: usize,
    pub n50: usize,
    pub n90: usize,
    pub largest_contig: usize,
    pub gc_content: f64,
    pub coverage_mean: f64,
    pub coverage_std: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrectionMetadata {
    pub algorithm: String,
    pub confidence_threshold: f64,
    pub context_window: usize,
    pub correction_time_ms: u64,
}

/// Complete assembly graph with petgraph integration
#[derive(Debug, Clone)]
pub struct AssemblyGraph {
    pub graph_fragment: GraphFragment,
    pub petgraph: Graph<u64, EdgeWeight, Directed>,
    pub contigs: Vec<Contig>,
    pub assembly_stats: AssemblyStats,
}

/// Assembly chunk representing a batch of reads and their graph contribution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyChunk {
    pub chunk_id: usize,
    pub reads: Vec<CorrectedRead>,
    pub k_size: usize,
    pub graph_fragment: GraphFragment,
    pub processing_stats: ProcessingStats,
}

impl AssemblyChunk {
    // Main implementation is below - removed duplicate
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrectedRead {
    pub id: usize,
    pub original: String,
    pub corrected: String,
    pub corrections: Vec<BaseCorrection>,
    pub quality_scores: Vec<u8>,
    pub correction_metadata: CorrectionMetadata,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BaseCorrection {
    pub position: usize,
    pub from: char,
    pub to: char,
    pub confidence: f64,
    pub correction_type: CorrectionType,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum CorrectionType {
    Substitution,
    Insertion,
    Deletion,
    QualityImprovement,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ProcessingStats {
    pub reads_processed: usize,
    pub kmers_extracted: usize,
    pub minimizers_found: usize,
    pub nodes_created: usize,
    pub edges_created: usize,
    pub processing_time_ms: u64,
    pub memory_usage_bytes: usize,
}

impl AssemblyChunk {
    pub fn new(chunk_id: usize, k_size: usize) -> Self {
        Self {
            chunk_id,
            reads: Vec::new(),
            k_size,
            graph_fragment: GraphFragment::new(chunk_id),
            processing_stats: ProcessingStats::default(),
        }
    }

    pub fn add_read(&mut self, read: CorrectedRead) -> Result<()> {
        let start_time = std::time::Instant::now();

        // Simple k-mer extraction instead of minimizers for more reliable results
        let sequence = &read.corrected;
        if sequence.len() < self.k_size {
            println!(
                "   Warning: Read {} too short ({} bp) for k-mer size {}",
                read.id,
                sequence.len(),
                self.k_size
            );
            return Ok(());
        }

        let mut kmers = Vec::new();
        let sequence_bytes = sequence.as_bytes();
        
        // Pre-allocate capacity for better performance
        kmers.reserve_exact(sequence_bytes.len().saturating_sub(self.k_size - 1));
        
        // Optimized k-mer extraction with fewer allocations
        for (i, window) in sequence_bytes.windows(self.k_size).enumerate() {
            // Quick byte validation before UTF-8 conversion
            let mut valid_nucleotides = true;
            for &byte in window {
                match byte {
                    b'A' | b'T' | b'G' | b'C' | b'a' | b't' | b'g' | b'c' | b'N' | b'n' |
                    b'Y' | b'y' | b'R' | b'r' | b'W' | b'w' | b'S' | b's' | b'K' | b'k' |
                    b'M' | b'm' | b'D' | b'd' | b'V' | b'v' | b'H' | b'h' | b'B' | b'b' => {},
                    _ => { valid_nucleotides = false; break; }
                }
            }
            
            if !valid_nucleotides {
                continue;
            }

            // Convert to string only after validation
            let kmer_str = match std::str::from_utf8(window) {
                Ok(s) => s,
                Err(_) => continue, // Skip invalid UTF-8, already logged elsewhere
            };

            match CanonicalKmer::new(kmer_str) {
                Ok(canonical_kmer) => {
                    kmers.push((i, canonical_kmer));
                }
                Err(_) => continue, // Skip invalid k-mers, reduced logging for performance
            }
        }

        println!(
            "   Read {}: extracted {} k-mers from {} bp sequence",
            read.id,
            kmers.len(),
            sequence.len()
        );
        self.processing_stats.minimizers_found += kmers.len();

        // Create nodes from k-mers
        for (pos, kmer) in &kmers {
            self.add_or_update_node(kmer, read.id, *pos)?;
        }

        // Create edges between consecutive k-mers
        for window in kmers.windows(2) {
            let (_, kmer1) = &window[0];
            let (_, kmer2) = &window[1];

            let edge = GraphEdge::new(
                kmer1.hash,
                kmer2.hash,
                self.k_size - 1, // Overlap length
            );
            self.graph_fragment.add_edge(edge);
            self.processing_stats.edges_created += 1;
        }

        self.reads.push(read);
        self.processing_stats.reads_processed += 1;
        self.processing_stats.processing_time_ms += start_time.elapsed().as_millis() as u64;

        Ok(())
    }

    fn add_or_update_node(
        &mut self,
        kmer: &CanonicalKmer,
        read_id: usize,
        position: usize,
    ) -> Result<()> {
        if let Some(existing_node) = self.graph_fragment.nodes.get_mut(&kmer.hash) {
            existing_node.add_read_position(read_id, position, Strand::Forward);
        } else {
            let mut new_node = GraphNode::new(kmer.clone(), self.k_size);
            new_node.add_read_position(read_id, position, Strand::Forward);
            self.graph_fragment.add_node(new_node);
            self.processing_stats.nodes_created += 1;
        }

        Ok(())
    }

    pub fn finalize(&mut self) {
        // Update node types and complexity scores
        for node in self.graph_fragment.nodes.values_mut() {
            node.update_node_type();
            // Calculate complexity based on coverage and read distribution
            let kmer_seq = node.kmer.sequence.clone();
            node.calculate_complexity(&kmer_seq);
        }

        // Update fragment statistics
        self.graph_fragment.update_coverage_stats();

        // Calculate memory usage
        self.processing_stats.memory_usage_bytes = self.estimate_memory_usage();
    }

    pub fn estimate_memory_usage(&self) -> usize {
        let nodes_size = self.graph_fragment.nodes.len() * std::mem::size_of::<GraphNode>();
        let edges_size = self.graph_fragment.edges.len() * std::mem::size_of::<GraphEdge>();
        let reads_size = self.reads.iter().map(|r| r.corrected.len()).sum::<usize>();

        nodes_size + edges_size + reads_size
    }
}

/// Graph update structure for incremental graph building
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GraphUpdate {
    pub nodes_to_add: Vec<(u64, GraphNode)>,
    pub edges_to_add: Vec<GraphEdge>,
    pub nodes_to_update: Vec<(u64, NodeUpdate)>,
    pub edges_to_update: Vec<(u64, u64, EdgeUpdate)>,
    pub read_id: usize,
    pub timestamp: std::time::SystemTime,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeUpdate {
    pub coverage_increment: u32,
    pub new_read_positions: Vec<ReadPosition>,
    pub complexity_update: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeUpdate {
    pub weight_increment: u32,
    pub new_supporting_reads: HashSet<usize>,
    pub confidence_update: Option<f64>,
}

impl GraphUpdate {
    pub fn new(read_id: usize) -> Self {
        Self {
            nodes_to_add: Vec::new(),
            edges_to_add: Vec::new(),
            nodes_to_update: Vec::new(),
            edges_to_update: Vec::new(),
            read_id,
            timestamp: std::time::SystemTime::now(),
        }
    }

    pub fn add_node(&mut self, hash: u64, node: GraphNode) {
        self.nodes_to_add.push((hash, node));
    }

    pub fn add_edge(&mut self, edge: GraphEdge) -> Result<()> {
        self.edges_to_add.push(edge);
        Ok(())
    }

    pub fn update_node(&mut self, hash: u64, update: NodeUpdate) {
        self.nodes_to_update.push((hash, update));
    }

    pub fn update_edge(&mut self, from_hash: u64, to_hash: u64, update: EdgeUpdate) {
        self.edges_to_update.push((from_hash, to_hash, update));
    }
}

/// Utility functions for sequence analysis
pub fn calculate_sequence_complexity(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }

    // Calculate Shannon entropy
    let mut counts = [0u32; 4]; // A, C, G, T
    for byte in sequence.bytes() {
        match byte {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {}
        }
    }

    let total = counts.iter().sum::<u32>() as f64;
    if total == 0.0 {
        return 0.0;
    }

    let entropy = counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / total;
            -p * p.log2()
        })
        .sum::<f64>();

    // Normalize entropy (max entropy for DNA is 2.0)
    entropy / 2.0
}

pub fn calculate_gc_content(sequence: &str) -> f64 {
    if sequence.is_empty() {
        return 0.0;
    }

    let gc_count = sequence
        .chars()
        .filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c')
        .count();

    gc_count as f64 / sequence.len() as f64
}

pub fn validate_dna_sequence(sequence: &str) -> Result<()> {
    for (i, c) in sequence.chars().enumerate() {
        if !matches!(c, 'A' | 'C' | 'G' | 'T' | 'a' | 'c' | 'g' | 't' | 'N' | 'n') {
            return Err(anyhow!("Invalid DNA character '{}' at position {}", c, i));
        }
    }
    Ok(())
}

impl Default for AssemblyGraph {
    fn default() -> Self {
        Self::new()
    }
}

impl AssemblyGraph {
    pub fn new() -> Self {
        Self {
            graph_fragment: GraphFragment::new(0),
            petgraph: petgraph::Graph::new(),
            contigs: Vec::new(),
            assembly_stats: AssemblyStats::default(),
        }
    }

    pub fn from_fragment(fragment: GraphFragment) -> Self {
        Self {
            graph_fragment: fragment,
            petgraph: petgraph::Graph::new(),
            contigs: Vec::new(),
            assembly_stats: AssemblyStats::default(),
        }
    }

    pub fn calculate_assembly_stats(&mut self) {
        if self.contigs.is_empty() {
            return;
        }

        // Sort contigs by length
        let mut lengths: Vec<usize> = self.contigs.iter().map(|c| c.length).collect();
        lengths.sort_by(|a, b| b.cmp(a));

        self.assembly_stats.total_length = lengths.iter().sum();
        self.assembly_stats.num_contigs = lengths.len();
        self.assembly_stats.largest_contig = lengths[0];

        // Calculate N50 and N90
        let half_length = self.assembly_stats.total_length / 2;
        let ninety_percent = (self.assembly_stats.total_length as f64 * 0.1) as usize;

        let mut cumulative = 0;
        for &length in &lengths {
            cumulative += length;
            if self.assembly_stats.n50 == 0 && cumulative >= half_length {
                self.assembly_stats.n50 = length;
            }
            if self.assembly_stats.n90 == 0 && cumulative >= ninety_percent {
                self.assembly_stats.n90 = length;
            }
        }

        // Calculate average coverage
        let total_coverage: f64 = self
            .contigs
            .iter()
            .map(|c| c.coverage * c.length as f64)
            .sum();
        self.assembly_stats.coverage_mean =
            total_coverage / self.assembly_stats.total_length as f64;

        // Calculate GC content
        let gc_count: usize = self
            .contigs
            .iter()
            .map(|c| {
                c.sequence
                    .chars()
                    .filter(|&ch| ch == 'G' || ch == 'C')
                    .count()
            })
            .sum();
        self.assembly_stats.gc_content = gc_count as f64 / self.assembly_stats.total_length as f64;
    }

    // Duplicate functions removed - using GraphFragment implementations

    /// Write contigs to FASTA format
    pub fn write_contigs_fasta<P: AsRef<std::path::Path>>(&self, path: P) -> Result<()> {
        use std::io::Write;
        let mut file = std::fs::File::create(path.as_ref())?;

        for (i, contig) in self.contigs.iter().enumerate() {
            writeln!(
                file,
                ">contig_{} length={} coverage={:.2}",
                i + 1,
                contig.length,
                contig.coverage
            )?;

            // Write sequence in lines of 80 characters
            for line in contig.sequence.as_bytes().chunks(80) {
                writeln!(file, "{}", std::str::from_utf8(line)?)?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_canonical_kmer() {
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        assert_eq!(kmer.sequence, "ATCG"); // ATCG < CGAT (reverse complement)

        let kmer2 = CanonicalKmer::new("CGAT").unwrap();
        assert_eq!(kmer2.sequence, "ATCG"); // Should be canonicalized to ATCG

        assert_eq!(kmer.hash, kmer2.hash); // Same hash for canonical form
    }

    #[test]
    fn test_minimizer_extraction() {
        let extractor = MinimizerExtractor::new(3, 5);
        let sequence = "ATCGATCG";
        let minimizers = extractor.extract_minimizers(sequence).unwrap();

        assert!(!minimizers.is_empty());

        // Check that minimizers are within bounds
        for min in &minimizers {
            assert!(min.position < sequence.len() - 3 + 1);
            assert_eq!(min.kmer.len(), 3);
        }
    }

    #[test]
    fn test_graph_fragment_operations() {
        let mut fragment = GraphFragment::new(0);

        // Add a test node
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let node = GraphNode::new(kmer.clone(), 4);
        fragment.add_node(node);

        assert_eq!(fragment.nodes.len(), 1);
        assert!(fragment.nodes.contains_key(&kmer.hash));

        // Add an edge
        let kmer2 = CanonicalKmer::new("TCGA").unwrap();
        let node2 = GraphNode::new(kmer2.clone(), 4);
        fragment.add_node(node2);

        let edge = GraphEdge::new(kmer.hash, kmer2.hash, 3);
        fragment.add_edge(edge);

        assert_eq!(fragment.edges.len(), 1);
        assert_eq!(fragment.nodes.len(), 2);
    }

    #[test]
    fn test_assembly_chunk() {
        let mut chunk = AssemblyChunk::new(0, 4);

        let read = CorrectedRead {
            id: 0,
            original: "ATCGATCG".to_string(),
            corrected: "ATCGATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; 8],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.8,
                context_window: 3,
                correction_time_ms: 0,
            },
        };

        chunk.add_read(read).unwrap();
        chunk.finalize();

        assert_eq!(chunk.processing_stats.reads_processed, 1);
        assert!(chunk.processing_stats.nodes_created > 0);
        assert!(!chunk.graph_fragment.nodes.is_empty());
    }

    #[test]
    fn test_sequence_complexity() {
        // Uniform sequence (low complexity)
        let complexity1 = calculate_sequence_complexity("AAAAAAAAAA");
        assert!(complexity1 < 0.1);

        // Random sequence (high complexity)
        let complexity2 = calculate_sequence_complexity("ATCGATCGAT");
        assert!(complexity2 > 0.8);

        // Mixed sequence
        let complexity3 = calculate_sequence_complexity("AAAATCGATC");
        assert!(complexity3 > complexity1 && complexity3 < complexity2);
    }

    #[test]
    fn test_gc_content() {
        assert_eq!(calculate_gc_content("ATCG"), 0.5);
        assert_eq!(calculate_gc_content("AAAA"), 0.0);
        assert_eq!(calculate_gc_content("CCGG"), 1.0);
        assert_eq!(calculate_gc_content(""), 0.0);
    }

    #[test]
    fn test_graph_merge() {
        let mut fragment1 = GraphFragment::new(0);
        let mut fragment2 = GraphFragment::new(1);

        // Add overlapping content
        let kmer = CanonicalKmer::new("ATCG").unwrap();
        let node1 = GraphNode::new(kmer.clone(), 4);
        let node2 = GraphNode::new(kmer.clone(), 4);

        fragment1.add_node(node1);
        fragment2.add_node(node2);

        fragment1.merge_with(fragment2).unwrap();

        // Should have merged the nodes
        assert_eq!(fragment1.nodes.len(), 1);
        let merged_node = fragment1.nodes.get(&kmer.hash).unwrap();
        assert_eq!(merged_node.coverage, 2); // Combined coverage
    }
}
