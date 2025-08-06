use std::collections::{HashMap, HashSet, VecDeque, BinaryHeap};
use ahash::{AHashMap, AHashSet};
use anyhow::{Result, anyhow};
use petgraph::{Graph, Directed, NodeIndex, EdgeIndex};
use petgraph::algo::{connected_components, tarjan_scc, dijkstra};
use rayon::prelude::*;
use serde::{Serialize, Deserialize};

use crate::core_data_structures::*;

/// Complete assembly graph builder with adaptive k-mer selection
pub struct AssemblyGraphBuilder {
    /// Base k-mer size
    base_k: usize,
    /// Maximum k-mer size for adaptive selection
    max_k: usize,
    /// Minimum coverage threshold for including k-mers
    min_coverage: u32,
    /// Complexity threshold for switching k-mer sizes
    complexity_threshold: f64,
    /// Thread pool for parallel processing
    thread_pool: rayon::ThreadPool,
}

impl AssemblyGraphBuilder {
    pub fn new(base_k: usize, max_k: usize, min_coverage: u32, num_threads: usize) -> Result<Self> {
        let thread_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .thread_name(|i| format!("assembly-{}", i))
            .build()?;
        
        Ok(Self {
            base_k,
            max_k,
            min_coverage,
            complexity_threshold: 0.7,
            thread_pool,
        })
    }
    
    /// Build assembly graph from corrected reads
    pub fn build_graph(&self, reads: &[CorrectedRead]) -> Result<AssemblyGraph> {
        println!("🧬 Building assembly graph from {} reads", reads.len());
        
        // Phase 1: Parallel k-mer extraction and chunking
        let chunks = self.thread_pool.install(|| {
            self.create_assembly_chunks(reads)
        })?;
        
        // Phase 2: Build graph fragments for each chunk
        let fragments: Result<Vec<_>> = chunks.par_iter()
            .map(|chunk| self.build_fragment(chunk))
            .collect();
        let fragments = fragments?;
        
        // Phase 3: Merge fragments into final graph
        let merged_graph = self.merge_fragments(fragments)?;
        
        // Phase 4: Apply graph simplification
        let simplified_graph = self.simplify_graph(merged_graph)?;
        
        Ok(simplified_graph)
    }
    
    /// Create assembly chunks with adaptive k-mer selection
    fn create_assembly_chunks(&self, reads: &[CorrectedRead]) -> Result<Vec<AssemblyChunk>> {
        let chunk_size = 1000; // Process 1000 reads per chunk
        let mut chunks = Vec::new();
        
        for (chunk_id, read_batch) in reads.chunks(chunk_size).enumerate() {
            // Analyze batch complexity to determine optimal k
            let batch_complexity = self.analyze_batch_complexity(read_batch);
            let optimal_k = self.select_optimal_k(batch_complexity);
            
            let mut chunk = AssemblyChunk::new(chunk_id, optimal_k);
            
            // Process reads in this chunk
            for read in read_batch {
                chunk.add_read(read.clone())?;
            }
            
            chunk.finalize();
            chunks.push(chunk);
        }
        
        Ok(chunks)
    }
    
    fn analyze_batch_complexity(&self, reads: &[CorrectedRead]) -> f64 {
        let total_complexity: f64 = reads.par_iter()
            .map(|read| calculate_sequence_complexity(&read.corrected))
            .sum();
        
        total_complexity / reads.len() as f64
    }
    
    fn select_optimal_k(&self, complexity: f64) -> usize {
        // Higher complexity → larger k-mer size
        let k_range = self.max_k - self.base_k;
        let k_adjustment = (complexity * k_range as f64) as usize;
        
        (self.base_k + k_adjustment).min(self.max_k)
    }
    
    /// Build graph fragment from assembly chunk
    fn build_fragment(&self, chunk: &AssemblyChunk) -> Result<GraphFragment> {
        let mut fragment = chunk.graph_fragment.clone();
        
        // Filter low-coverage nodes
        self.filter_low_coverage_nodes(&mut fragment)?;
        
        // Detect and mark structural features
        self.detect_structural_features(&mut fragment)?;
        
        // Calculate edge weights and confidences
        self.calculate_edge_weights(&mut fragment)?;
        
        Ok(fragment)
    }
    
    fn filter_low_coverage_nodes(&self, fragment: &mut GraphFragment) -> Result<()> {
        let mut nodes_to_remove = Vec::new();
        
        for (&hash, node) in &fragment.nodes {
            if node.coverage < self.min_coverage {
                nodes_to_remove.push(hash);
            }
        }
        
        // Remove low-coverage nodes and associated edges
        for hash in nodes_to_remove {
            fragment.nodes.remove(&hash);
            fragment.edges.retain(|edge| edge.from_hash != hash && edge.to_hash != hash);
        }
        
        println!("   Filtered {} low-coverage nodes", fragment.nodes.len());
        Ok(())
    }
    
    fn detect_structural_features(&self, fragment: &mut GraphFragment) -> Result<()> {
        // Find tips (dead ends)
        let tips = fragment.find_tips();
        for &tip_hash in &tips {
            if let Some(node) = fragment.nodes.get_mut(&tip_hash) {
                node.node_type = NodeType::Tip;
            }
        }
        
        // Find bubbles (alternative paths)
        let bubbles = fragment.find_bubbles();
        for bubble in &bubbles {
            // Mark bubble nodes
            for path in &bubble.alternative_paths {
                for &node_hash in path {
                    if let Some(node) = fragment.nodes.get_mut(&node_hash) {
                        node.node_type = NodeType::Bubble;
                    }
                }
            }
        }
        
        println!("   Detected {} tips and {} bubbles", tips.len(), bubbles.len());
        Ok(())
    }
    
    fn calculate_edge_weights(&self, fragment: &mut GraphFragment) -> Result<()> {
        for edge in &mut fragment.edges {
            // Weight based on supporting reads and coverage
            let from_coverage = fragment.nodes.get(&edge.from_hash)
                .map(|n| n.coverage)
                .unwrap_or(1);
            let to_coverage = fragment.nodes.get(&edge.to_hash)
                .map(|n| n.coverage)
                .unwrap_or(1);
            
            // Confidence based on minimum coverage and number of supporting reads
            let min_coverage = from_coverage.min(to_coverage);
            let support_score = edge.supporting_reads.len() as f64;
            
            edge.confidence = (support_score * min_coverage as f64).sqrt() / 100.0;
            edge.confidence = edge.confidence.min(1.0).max(0.1);
        }
        
        Ok(())
    }
    
    /// Merge multiple graph fragments into a single assembly graph
    fn merge_fragments(&self, fragments: Vec<GraphFragment>) -> Result<AssemblyGraph> {
        println!("🔗 Merging {} graph fragments", fragments.len());
        
        let mut merged_fragment = GraphFragment::new(0);
        
        // Merge all fragments
        for fragment in fragments {
            merged_fragment.merge_with(fragment)?;
        }
        
        // Convert to petgraph for advanced algorithms
        let petgraph = self.convert_to_petgraph(&merged_fragment)?;
        
        // Create final assembly graph
        let assembly_graph = AssemblyGraph {
            graph_fragment: merged_fragment,
            petgraph,
            contigs: Vec::new(),
            assembly_stats: AssemblyStats::default(),
        };
        
        Ok(assembly_graph)
    }
    
    fn convert_to_petgraph(&self, fragment: &GraphFragment) -> Result<Graph<u64, EdgeWeight, Directed>> {
        let mut graph = Graph::new();
        let mut node_map = AHashMap::new();
        
        // Add nodes
        for &node_hash in fragment.nodes.keys() {
            let node_idx = graph.add_node(node_hash);
            node_map.insert(node_hash, node_idx);
        }
        
        // Add edges
        for edge in &fragment.edges {
            if let (Some(&from_idx), Some(&to_idx)) = 
                (node_map.get(&edge.from_hash), node_map.get(&edge.to_hash)) 
            {
                let edge_weight = EdgeWeight {
                    weight: edge.weight,
                    confidence: edge.confidence,
                };
                graph.add_edge(from_idx, to_idx, edge_weight);
            }
        }
        
        println!("   Created petgraph with {} nodes and {} edges", 
            graph.node_count(), graph.edge_count());
        
        Ok(graph)
    }
    
    /// Simplify graph by removing redundant structures
    fn simplify_graph(&self, mut assembly_graph: AssemblyGraph) -> Result<AssemblyGraph> {
        println!("✂️  Simplifying assembly graph");
        
        // Remove tips (short dead-end branches)
        self.remove_tips(&mut assembly_graph)?;
        
        // Pop bubbles (remove alternative paths in small bubbles)
        self.pop_bubbles(&mut assembly_graph)?;
        
        // Merge linear chains
        self.merge_linear_chains(&mut assembly_graph)?;
        
        // Remove low-confidence edges
        self.remove_low_confidence_edges(&mut assembly_graph)?;
        
        Ok(assembly_graph)
    }
    
    fn remove_tips(&self, assembly_graph: &mut AssemblyGraph) -> Result<()> {
        let max_tip_length = 100; // Maximum length of tips to remove
        let min_coverage_ratio = 0.1; // Tips must have very low coverage relative to neighbors
        
        let tips = assembly_graph.graph_fragment.find_tips();
        let mut removed_count = 0;
        
        for &tip_hash in &tips {
            if let Some(tip_node) = assembly_graph.graph_fragment.nodes.get(&tip_hash) {
                // Check if tip should be removed
                if tip_node.kmer.len() <= max_tip_length {
                    // Check coverage relative to neighbors
                    let neighbor_coverage = self.get_neighbor_coverage(&assembly_graph.graph_fragment, tip_hash);
                    let coverage_ratio = tip_node.coverage as f64 / neighbor_coverage.max(1.0);
                    
                    if coverage_ratio < min_coverage_ratio {
                        // Remove tip
                        assembly_graph.graph_fragment.nodes.remove(&tip_hash);
                        assembly_graph.graph_fragment.edges.retain(|edge| 
                            edge.from_hash != tip_hash && edge.to_hash != tip_hash
                        );
                        removed_count += 1;
                    }
                }
            }
        }
        
        println!("   Removed {} tips", removed_count);
        Ok(())
    }
    
    fn get_neighbor_coverage(&self, fragment: &GraphFragment, node_hash: u64) -> f64 {
        let mut neighbor_coverages = Vec::new();
        
        for edge in &fragment.edges {
            if edge.from_hash == node_hash {
                if let Some(neighbor) = fragment.nodes.get(&edge.to_hash) {
                    neighbor_coverages.push(neighbor.coverage as f64);
                }
            } else if edge.to_hash == node_hash {
                if let Some(neighbor) = fragment.nodes.get(&edge.from_hash) {
                    neighbor_coverages.push(neighbor.coverage as f64);
                }
            }
        }
        
        if neighbor_coverages.is_empty() {
            0.0
        } else {
            neighbor_coverages.iter().sum::<f64>() / neighbor_coverages.len() as f64
        }
    }
    
    fn pop_bubbles(&self, assembly_graph: &mut AssemblyGraph) -> Result<()> {
        let bubbles = assembly_graph.graph_fragment.find_bubbles();
        let mut popped_count = 0;
        
        for bubble in bubbles {
            if bubble.bubble_type == BubbleType::Simple && bubble.alternative_paths.len() == 2 {
                // Choose the path with higher coverage
                let path1_coverage = self.calculate_path_coverage(&assembly_graph.graph_fragment, &bubble.alternative_paths[0]);
                let path2_coverage = self.calculate_path_coverage(&assembly_graph.graph_fragment, &bubble.alternative_paths[1]);
                
                let path_to_remove = if path1_coverage > path2_coverage {
                    &bubble.alternative_paths[1]
                } else {
                    &bubble.alternative_paths[0]
                };
                
                // Remove the lower-coverage path
                for &node_hash in path_to_remove {
                    assembly_graph.graph_fragment.nodes.remove(&node_hash);
                }
                
                // Remove associated edges
                assembly_graph.graph_fragment.edges.retain(|edge| {
                    !path_to_remove.contains(&edge.from_hash) && !path_to_remove.contains(&edge.to_hash)
                });
                
                popped_count += 1;
            }
        }
        
        println!("   Popped {} bubbles", popped_count);
        Ok(())
    }
    
    fn calculate_path_coverage(&self, fragment: &GraphFragment, path: &[u64]) -> f64 {
        let total_coverage: u32 = path.iter()
            .filter_map(|&hash| fragment.nodes.get(&hash))
            .map(|node| node.coverage)
            .sum();
        
        if path.is_empty() {
            0.0
        } else {
            total_coverage as f64 / path.len() as f64
        }
    }
    
    fn merge_linear_chains(&self, assembly_graph: &mut AssemblyGraph) -> Result<()> {
        // Find linear chains (nodes with exactly one incoming and one outgoing edge)
        let mut linear_chains = self.find_linear_chains(&assembly_graph.graph_fragment);
        let mut merged_count = 0;
        
        for chain in linear_chains {
            if chain.len() > 1 {
                // Merge chain into a single node
                let merged_node = self.merge_chain_nodes(&assembly_graph.graph_fragment, &chain)?;
                
                // Remove original nodes
                for &node_hash in &chain[1..] {
                    assembly_graph.graph_fragment.nodes.remove(&node_hash);
                }
                
                // Update first node with merged content
                if let Some(first_node) = assembly_graph.graph_fragment.nodes.get_mut(&chain[0]) {
                    *first_node = merged_node;
                }
                
                // Remove internal edges
                assembly_graph.graph_fragment.edges.retain(|edge| {
                    !chain[1..].contains(&edge.from_hash) && !chain[1..].contains(&edge.to_hash)
                });
                
                merged_count += 1;
            }
        }
        
        println!("   Merged {} linear chains", merged_count);
        Ok(())
    }
    
    fn find_linear_chains(&self, fragment: &GraphFragment) -> Vec<Vec<u64>> {
        let adjacency = fragment.get_adjacency_list();
        let mut in_degree = AHashMap::new();
        let mut visited = AHashSet::new();
        let mut chains = Vec::new();
        
        // Calculate in-degrees
        for neighbors in adjacency.values() {
            for &neighbor in neighbors {
                *in_degree.entry(neighbor).or_insert(0) += 1;
            }
        }
        
        // Find chain starts (nodes with in-degree != 1 or out-degree != 1)
        for (&node, neighbors) in &adjacency {
            if visited.contains(&node) {
                continue;
            }
            
            let in_deg = in_degree.get(&node).copied().unwrap_or(0);
            let out_deg = neighbors.len();
            
            // Start a new chain if this is not a simple linear node
            if in_deg != 1 || out_deg != 1 {
                let chain = self.trace_linear_chain(node, &adjacency, &in_degree, &mut visited);
                if !chain.is_empty() {
                    chains.push(chain);
                }
            }
        }
        
        chains
    }
    
    fn trace_linear_chain(
        &self,
        start: u64,
        adjacency: &AHashMap<u64, Vec<u64>>,
        in_degree: &AHashMap<u64, usize>,
        visited: &mut AHashSet<u64>,
    ) -> Vec<u64> {
        let mut chain = vec![start];
        let mut current = start;
        visited.insert(start);
        
        // Follow the linear path
        while let Some(neighbors) = adjacency.get(&current) {
            if neighbors.len() == 1 {
                let next = neighbors[0];
                let next_in_degree = in_degree.get(&next).copied().unwrap_or(0);
                
                if next_in_degree == 1 && !visited.contains(&next) {
                    chain.push(next);
                    visited.insert(next);
                    current = next;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        
        chain
    }
    
    fn merge_chain_nodes(&self, fragment: &GraphFragment, chain: &[u64]) -> Result<GraphNode> {
        if chain.is_empty() {
            return Err(anyhow!("Cannot merge empty chain"));
        }
        
        let first_node = fragment.nodes.get(&chain[0])
            .ok_or_else(|| anyhow!("First node in chain not found"))?;
        
        let mut merged_sequence = first_node.kmer.sequence.clone();
        let mut total_coverage = first_node.coverage;
        let mut all_read_positions = first_node.read_positions.clone();
        
        // Merge subsequent nodes
        for &node_hash in &chain[1..] {
            if let Some(node) = fragment.nodes.get(&node_hash) {
                // Append k-mer sequence (overlapping by k-1)
                let overlap = first_node.kmer_size - 1;
                if node.kmer.sequence.len() > overlap {
                    merged_sequence.push_str(&node.kmer.sequence[overlap..]);
                }
                
                total_coverage += node.coverage;
                all_read_positions.extend(node.read_positions.iter().cloned());
            }
        }
        
        // Create merged canonical k-mer
        let merged_kmer = CanonicalKmer::new(&merged_sequence)?;
        
        let mut merged_node = GraphNode::new(merged_kmer, merged_sequence.len());
        merged_node.coverage = total_coverage;
        merged_node.read_positions = all_read_positions;
        merged_node.complexity_score = calculate_sequence_complexity(&merged_sequence);
        merged_node.update_node_type();
        
        Ok(merged_node)
    }
    
    fn remove_low_confidence_edges(&self, assembly_graph: &mut AssemblyGraph) -> Result<()> {
        let confidence_threshold = 0.3;
        let initial_count = assembly_graph.graph_fragment.edges.len();
        
        assembly_graph.graph_fragment.edges.retain(|edge| edge.confidence >= confidence_threshold);
        
        let removed_count = initial_count - assembly_graph.graph_fragment.edges.len();
        println!("   Removed {} low-confidence edges", removed_count);
        
        Ok(())
    }
}

/// Complete assembly graph structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssemblyGraph {
    pub graph_fragment: GraphFragment,
    #[serde(skip)] // petgraph doesn't implement Serialize
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

impl AssemblyGraph {
    /// Generate contigs from the assembly graph
    pub fn generate_contigs(&mut self) -> Result<()> {
        println!("📝 Generating contigs from assembly graph");
        
        // Find strongly connected components
        let sccs = tarjan_scc(&self.petgraph);
        
        let mut contigs = Vec::new();
        let mut contig_id = 0;
        
        for component in sccs {
            if component.len() == 1 {
                // Linear contig
                if let Some(contig) = self.build_linear_contig(contig_id, &component)? {
                    contigs.push(contig);
                    contig_id += 1;
                }
            } else if component.len() > 1 {
                // Potential circular contig or complex structure
                if let Some(contig) = self.build_complex_contig(contig_id, &component)? {
                    contigs.push(contig);
                    contig_id += 1;
                }
            }
        }
        
        // Sort contigs by length (descending)
        contigs.sort_by(|a, b| b.length.cmp(&a.length));
        
        // Update contig IDs after sorting
        for (i, contig) in contigs.iter_mut().enumerate() {
            contig.id = i;
        }
        
        self.contigs = contigs;
        self.calculate_assembly_stats();
        
        println!("   Generated {} contigs", self.contigs.len());
        Ok(())
    }
    
    fn build_linear_contig(&self, contig_id: usize, component: &[NodeIndex]) -> Result<Option<Contig>> {
        if component.is_empty() {
            return Ok(None);
        }
        
        let node_idx = component[0];
        if let Some(&node_hash) = self.petgraph.node_weight(node_idx) {
            if let Some(node) = self.graph_fragment.nodes.get(&node_hash) {
                let contig = Contig {
                    id: contig_id,
                    sequence: node.kmer.sequence.clone(),
                    coverage: node.coverage as f64,
                    length: node.kmer.sequence.len(),
                    node_path: vec![node_hash],
                    contig_type: ContigType::Linear,
                };
                
                return Ok(Some(contig));
            }
        }
        
        Ok(None)
    }
    
    fn build_complex_contig(&self, contig_id: usize, component: &[NodeIndex]) -> Result<Option<Contig>> {
        // Try to find an Eulerian path through the component
        if let Some(path) = self.find_eulerian_path(component)? {
            let sequence = self.reconstruct_sequence_from_path(&path)?;
            let coverage = self.calculate_path_coverage_from_hashes(&path);
            
            let contig_type = if self.is_circular_path(&path) {
                ContigType::Circular
            } else {
                ContigType::Scaffold
            };
            
            let contig = Contig {
                id: contig_id,
                sequence,
                coverage,
                length: path.len() * self.graph_fragment.nodes.values().next()
                    .map(|n| n.kmer_size).unwrap_or(0),
                node_path: path,
                contig_type,
            };
            
            return Ok(Some(contig));
        }
        
        Ok(None)
    }
    
    fn find_eulerian_path(&self, component: &[NodeIndex]) -> Result<Option<Vec<u64>>> {
        // Simplified Eulerian path finding
        // In a real implementation, you'd use a proper algorithm like Hierholzer's
        
        if component.is_empty() {
            return Ok(None);
        }
        
        // For now, just create a simple path through the component
        let mut path = Vec::new();
        
        for &node_idx in component {
            if let Some(&node_hash) = self.petgraph.node_weight(node_idx) {
                path.push(node_hash);
            }
        }
        
        Ok(if path.is_empty() { None } else { Some(path) })
    }
    
    fn reconstruct_sequence_from_path(&self, path: &[u64]) -> Result<String> {
        if path.is_empty() {
            return Ok(String::new());
        }
        
        let first_node = self.graph_fragment.nodes.get(&path[0])
            .ok_or_else(|| anyhow!("First node in path not found"))?;
        
        let mut sequence = first_node.kmer.sequence.clone();
        
        for &node_hash in &path[1..] {
            if let Some(node) = self.graph_fragment.nodes.get(&node_hash) {
                // Assume k-1 overlap between consecutive k-mers
                let overlap = first_node.kmer_size - 1;
                if node.kmer.sequence.len() > overlap {
                    sequence.push_str(&node.kmer.sequence[overlap..]);
                }
            }
        }
        
        Ok(sequence)
    }
    
    fn calculate_path_coverage_from_hashes(&self, path: &[u64]) -> f64 {
        let total_coverage: u32 = path.iter()
            .filter_map(|&hash| self.graph_fragment.nodes.get(&hash))
            .map(|node| node.coverage)
            .sum();
        
        if path.is_empty() {
            0.0
        } else {
            total_coverage as f64 / path.len() as f64
        }
    }
    
    fn is_circular_path(&self, path: &[u64]) -> bool {
        if path.len() < 3 {
            return false;
        }
        
        // Check if last node connects back to first node
        let first = path[0];
        let last = path[path.len() - 1];
        
        self.graph_fragment.edges.iter().any(|edge| 
            edge.from_hash == last && edge.to_hash == first
        )
    }
    
    fn calculate_assembly_stats(&mut self) {
        if self.contigs.is_empty() {
            return;
        }
        
        // Basic statistics
        self.assembly_stats.total_contigs = self.contigs.len();
        self.assembly_stats.total_length = self.contigs.iter().map(|c| c.length).sum();
        self.assembly_stats.longest_contig = self.contigs.iter().map(|c| c.length).max().unwrap_or(0);
        
        // Calculate N50 and N90
        let mut sorted_lengths: Vec<usize> = self.contigs.iter().map(|c| c.length).collect();
        sorted_lengths.sort_unstable_by(|a, b| b.cmp(a)); // Descending order
        
        let total_length = self.assembly_stats.total_length;
        let mut cumulative_length = 0;
        
        for &length in &sorted_lengths {
            cumulative_length += length;
            
            if self.assembly_stats.n50 == 0 && cumulative_length >= total_length / 2 {
                self.assembly_stats.n50 = length;
            }
            
            if self.assembly_stats.n90 == 0 && cumulative_length >= (total_length * 9) / 10 {
                self.assembly_stats.n90 = length;
                break;
            }
        }
        
        // Calculate mean coverage
        let total_coverage: f64 = self.contigs.iter()
            .map(|c| c.coverage * c.length as f64)
            .sum();
        self.assembly_stats.mean_coverage = total_coverage / total_length as f64;
        
        // Calculate GC content
        let total_gc: usize = self.contigs.iter()
            .map(|c| c.sequence.chars().filter(|&ch| ch == 'G' || ch == 'C').count())
            .sum();
        self.assembly_stats.gc_content = total_gc as f64 / total_length as f64;
        
        println!("📊 Assembly Statistics:");
        println!("   Total contigs: {}", self.assembly_stats.total_contigs);
        println!("   Total length: {} bp", self.assembly_stats.total_length);
        println!("   Longest contig: {} bp", self.assembly_stats.longest_contig);
        println!("   N50: {} bp", self.assembly_stats.n50);
        println!("   N90: {} bp", self.assembly_stats.n90);
        println!("   Mean coverage: {:.2}x", self.assembly_stats.mean_coverage);
        println!("   GC content: {:.2}%", self.assembly_stats.gc_content * 100.0);
    }
    
    /// Export contigs to FASTA format
    pub fn write_contigs_fasta(&self, output_path: &std::path::Path) -> Result<()> {
        use std::io::Write;
        
        let mut file = std::fs::File::create(output_path)?;
        
        for contig in &self.contigs {
            writeln!(file, ">contig_{} length={} coverage={:.2}", 
                contig.id, contig.length, contig.coverage)?;
            
            // Write sequence in 80-character lines
            for chunk in contig.sequence.as_bytes().chunks(80) {
                writeln!(file, "{}", std::str::from_utf8(chunk)?)?;
            }
        }
        
        println!("✅ Wrote {} contigs to {}", self.contigs.len(), output_path.display());
        Ok(())
    }
    
    /// Export assembly graph to GFA format
    pub fn write_graph_gfa(&self, output_path: &std::path::Path) -> Result<()> {
        use std::io::Write;
        
        let mut file = std::fs::File::create(output_path)?;
        
        // Write GFA header
        writeln!(file, "H\tVN:Z:1.0")?;
        
        // Write segments (nodes)
        for (hash, node) in &self.graph_fragment.nodes {
            writeln!(file, "S\t{}\t{}\tLN:i:{}\tRC:i:{}",
                hash, node.kmer.sequence, node.kmer.sequence.len(), node.coverage)?;
        }
        
        // Write links (edges)
        for edge in &self.graph_fragment.edges {
            writeln!(file, "L\t{}\t+\t{}\t+\t{}M\tRC:i:{}",
                edge.from_hash, edge.to_hash, edge.overlap_length, edge.weight)?;
        }
        
        println!("✅ Wrote assembly graph to {}", output_path.display());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core_data_structures::*;
    
    #[test]
    fn test_assembly_graph_builder() {
        let builder = AssemblyGraphBuilder::new(15, 25, 2, 2).unwrap();
        
        // Create test reads
        let reads = vec![
            CorrectedRead {
                id: 0,
                original: "ATCGATCGATCGATCG".to_string(),
                corrected: "ATCGATCGATCGATCG".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30; 16],
            },
            CorrectedRead {
                id: 1,
                original: "TCGATCGATCGATCGA".to_string(),
                corrected: "TCGATCGATCGATCGA".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30; 16],
            },
        ];
        
        let assembly_graph = builder.build_graph(&reads).unwrap();
        
        assert!(assembly_graph.graph_fragment.nodes.len() > 0);
        assert!(assembly_graph.graph_fragment.edges.len() > 0);
    }
    
    #[test]
    fn test_contig_generation() {
        let mut assembly_graph = create_test_assembly_graph();
        assembly_graph.generate_contigs().unwrap();
        
        assert!(assembly_graph.contigs.len() > 0);
        assert!(assembly_graph.assembly_stats.total_length > 0);
        assert!(assembly_graph.assembly_stats.n50 > 0);
    }
    
    #[test]
    fn test_graph_simplification() {
        let builder = AssemblyGraphBuilder::new(15, 25, 1, 1).unwrap();
        let mut assembly_graph = create_test_assembly_graph();
        
        let initial_nodes = assembly_graph.graph_fragment.nodes.len();
        let initial_edges = assembly_graph.graph_fragment.edges.len();
        
        let simplified = builder.simplify_graph(assembly_graph).unwrap();
        
        // Simplification should reduce complexity
        assert!(simplified.graph_fragment.nodes.len() <= initial_nodes);
        assert!(simplified.graph_fragment.edges.len() <= initial_edges);
    }
    
    fn create_test_assembly_graph() -> AssemblyGraph {
        let mut fragment = GraphFragment::new(0);
        
        // Create test nodes
        for i in 0..5 {
            let sequence = format!("ATCG{}", "ATCG".repeat(i));
            let kmer = CanonicalKmer::new(&sequence).unwrap();
            let node = GraphNode::new(kmer, sequence.len());
            fragment.add_node(node);
        }
        
        // Create test edges
        let node_hashes: Vec<u64> = fragment.nodes.keys().copied().collect();
        for i in 0..node_hashes.len()-1 {
            let edge = GraphEdge::new(node_hashes[i], node_hashes[i+1], 3);
            fragment.add_edge(edge);
        }
        
        let petgraph = Graph::new();
        
        AssemblyGraph {
            graph_fragment: fragment,
            petgraph,
            contigs: Vec::new(),
            assembly_stats: AssemblyStats::default(),
        }
    }
}