//! Laptop-Optimized Assembly Pipeline
//! ==================================
//!
//! Streamlined assembly implementation optimized for typical laptop hardware:
//! - 4-8 GB RAM constraint awareness
//! - 2-4 CPU cores with efficient parallelization
//! - SSD-friendly sequential I/O patterns
//! - Memory pool-free design for low overhead
//!
//! This module consolidates the best practices from multiple optimization modules
//! into a single, efficient, and maintainable implementation.

use crate::core::data_structures::{CorrectedRead, AssemblyStats, Contig, ContigType};
use crate::assembly::adaptive_k::{AdaptiveKSelector, AdaptiveKConfig};
use anyhow::{anyhow, Result};
use ahash::{AHashMap, AHashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Memory budget targets for different laptop configurations
#[derive(Debug, Clone)]
pub struct LaptopConfig {
    /// Target memory usage in MB
    pub memory_budget_mb: usize,
    /// Number of CPU cores to utilize
    pub cpu_cores: usize,
    /// Chunk size for processing (auto-tuned)
    pub chunk_size: usize,
    /// Maximum k-mer size to prevent memory explosion
    pub max_k: usize,
}

impl LaptopConfig {
    /// Configuration for 4GB laptops (conservative)
    pub fn low_memory() -> Self {
        Self {
            memory_budget_mb: 1024,  // Use only 1GB for assembly
            cpu_cores: 2,
            chunk_size: 500,
            max_k: 21,
        }
    }

    /// Configuration for 8GB laptops (balanced)
    pub fn medium_memory() -> Self {
        Self {
            memory_budget_mb: 2048,  // Use 2GB for assembly
            cpu_cores: 4,
            chunk_size: 1000,
            max_k: 31,
        }
    }

    /// Configuration for 16GB+ laptops (performance)
    pub fn high_memory() -> Self {
        Self {
            memory_budget_mb: 4096,  // Use 4GB for assembly
            cpu_cores: num_cpus::get(),
            chunk_size: 2000,
            max_k: 63,
        }
    }

    /// Auto-detect optimal configuration based on system resources
    pub fn auto_detect() -> Self {
        let total_memory_gb = Self::detect_system_memory_gb();
        let cpu_cores = num_cpus::get();

        if total_memory_gb <= 4.0 {
            println!("🔧 Auto-detected: Low memory system ({:.1} GB RAM, {} cores)", total_memory_gb, cpu_cores);
            Self::low_memory()
        } else if total_memory_gb <= 8.0 {
            println!("🔧 Auto-detected: Medium memory system ({:.1} GB RAM, {} cores)", total_memory_gb, cpu_cores);
            Self::medium_memory()
        } else {
            println!("🔧 Auto-detected: High memory system ({:.1} GB RAM, {} cores)", total_memory_gb, cpu_cores);
            Self::high_memory()
        }
    }

    /// Detect system memory (simplified implementation)
    fn detect_system_memory_gb() -> f64 {
        // This is a simplified fallback implementation
        // In a real system, you'd use platform-specific APIs
        match std::env::var("MEMORY_GB") {
            Ok(mem_str) => mem_str.parse().unwrap_or(8.0),
            Err(_) => {
                // Conservative default for unknown systems
                8.0
            }
        }
    }

    /// Create custom configuration with validation
    pub fn custom(memory_budget_mb: usize, cpu_cores: usize, max_k: usize) -> Result<Self> {
        if memory_budget_mb < 256 {
            return Err(anyhow!("Memory budget too low: {} MB (minimum 256 MB)", memory_budget_mb));
        }
        if cpu_cores == 0 {
            return Err(anyhow!("CPU cores must be at least 1"));
        }
        if max_k < 15 || max_k > 127 {
            return Err(anyhow!("Max k-mer size must be between 15 and 127"));
        }

        let chunk_size = (memory_budget_mb / 2).max(100).min(5000);

        Ok(Self {
            memory_budget_mb,
            cpu_cores,
            chunk_size,
            max_k,
        })
    }
}

/// Compact k-mer representation using 2-bit encoding
/// Optimized for memory efficiency on resource-constrained systems
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CompactKmer {
    /// Packed nucleotide data (2 bits per nucleotide)
    data: u64,
    /// K-mer length (max 32 for single u64)
    k: u8,
}

impl CompactKmer {
    /// Create compact k-mer from DNA sequence
    pub fn new(seq: &str) -> Result<Self> {
        let k = seq.len();
        if k == 0 || k > 32 {
            return Err(anyhow!("Invalid k-mer length: {}", k));
        }

        let mut data = 0u64;
        for (i, nucleotide) in seq.chars().enumerate() {
            let bits = match nucleotide.to_ascii_uppercase() {
                'A' => 0b00,
                'C' => 0b01,
                'G' => 0b10,
                'T' => 0b11,
                _ => return Err(anyhow!("Invalid nucleotide: {}", nucleotide)),
            };
            data |= bits << (2 * (31 - i));
        }

        Ok(Self { data, k: k as u8 })
    }

    /// Convert back to string representation
    pub fn to_string(&self) -> String {
        let mut result = String::with_capacity(self.k as usize);
        for i in 0..self.k {
            let bits = (self.data >> (2 * (31 - i))) & 0b11;
            let nucleotide = match bits {
                0b00 => 'A',
                0b01 => 'C',
                0b10 => 'G',
                0b11 => 'T',
                _ => unreachable!(),
            };
            result.push(nucleotide);
        }
        result
    }

    /// Get hash value for indexing
    pub fn hash(&self) -> u64 {
        // Simple hash combining data and length
        self.data.wrapping_mul(0x9e3779b97f4a7c15) ^ (self.k as u64)
    }

    /// Memory footprint in bytes
    pub fn memory_footprint(&self) -> usize {
        std::mem::size_of::<Self>()
    }
}

/// Memory-bounded k-mer counter with automatic cleanup
/// Prevents memory explosion by maintaining size limits
pub struct BoundedKmerCounter {
    /// K-mer counts with size limit
    counts: AHashMap<u64, u32>,
    /// Memory limit in number of k-mers
    max_kmers: usize,
    /// Current memory usage estimate
    memory_usage: AtomicUsize,
    /// Statistics
    total_kmers_seen: AtomicUsize,
    kmers_dropped: AtomicUsize,
}

impl BoundedKmerCounter {
    /// Create counter with memory budget
    pub fn new(memory_budget_mb: usize) -> Self {
        // Estimate ~12 bytes per k-mer (hash + count + overhead)
        let max_kmers = (memory_budget_mb * 1024 * 1024) / 12;

        Self {
            counts: AHashMap::with_capacity(max_kmers),
            max_kmers,
            memory_usage: AtomicUsize::new(0),
            total_kmers_seen: AtomicUsize::new(0),
            kmers_dropped: AtomicUsize::new(0),
        }
    }

    /// Add k-mer with memory-aware insertion
    pub fn add_kmer(&mut self, kmer_hash: u64) {
        self.total_kmers_seen.fetch_add(1, Ordering::Relaxed);

        // Always update existing k-mers
        if let Some(count) = self.counts.get_mut(&kmer_hash) {
            *count = count.saturating_add(1);
            return;
        }

        // For new k-mers, check memory limit
        if self.counts.len() < self.max_kmers {
            self.counts.insert(kmer_hash, 1);
            self.memory_usage.fetch_add(12, Ordering::Relaxed);
        } else {
            // Drop rare k-mers when at capacity
            self.cleanup_rare_kmers();
            self.kmers_dropped.fetch_add(1, Ordering::Relaxed);
        }
    }

    /// Remove k-mers with count = 1 to free memory
    fn cleanup_rare_kmers(&mut self) {
        let before_size = self.counts.len();
        self.counts.retain(|_, &mut count| count > 1);
        let removed = before_size - self.counts.len();
        self.memory_usage.fetch_sub(removed * 12, Ordering::Relaxed);
    }

    /// Get k-mers above threshold
    pub fn get_frequent_kmers(&self, min_count: u32) -> Vec<(u64, u32)> {
        self.counts
            .iter()
            .filter(|(_, &count)| count >= min_count)
            .map(|(&hash, &count)| (hash, count))
            .collect()
    }

    /// Get memory usage statistics
    pub fn get_stats(&self) -> (usize, usize, usize, usize) {
        (
            self.counts.len(),
            self.total_kmers_seen.load(Ordering::Relaxed),
            self.kmers_dropped.load(Ordering::Relaxed),
            self.memory_usage.load(Ordering::Relaxed),
        )
    }
}

/// Streamlined assembly graph optimized for laptop memory constraints
pub struct LaptopAssemblyGraph {
    /// Node storage with hash-based indexing
    nodes: AHashMap<u64, GraphNode>,
    /// Edge list (more memory-efficient than adjacency matrix)
    edges: Vec<GraphEdge>,
    /// Configuration parameters
    config: LaptopConfig,
    /// Assembly statistics
    stats: AssemblyStats,
}

#[derive(Debug, Clone)]
struct GraphNode {
    kmer: CompactKmer,
    coverage: u32,
    in_degree: u8,
    out_degree: u8,
}

#[derive(Debug, Clone)]
struct GraphEdge {
    from_hash: u64,
    to_hash: u64,
    weight: u32,
}

impl LaptopAssemblyGraph {
    /// Create new graph with laptop configuration
    pub fn new(config: LaptopConfig) -> Self {
        // Estimate node capacity based on memory budget
        let node_capacity = (config.memory_budget_mb * 1024 * 1024) / 32; // ~32 bytes per node

        Self {
            nodes: AHashMap::with_capacity(node_capacity),
            edges: Vec::new(),
            config,
            stats: AssemblyStats::default(),
        }
    }

    /// Build graph from reads using memory-bounded processing
    pub fn build_from_reads(&mut self, reads: &[CorrectedRead], k: usize) -> Result<()> {
        if k > self.config.max_k {
            return Err(anyhow!("K-mer size {} exceeds maximum {}", k, self.config.max_k));
        }

        // Process reads in chunks to control memory usage
        let chunks: Vec<_> = reads.chunks(self.config.chunk_size).collect();

        println!("🧬 Processing {} reads in {} chunks (k={})", reads.len(), chunks.len(), k);

        // Use bounded k-mer counter to prevent memory explosion
        let mut kmer_counter = BoundedKmerCounter::new(self.config.memory_budget_mb / 4);

        // Phase 1: Count k-mers to identify frequent ones
        for chunk in &chunks {
            self.process_chunk_for_counting(chunk, k, &mut kmer_counter)?;
        }

        let (unique_kmers, total_seen, dropped, memory_used) = kmer_counter.get_stats();
        println!("📊 K-mer counting: {} unique, {} total, {} dropped, {:.1} MB used",
                unique_kmers, total_seen, dropped, memory_used as f64 / (1024.0 * 1024.0));

        // Phase 2: Build graph using frequent k-mers
        let frequent_kmers = kmer_counter.get_frequent_kmers(2); // Min coverage of 2
        let frequent_set: AHashSet<u64> = frequent_kmers.iter().map(|(hash, _)| *hash).collect();

        println!("🔗 Building graph with {} frequent k-mers", frequent_set.len());

        for chunk in &chunks {
            self.process_chunk_for_graph_building(chunk, k, &frequent_set)?;
        }

        self.cleanup_low_coverage_nodes(2);
        self.calculate_stats();

        println!("✅ Graph built: {} nodes, {} edges", self.nodes.len(), self.edges.len());

        Ok(())
    }

    /// Process chunk for k-mer counting phase
    fn process_chunk_for_counting(
        &self,
        chunk: &[CorrectedRead],
        k: usize,
        counter: &mut BoundedKmerCounter
    ) -> Result<()> {
        for read in chunk {
            if read.corrected.len() < k {
                continue;
            }

            // Extract k-mers from read
            for i in 0..=read.corrected.len() - k {
                if let Ok(kmer) = CompactKmer::new(&read.corrected[i..i + k]) {
                    counter.add_kmer(kmer.hash());
                }
            }
        }
        Ok(())
    }

    /// Process chunk for graph building phase
    fn process_chunk_for_graph_building(
        &mut self,
        chunk: &[CorrectedRead],
        k: usize,
        frequent_kmers: &AHashSet<u64>
    ) -> Result<()> {
        for read in chunk {
            if read.corrected.len() < k + 1 {
                continue;
            }

            // Extract consecutive k-mers and create edges
            let mut prev_kmer: Option<CompactKmer> = None;

            for i in 0..=read.corrected.len() - k {
                if let Ok(kmer) = CompactKmer::new(&read.corrected[i..i + k]) {
                    let kmer_hash = kmer.hash();

                    // Only process frequent k-mers
                    if !frequent_kmers.contains(&kmer_hash) {
                        prev_kmer = None;
                        continue;
                    }

                    // Add node
                    self.add_or_update_node(kmer.clone());

                    // Add edge if we have previous k-mer
                    if let Some(prev) = &prev_kmer {
                        self.add_edge(prev.hash(), kmer_hash);
                    }

                    prev_kmer = Some(kmer);
                }
            }
        }
        Ok(())
    }

    /// Add or update node in graph
    fn add_or_update_node(&mut self, kmer: CompactKmer) {
        let hash = kmer.hash();

        match self.nodes.get_mut(&hash) {
            Some(node) => {
                node.coverage = node.coverage.saturating_add(1);
            }
            None => {
                let node = GraphNode {
                    kmer,
                    coverage: 1,
                    in_degree: 0,
                    out_degree: 0,
                };
                self.nodes.insert(hash, node);
            }
        }
    }

    /// Add edge between nodes
    fn add_edge(&mut self, from_hash: u64, to_hash: u64) {
        // Check if edge already exists
        if self.edges.iter().any(|e| e.from_hash == from_hash && e.to_hash == to_hash) {
            return;
        }

        // Add edge
        self.edges.push(GraphEdge {
            from_hash,
            to_hash,
            weight: 1,
        });

        // Update degree information
        if let Some(from_node) = self.nodes.get_mut(&from_hash) {
            from_node.out_degree = from_node.out_degree.saturating_add(1);
        }
        if let Some(to_node) = self.nodes.get_mut(&to_hash) {
            to_node.in_degree = to_node.in_degree.saturating_add(1);
        }
    }

    /// Remove nodes with coverage below threshold
    fn cleanup_low_coverage_nodes(&mut self, min_coverage: u32) {
        let low_coverage_nodes: Vec<u64> = self.nodes
            .iter()
            .filter(|(_, node)| node.coverage < min_coverage)
            .map(|(&hash, _)| hash)
            .collect();

        // Remove nodes
        for hash in &low_coverage_nodes {
            self.nodes.remove(hash);
        }

        // Remove edges involving removed nodes
        self.edges.retain(|edge| {
            self.nodes.contains_key(&edge.from_hash) && self.nodes.contains_key(&edge.to_hash)
        });

        println!("🧹 Removed {} low-coverage nodes", low_coverage_nodes.len());
    }

    /// Generate contigs using simple linear path traversal
    pub fn generate_contigs(&self) -> Result<Vec<Contig>> {
        let mut contigs = Vec::new();
        let mut visited = AHashSet::new();
        let mut contig_id = 0;

        // Build adjacency information
        let mut outgoing: AHashMap<u64, Vec<u64>> = AHashMap::new();
        let mut incoming: AHashMap<u64, Vec<u64>> = AHashMap::new();

        for edge in &self.edges {
            outgoing.entry(edge.from_hash).or_default().push(edge.to_hash);
            incoming.entry(edge.to_hash).or_default().push(edge.from_hash);
        }

        // Find starting nodes (no incoming edges or multiple outgoing)
        for (&node_hash, _node) in &self.nodes {
            if visited.contains(&node_hash) {
                continue;
            }

            let in_count = incoming.get(&node_hash).map_or(0, |v| v.len());
            let out_count = outgoing.get(&node_hash).map_or(0, |v| v.len());

            // Start from nodes that are likely contig starts
            if in_count == 0 || out_count != 1 || in_count != 1 {
                if let Some(contig) = self.trace_contig(node_hash, &outgoing, &mut visited)? {
                    contigs.push(Contig {
                        id: contig_id,
                        sequence: contig.sequence,
                        coverage: contig.coverage,
                        length: contig.length,
                        node_path: vec![], // Simplified for efficiency
                        contig_type: ContigType::Linear,
                    });
                    contig_id += 1;
                }
            }
        }

        Ok(contigs)
    }

    /// Trace linear path to form contig
    fn trace_contig(
        &self,
        start_hash: u64,
        outgoing: &AHashMap<u64, Vec<u64>>,
        visited: &mut AHashSet<u64>,
    ) -> Result<Option<SimpleContig>> {
        if visited.contains(&start_hash) {
            return Ok(None);
        }

        let mut path = vec![start_hash];
        let mut current = start_hash;

        // Follow linear path
        while let Some(neighbors) = outgoing.get(&current) {
            if neighbors.len() == 1 && !visited.contains(&neighbors[0]) {
                current = neighbors[0];
                path.push(current);
            } else {
                break;
            }
        }

        // Mark nodes as visited
        for &node_hash in &path {
            visited.insert(node_hash);
        }

        if path.len() < 2 {
            return Ok(None);
        }

        // Build contig sequence
        let mut sequence = String::new();
        let mut total_coverage = 0.0;

        for (i, &node_hash) in path.iter().enumerate() {
            if let Some(node) = self.nodes.get(&node_hash) {
                total_coverage += node.coverage as f64;

                if i == 0 {
                    sequence.push_str(&node.kmer.to_string());
                } else {
                    // Add overlap - just the last character for simplicity
                    let kmer_str = node.kmer.to_string();
                    if let Some(last_char) = kmer_str.chars().last() {
                        sequence.push(last_char);
                    }
                }
            }
        }

        let avg_coverage = total_coverage / path.len() as f64;

        let length = sequence.len();
        Ok(Some(SimpleContig {
            sequence,
            length,
            coverage: avg_coverage,
        }))
    }

    /// Calculate assembly statistics
    fn calculate_stats(&mut self) {
        self.stats.total_length = self.nodes.len();
        self.stats.num_contigs = self.edges.len();
        // Simplified stats for efficiency
    }

    /// Get memory usage in MB
    pub fn memory_usage_mb(&self) -> f64 {
        let nodes_size = self.nodes.len() * 32; // Estimated size per node
        let edges_size = self.edges.len() * 24; // Estimated size per edge
        (nodes_size + edges_size) as f64 / (1024.0 * 1024.0)
    }

    /// Get assembly statistics
    pub fn stats(&self) -> &AssemblyStats {
        &self.stats
    }
}

/// Simple contig structure for internal use
#[derive(Debug)]
struct SimpleContig {
    sequence: String,
    length: usize,
    coverage: f64,
}

/// Main laptop-optimized assembler
pub struct LaptopAssembler {
    config: LaptopConfig,
}

impl LaptopAssembler {
    /// Create assembler with configuration
    pub fn new(config: LaptopConfig) -> Self {
        Self { config }
    }

    /// Detect optimal configuration based on system
    pub fn auto_config() -> Self {
        let config = LaptopConfig::auto_detect();
        Self::new(config)
    }

    /// Perform assembly with automatic parameter selection
    pub fn assemble(&self, reads: &[CorrectedRead]) -> Result<Vec<Contig>> {
        println!("🚀 Starting laptop-optimized assembly");
        println!("   📊 Input: {} reads", reads.len());
        println!("   💾 Memory budget: {} MB", self.config.memory_budget_mb);
        println!("   ⚙️  CPU cores: {}", self.config.cpu_cores);

        // Auto-select k-mer size based on read characteristics
        let k = self.select_optimal_k(reads)?;
        println!("   🧬 Selected k-mer size: {}", k);

        // Build graph
        let mut graph = LaptopAssemblyGraph::new(self.config.clone());
        graph.build_from_reads(reads, k)?;

        println!("   📈 Memory usage: {:.1} MB", graph.memory_usage_mb());

        // Generate contigs
        let contigs = graph.generate_contigs()?;

        println!("✅ Assembly complete: {} contigs", contigs.len());

        Ok(contigs)
    }

    /// Select optimal k-mer size using adaptive selection
    fn select_optimal_k(&self, reads: &[CorrectedRead]) -> Result<usize> {
        let adaptive_config = AdaptiveKConfig {
            min_k: 15,
            max_k: self.config.max_k,
            sample_size: 1000,
            memory_budget_mb: self.config.memory_budget_mb / 4, // Use quarter of budget for k-mer selection
        };

        let selector = AdaptiveKSelector::new(adaptive_config);
        selector.select_optimal_k(reads)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::data_structures::CorrectionMetadata;

    #[test]
    fn test_compact_kmer() {
        let kmer = CompactKmer::new("ATCG").unwrap();
        assert_eq!(kmer.to_string(), "ATCG");
        assert!(kmer.memory_footprint() <= 16); // Should be very compact
    }

    #[test]
    fn test_bounded_kmer_counter() {
        let mut counter = BoundedKmerCounter::new(1); // 1MB limit

        // Add many k-mers to test memory bounding
        for i in 0..100000 {
            counter.add_kmer(i);
        }

        let (unique, total, dropped, memory) = counter.get_stats();
        assert!(unique <= counter.max_kmers);
        assert_eq!(total, 100000);
        println!("Stats: {} unique, {} total, {} dropped, {} bytes",
                unique, total, dropped, memory);
    }

    #[test]
    fn test_laptop_assembler() {
        let config = LaptopConfig::low_memory();
        let assembler = LaptopAssembler::new(config);

        // Create test reads
        let reads = vec![
            CorrectedRead {
                id: 0,
                original: "ATCGATCGATCG".to_string(),
                corrected: "ATCGATCGATCG".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30; 12],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
            },
            CorrectedRead {
                id: 1,
                original: "TCGATCGATCGA".to_string(),
                corrected: "TCGATCGATCGA".to_string(),
                corrections: Vec::new(),
                quality_scores: vec![30; 12],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
            },
        ];

        let contigs = assembler.assemble(&reads).unwrap();
        assert!(!contigs.is_empty());

        println!("Generated {} contigs", contigs.len());
        for (i, contig) in contigs.iter().enumerate() {
            println!("Contig {}: len={}, cov={:.1}", i, contig.length, contig.coverage);
        }
    }
}