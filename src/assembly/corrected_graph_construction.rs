//! Corrected De Bruijn Graph Construction for Metagenomic Assembly
//!
//! This module implements biologically-sound graph construction algorithms that:
//! 1. Merge overlapping k-mers from different reads representing same genomic regions
//! 2. Use coverage thresholds for error correction and quality control
//! 3. Handle sequencing errors and ambiguous bases appropriately
//! 4. Create proper many-to-one read:contig relationships

use ahash::{AHashMap, AHashSet};
use anyhow::{anyhow, Result};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

use crate::core::data_structures::{CorrectedRead, CanonicalKmer, GraphNode, GraphEdge};

/// Biologically-aware k-mer node with coverage tracking
#[derive(Debug, Clone)]
pub struct BiologicalKmerNode {
    pub kmer_hash: u64,
    pub sequence: String,
    pub coverage: u32,
    pub quality_scores: Vec<f64>, // Track quality from multiple reads
    pub supporting_reads: Vec<usize>, // Read IDs that contribute this k-mer
    pub error_corrected: bool,
}

impl BiologicalKmerNode {
    pub fn new(kmer_hash: u64, sequence: String, read_id: usize, quality: f64) -> Self {
        Self {
            kmer_hash,
            sequence,
            coverage: 1,
            quality_scores: vec![quality],
            supporting_reads: vec![read_id],
            error_corrected: false,
        }
    }
    
    /// Add support from another read
    pub fn add_read_support(&mut self, read_id: usize, quality: f64) {
        self.coverage += 1;
        self.quality_scores.push(quality);
        self.supporting_reads.push(read_id);
    }
    
    /// Calculate mean quality score across all supporting reads
    pub fn mean_quality(&self) -> f64 {
        if self.quality_scores.is_empty() {
            0.0
        } else {
            self.quality_scores.iter().sum::<f64>() / self.quality_scores.len() as f64
        }
    }
    
    /// Check if k-mer passes biological quality thresholds
    pub fn passes_biological_thresholds(&self, min_coverage: u32, min_quality: f64) -> bool {
        self.coverage >= min_coverage && self.mean_quality() >= min_quality
    }
}

/// Edge representing valid k-mer transitions with biological validation
#[derive(Debug, Clone)]
pub struct BiologicalEdge {
    pub from_kmer: u64,
    pub to_kmer: u64,
    pub support_count: u32,
    pub transition_quality: f64,
    pub supporting_reads: AHashSet<usize>,
}

impl BiologicalEdge {
    pub fn new(from_kmer: u64, to_kmer: u64, read_id: usize, quality: f64) -> Self {
        let mut supporting_reads = AHashSet::new();
        supporting_reads.insert(read_id);
        
        Self {
            from_kmer,
            to_kmer,
            support_count: 1,
            transition_quality: quality,
            supporting_reads,
        }
    }
    
    pub fn add_support(&mut self, read_id: usize, quality: f64) {
        if self.supporting_reads.insert(read_id) {
            self.support_count += 1;
            // Update quality as weighted average
            self.transition_quality = (self.transition_quality + quality) / 2.0;
        }
    }
}

/// Corrected De Bruijn graph with proper biological semantics
pub struct CorrectedDeBruijnGraph {
    k: usize,
    nodes: AHashMap<u64, BiologicalKmerNode>,
    edges: AHashMap<(u64, u64), BiologicalEdge>,
    min_coverage: u32,
    min_quality: f64,
    error_correction_enabled: bool,
    
    // Biological parameters
    max_error_rate: f64,
    gc_content_range: (f64, f64), // Valid GC content range
    complexity_threshold: f64,
}

impl CorrectedDeBruijnGraph {
    /// Create new biologically-aware De Bruijn graph
    pub fn new(k: usize, min_coverage: u32, min_quality: f64) -> Self {
        Self {
            k,
            nodes: AHashMap::new(),
            edges: AHashMap::new(),
            min_coverage,
            min_quality,
            error_correction_enabled: true,
            max_error_rate: 0.05, // 5% maximum error rate
            gc_content_range: (0.2, 0.8), // 20-80% GC content
            complexity_threshold: 0.3, // Minimum sequence complexity
        }
    }
    
    /// Process all reads to build unified De Bruijn graph
    pub fn build_from_reads(&mut self, reads: &[CorrectedRead]) -> Result<()> {
        println!("🔬 Building corrected De Bruijn graph from {} reads", reads.len());
        
        // Phase 1: Add all k-mers with coverage tracking
        for (read_id, read) in reads.iter().enumerate() {
            self.add_read_to_graph(read, read_id)?;
        }
        
        println!("   Initial k-mers: {}", self.nodes.len());
        
        // Phase 2: Error correction based on coverage
        if self.error_correction_enabled {
            self.perform_coverage_based_error_correction()?;
        }
        
        // Phase 3: Remove low-quality k-mers and edges
        self.filter_by_biological_criteria()?;
        
        println!("   Final k-mers after filtering: {}", self.nodes.len());
        println!("   Edges: {}", self.edges.len());
        
        Ok(())
    }
    
    /// Add single read to graph with proper k-mer merging
    fn add_read_to_graph(&mut self, read: &CorrectedRead, read_id: usize) -> Result<()> {
        if read.corrected.len() < self.k {
            return Ok(()); // Skip short reads
        }
        
        let sequence_bytes = read.corrected.as_bytes();
        let qualities = &read.quality_scores;
        let mut prev_kmer_hash: Option<u64> = None;
        
        // Extract all k-mers from this read
        for (i, window) in sequence_bytes.windows(self.k).enumerate() {
            if let Ok(kmer_str) = std::str::from_utf8(window) {
                // Validate k-mer quality before processing
                if self.is_valid_biological_kmer(kmer_str, qualities, i)? {
                    if let Ok(kmer) = CanonicalKmer::new(kmer_str) {
                        let kmer_hash = kmer.hash;
                        
                        // Calculate average quality for this k-mer
                        let kmer_quality = if i + self.k <= qualities.len() {
                            qualities[i..i + self.k].iter().sum::<u8>() as f64 / self.k as f64
                        } else {
                            30.0 // Default quality
                        };
                        
                        // Add or update k-mer node
                        self.add_or_update_kmer_node(kmer_hash, kmer_str.to_string(), read_id, kmer_quality);
                        
                        // Add edge between consecutive k-mers
                        if let Some(prev_hash) = prev_kmer_hash {
                            self.add_or_update_edge(prev_hash, kmer_hash, read_id, kmer_quality);
                        }
                        
                        prev_kmer_hash = Some(kmer_hash);
                    }
                }
            }
        }
        
        Ok(())
    }
    
    /// Validate k-mer from biological perspective
    fn is_valid_biological_kmer(&self, kmer_str: &str, qualities: &[u8], position: usize) -> Result<bool> {
        // Check 1: No invalid characters
        for c in kmer_str.chars() {
            match c.to_ascii_uppercase() {
                'A' | 'T' | 'G' | 'C' => {},
                'N' => return Ok(false), // Skip ambiguous bases for now
                _ => return Ok(false),
            }
        }
        
        // Check 2: Minimum quality threshold
        if position + self.k <= qualities.len() {
            let avg_quality = qualities[position..position + self.k].iter().sum::<u8>() as f64 / self.k as f64;
            if avg_quality < self.min_quality {
                return Ok(false);
            }
        }
        
        // Check 3: GC content in reasonable range
        let gc_content = self.calculate_gc_content(kmer_str);
        if gc_content < self.gc_content_range.0 || gc_content > self.gc_content_range.1 {
            return Ok(false);
        }
        
        // Check 4: Sequence complexity (avoid homopolymers)
        let complexity = self.calculate_sequence_complexity(kmer_str);
        if complexity < self.complexity_threshold {
            return Ok(false);
        }
        
        Ok(true)
    }
    
    /// Calculate GC content
    fn calculate_gc_content(&self, sequence: &str) -> f64 {
        let total = sequence.len() as f64;
        let gc_count = sequence.chars()
            .filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c')
            .count() as f64;
        gc_count / total
    }
    
    /// Calculate sequence complexity using Shannon entropy
    fn calculate_sequence_complexity(&self, sequence: &str) -> f64 {
        let mut counts = [0usize; 4]; // A, T, G, C
        let mut total = 0;
        
        for c in sequence.chars() {
            match c.to_ascii_uppercase() {
                'A' => { counts[0] += 1; total += 1; }
                'T' => { counts[1] += 1; total += 1; }
                'G' => { counts[2] += 1; total += 1; }
                'C' => { counts[3] += 1; total += 1; }
                _ => {}
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
        
        entropy / 2.0 // Normalize by maximum entropy
    }
    
    /// Add or update k-mer node with coverage tracking
    fn add_or_update_kmer_node(&mut self, kmer_hash: u64, sequence: String, read_id: usize, quality: f64) {
        match self.nodes.get_mut(&kmer_hash) {
            Some(node) => {
                // K-mer already exists - increment coverage
                node.add_read_support(read_id, quality);
            }
            None => {
                // New k-mer - create node
                let node = BiologicalKmerNode::new(kmer_hash, sequence, read_id, quality);
                self.nodes.insert(kmer_hash, node);
            }
        }
    }
    
    /// Add or update edge with support tracking
    fn add_or_update_edge(&mut self, from_kmer: u64, to_kmer: u64, read_id: usize, quality: f64) {
        let edge_key = (from_kmer, to_kmer);
        
        match self.edges.get_mut(&edge_key) {
            Some(edge) => {
                edge.add_support(read_id, quality);
            }
            None => {
                let edge = BiologicalEdge::new(from_kmer, to_kmer, read_id, quality);
                self.edges.insert(edge_key, edge);
            }
        }
    }
    
    /// Perform coverage-based error correction
    fn perform_coverage_based_error_correction(&mut self) -> Result<()> {
        println!("🧬 Performing coverage-based error correction");
        
        let mut corrected_count = 0;
        let low_coverage_kmers: Vec<u64> = self.nodes
            .iter()
            .filter(|(_, node)| node.coverage < self.min_coverage)
            .map(|(&hash, _)| hash)
            .collect();
        
        // Try to correct low-coverage k-mers by finding high-coverage neighbors
        for low_cov_kmer in low_coverage_kmers {
            if let Some(corrected_kmer) = self.find_error_correction_candidate(low_cov_kmer) {
                self.merge_kmers(low_cov_kmer, corrected_kmer)?;
                corrected_count += 1;
            }
        }
        
        println!("   Corrected {} k-mers through coverage analysis", corrected_count);
        Ok(())
    }
    
    /// Find high-coverage k-mer that could be error-corrected version
    fn find_error_correction_candidate(&self, low_cov_kmer: u64) -> Option<u64> {
        if let Some(low_node) = self.nodes.get(&low_cov_kmer) {
            // Generate single-base mutations and check for high-coverage matches
            for high_cov_kmer in self.nodes.keys() {
                if let Some(high_node) = self.nodes.get(high_cov_kmer) {
                    if high_node.coverage >= self.min_coverage * 3 && // Much higher coverage
                       self.is_single_error_apart(&low_node.sequence, &high_node.sequence) {
                        return Some(*high_cov_kmer);
                    }
                }
            }
        }
        None
    }
    
    /// Check if two k-mers differ by single substitution
    fn is_single_error_apart(&self, seq1: &str, seq2: &str) -> bool {
        if seq1.len() != seq2.len() {
            return false;
        }
        
        let mut diff_count = 0;
        for (c1, c2) in seq1.chars().zip(seq2.chars()) {
            if c1 != c2 {
                diff_count += 1;
                if diff_count > 1 {
                    return false;
                }
            }
        }
        diff_count == 1
    }
    
    /// Merge low-coverage k-mer into high-coverage k-mer
    fn merge_kmers(&mut self, low_cov_kmer: u64, high_cov_kmer: u64) -> Result<()> {
        if let (Some(low_node), Some(high_node)) = (
            self.nodes.remove(&low_cov_kmer), 
            self.nodes.get_mut(&high_cov_kmer)
        ) {
            // Transfer coverage and read support
            high_node.coverage += low_node.coverage;
            high_node.quality_scores.extend(low_node.quality_scores);
            high_node.supporting_reads.extend(low_node.supporting_reads);
            high_node.error_corrected = true;
            
            // Update edges to point to corrected k-mer
            let edges_to_update: Vec<_> = self.edges.keys()
                .filter(|&(from, to)| *from == low_cov_kmer || *to == low_cov_kmer)
                .copied()
                .collect();
            
            for (from, to) in edges_to_update {
                if let Some(edge) = self.edges.remove(&(from, to)) {
                    let new_key = if from == low_cov_kmer {
                        (high_cov_kmer, to)
                    } else {
                        (from, high_cov_kmer)
                    };
                    
                    // Merge with existing edge if it exists
                    match self.edges.get_mut(&new_key) {
                        Some(existing_edge) => {
                            existing_edge.support_count += edge.support_count;
                            existing_edge.supporting_reads.extend(edge.supporting_reads);
                        }
                        None => {
                            self.edges.insert(new_key, edge);
                        }
                    }
                }
            }
        }
        Ok(())
    }
    
    /// Filter k-mers and edges by biological criteria
    fn filter_by_biological_criteria(&mut self) -> Result<()> {
        println!("🧬 Filtering by biological quality criteria");
        
        let initial_nodes = self.nodes.len();
        let initial_edges = self.edges.len();
        
        // Remove low-quality k-mers
        self.nodes.retain(|_, node| {
            node.passes_biological_thresholds(self.min_coverage, self.min_quality)
        });
        
        // Remove edges involving removed k-mers
        let valid_kmers: AHashSet<u64> = self.nodes.keys().copied().collect();
        self.edges.retain(|(from, to), edge| {
            valid_kmers.contains(from) && 
            valid_kmers.contains(to) && 
            edge.support_count >= self.min_coverage
        });
        
        let filtered_nodes = initial_nodes - self.nodes.len();
        let filtered_edges = initial_edges - self.edges.len();
        
        println!("   Filtered {} low-quality k-mers", filtered_nodes);
        println!("   Filtered {} low-support edges", filtered_edges);
        
        Ok(())
    }
    
    /// Generate contigs from connected components with proper biological assembly
    pub fn generate_biological_contigs(&self) -> Result<Vec<BiologicalContig>> {
        println!("🧬 Generating biologically-meaningful contigs");
        
        let mut contigs = Vec::new();
        let mut visited = AHashSet::new();
        
        // Find connected components and assemble into contigs
        for &kmer_hash in self.nodes.keys() {
            if !visited.contains(&kmer_hash) {
                if let Some(contig) = self.assemble_contig_from_seed(kmer_hash, &mut visited)? {
                    if contig.length >= 100 { // Minimum contig length
                        contigs.push(contig);
                    }
                }
            }
        }
        
        // Sort contigs by length (longest first)
        contigs.sort_by_key(|c| std::cmp::Reverse(c.length));
        
        println!("✅ Generated {} biological contigs", contigs.len());
        Ok(contigs)
    }
    
    /// Assemble single contig starting from seed k-mer
    fn assemble_contig_from_seed(
        &self, 
        seed_kmer: u64, 
        visited: &mut AHashSet<u64>
    ) -> Result<Option<BiologicalContig>> {
        let mut path = Vec::new();
        let mut current_kmer = seed_kmer;
        let mut total_coverage = 0.0;
        
        visited.insert(seed_kmer);
        path.push(seed_kmer);
        
        if let Some(seed_node) = self.nodes.get(&seed_kmer) {
            total_coverage += seed_node.coverage as f64;
        }
        
        // Extend in forward direction
        while let Some(next_kmer) = self.get_best_successor(current_kmer, visited) {
            visited.insert(next_kmer);
            path.push(next_kmer);
            
            if let Some(node) = self.nodes.get(&next_kmer) {
                total_coverage += node.coverage as f64;
            }
            
            current_kmer = next_kmer;
        }
        
        // Extend in reverse direction from seed
        current_kmer = seed_kmer;
        let mut reverse_path = Vec::new();
        
        while let Some(prev_kmer) = self.get_best_predecessor(current_kmer, visited) {
            visited.insert(prev_kmer);
            reverse_path.push(prev_kmer);
            
            if let Some(node) = self.nodes.get(&prev_kmer) {
                total_coverage += node.coverage as f64;
            }
            
            current_kmer = prev_kmer;
        }
        
        // Combine reverse and forward paths
        reverse_path.reverse();
        reverse_path.extend(path);
        path = reverse_path;
        
        if path.is_empty() {
            return Ok(None);
        }
        
        // Construct contig sequence
        let sequence = self.reconstruct_sequence_from_path(&path)?;
        let mean_coverage = total_coverage / path.len() as f64;
        
        Ok(Some(BiologicalContig {
            id: path[0] as usize, // Use first k-mer hash as ID
            kmer_path: path.clone(),
            sequence,
            length: sequence.len(),
            coverage: mean_coverage,
            supporting_reads: self.get_supporting_reads_for_path(&path),
            gc_content: self.calculate_gc_content(&sequence),
        }))
    }
    
    /// Find best successor k-mer based on coverage and quality
    fn get_best_successor(&self, kmer: u64, visited: &AHashSet<u64>) -> Option<u64> {
        let mut best_successor = None;
        let mut best_score = 0.0;
        
        for ((from, to), edge) in &self.edges {
            if *from == kmer && !visited.contains(to) {
                if let Some(node) = self.nodes.get(to) {
                    let score = (edge.support_count as f64) * node.mean_quality();
                    if score > best_score {
                        best_score = score;
                        best_successor = Some(*to);
                    }
                }
            }
        }
        
        best_successor
    }
    
    /// Find best predecessor k-mer based on coverage and quality
    fn get_best_predecessor(&self, kmer: u64, visited: &AHashSet<u64>) -> Option<u64> {
        let mut best_predecessor = None;
        let mut best_score = 0.0;
        
        for ((from, to), edge) in &self.edges {
            if *to == kmer && !visited.contains(from) {
                if let Some(node) = self.nodes.get(from) {
                    let score = (edge.support_count as f64) * node.mean_quality();
                    if score > best_score {
                        best_score = score;
                        best_predecessor = Some(*from);
                    }
                }
            }
        }
        
        best_predecessor
    }
    
    /// Reconstruct DNA sequence from k-mer path
    fn reconstruct_sequence_from_path(&self, path: &[u64]) -> Result<String> {
        if path.is_empty() {
            return Ok(String::new());
        }
        
        let first_kmer = self.nodes.get(&path[0])
            .ok_or_else(|| anyhow!("First k-mer not found in graph"))?;
        let mut sequence = first_kmer.sequence.clone();
        
        // Add last character from each subsequent k-mer
        for &kmer_hash in &path[1..] {
            let node = self.nodes.get(&kmer_hash)
                .ok_or_else(|| anyhow!("K-mer not found in graph: {}", kmer_hash))?;
            
            if let Some(last_char) = node.sequence.chars().last() {
                sequence.push(last_char);
            }
        }
        
        Ok(sequence)
    }
    
    /// Get all reads supporting a contig path
    fn get_supporting_reads_for_path(&self, path: &[u64]) -> AHashSet<usize> {
        let mut supporting_reads = AHashSet::new();
        
        for &kmer_hash in path {
            if let Some(node) = self.nodes.get(&kmer_hash) {
                supporting_reads.extend(&node.supporting_reads);
            }
        }
        
        supporting_reads
    }
}

/// Biologically-meaningful contig with proper assembly metadata
#[derive(Debug, Clone)]
pub struct BiologicalContig {
    pub id: usize,
    pub kmer_path: Vec<u64>,
    pub sequence: String,
    pub length: usize,
    pub coverage: f64,
    pub supporting_reads: AHashSet<usize>,
    pub gc_content: f64,
}

impl BiologicalContig {
    /// Get number of reads that contributed to this contig
    pub fn read_support_count(&self) -> usize {
        self.supporting_reads.len()
    }
    
    /// Check if contig meets biological quality criteria
    pub fn is_high_quality(&self, min_length: usize, min_coverage: f64) -> bool {
        self.length >= min_length && 
        self.coverage >= min_coverage && 
        self.gc_content >= 0.2 && 
        self.gc_content <= 0.8
    }
}