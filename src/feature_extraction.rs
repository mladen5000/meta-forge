use std::collections::{HashMap, HashSet};
use ahash::{AHashMap, AHashSet};
use anyhow::{Result, anyhow};
use ndarray::{Array1, Array2};
use serde::{Serialize, Deserialize};
use rayon::prelude::*;

use crate::core_data_structures::*;
use crate::assembly_graph_construction::*;

/// Comprehensive feature extraction for sequences and assembly graphs
pub struct AdvancedFeatureExtractor {
    /// Configuration for feature extraction
    config: FeatureConfig,
    /// Codon usage tables for different organisms
    codon_tables: CodonUsageTables,
    /// Pre-computed k-mer signatures
    kmer_signatures: KmerSignatures,
    /// Sequence pattern recognizers
    pattern_recognizers: PatternRecognizers,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureConfig {
    /// Sequence features
    pub include_composition: bool,
    pub include_codon_usage: bool,
    pub include_patterns: bool,
    pub include_complexity: bool,
    
    /// Graph features
    pub include_topology: bool,
    pub include_centrality: bool,
    pub include_clustering: bool,
    
    /// k-mer based features
    pub kmer_sizes: Vec<usize>,
    pub max_kmers: usize,
    
    /// Output dimensions
    pub sequence_feature_dim: usize,
    pub graph_feature_dim: usize,
    pub kmer_feature_dim: usize,
}

impl Default for FeatureConfig {
    fn default() -> Self {
        Self {
            include_composition: true,
            include_codon_usage: true,
            include_patterns: true,
            include_complexity: true,
            include_topology: true,
            include_centrality: true,
            include_clustering: true,
            kmer_sizes: vec![3, 4, 5, 6],
            max_kmers: 10000,
            sequence_feature_dim: 100,
            graph_feature_dim: 50,
            kmer_feature_dim: 64,
        }
    }
}

/// Comprehensive feature vector for machine learning
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureVector {
    /// Sequence-based features
    pub sequence_features: Array1<f64>,
    /// Graph topology features
    pub graph_features: Array1<f64>,
    /// k-mer based features
    pub kmer_features: Array1<f64>,
    /// Combined feature vector
    pub combined_features: Array1<f64>,
    /// Feature metadata
    pub metadata: FeatureMetadata,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FeatureMetadata {
    pub sequence_length: usize,
    pub gc_content: f64,
    pub complexity_score: f64,
    pub feature_names: Vec<String>,
    pub extraction_time_ms: u64,
}

impl AdvancedFeatureExtractor {
    pub fn new(config: FeatureConfig) -> Result<Self> {
        let codon_tables = CodonUsageTables::load_standard_tables()?;
        let kmer_signatures = KmerSignatures::new(&config.kmer_sizes)?;
        let pattern_recognizers = PatternRecognizers::new()?;
        
        Ok(Self {
            config,
            codon_tables,
            kmer_signatures,
            pattern_recognizers,
        })
    }
    
    /// Extract comprehensive features from a DNA sequence
    pub fn extract_sequence_features(&self, sequence: &str) -> Result<FeatureVector> {
        let start_time = std::time::Instant::now();
        let mut feature_names = Vec::new();
        let mut features = Vec::new();
        
        // Basic composition features
        if self.config.include_composition {
            let comp_features = self.extract_composition_features(sequence)?;
            features.extend(comp_features);
            feature_names.extend(self.get_composition_feature_names());
        }
        
        // Codon usage features
        if self.config.include_codon_usage {
            let codon_features = self.extract_codon_usage_features(sequence)?;
            features.extend(codon_features);
            feature_names.extend(self.get_codon_feature_names());
        }
        
        // Pattern recognition features
        if self.config.include_patterns {
            let pattern_features = self.extract_pattern_features(sequence)?;
            features.extend(pattern_features);
            feature_names.extend(self.get_pattern_feature_names());
        }
        
        // Complexity features
        if self.config.include_complexity {
            let complexity_features = self.extract_complexity_features(sequence)?;
            features.extend(complexity_features);
            feature_names.extend(self.get_complexity_feature_names());
        }
        
        // k-mer features
        let kmer_features = self.extract_kmer_features(sequence)?;
        feature_names.extend(self.get_kmer_feature_names());
        
        // Pad or truncate to desired dimensions
        self.normalize_feature_dimensions(&mut features, self.config.sequence_feature_dim);
        
        let sequence_features = Array1::from_vec(features);
        let kmer_features_array = Array1::from_vec(kmer_features);
        let combined_features = self.combine_features(&sequence_features, &Array1::zeros(0), &kmer_features_array)?;
        
        let metadata = FeatureMetadata {
            sequence_length: sequence.len(),
            gc_content: calculate_gc_content(sequence),
            complexity_score: calculate_sequence_complexity(sequence),
            feature_names,
            extraction_time_ms: start_time.elapsed().as_millis() as u64,
        };
        
        Ok(FeatureVector {
            sequence_features,
            graph_features: Array1::zeros(0),
            kmer_features: kmer_features_array,
            combined_features,
            metadata,
        })
    }
    
    /// Extract features from assembly graph
    pub fn extract_graph_features(&self, graph: &AssemblyGraph) -> Result<Array1<f64>> {
        let mut features = Vec::new();
        
        if self.config.include_topology {
            features.extend(self.extract_topology_features(graph)?);
        }
        
        if self.config.include_centrality {
            features.extend(self.extract_centrality_features(graph)?);
        }
        
        if self.config.include_clustering {
            features.extend(self.extract_clustering_features(graph)?);
        }
        
        self.normalize_feature_dimensions(&mut features, self.config.graph_feature_dim);
        Ok(Array1::from_vec(features))
    }
    
    /// Extract features for a specific node in the graph
    pub fn extract_node_features(&self, node: &GraphNode, graph: &AssemblyGraph) -> Result<Array1<f64>> {
        let mut features = Vec::new();
        
        // Node-specific features
        features.push(node.coverage as f64);
        features.push(node.complexity_score);
        features.push(node.kmer.sequence.len() as f64);
        features.push(node.read_positions.len() as f64);
        features.push(calculate_gc_content(&node.kmer.sequence));
        
        // Node type encoding
        let node_type_encoding = match node.node_type {
            NodeType::Unique => [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            NodeType::Repetitive => [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            NodeType::LowCoverage => [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            NodeType::HighCoverage => [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            NodeType::Tip => [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            NodeType::Bubble => [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        };
        features.extend(node_type_encoding);
        
        // Local graph topology
        let local_features = self.extract_local_topology_features(node, graph)?;
        features.extend(local_features);
        
        // Sequence composition for this node
        let comp_features = self.extract_composition_features(&node.kmer.sequence)?;
        features.extend(comp_features.into_iter().take(10)); // Take first 10 composition features
        
        Ok(Array1::from_vec(features))
    }
    
    /// Basic nucleotide composition features
    fn extract_composition_features(&self, sequence: &str) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        let total_len = sequence.len() as f64;
        
        if total_len == 0.0 {
            return Ok(vec![0.0; 20]); // Return zeros for empty sequences
        }
        
        // Single nucleotide frequencies
        let mut nucleotide_counts = [0u32; 4]; // A, C, G, T
        for byte in sequence.bytes() {
            match byte {
                b'A' | b'a' => nucleotide_counts[0] += 1,
                b'C' | b'c' => nucleotide_counts[1] += 1,
                b'G' | b'g' => nucleotide_counts[2] += 1,
                b'T' | b't' => nucleotide_counts[3] += 1,
                _ => {} // Skip ambiguous bases
            }
        }
        
        // Add normalized single nucleotide frequencies
        for &count in &nucleotide_counts {
            features.push(count as f64 / total_len);
        }
        
        // Dinucleotide frequencies
        let dinucleotides = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                            "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"];
        
        let mut dinuc_counts = vec![0u32; 16];
        for window in sequence.as_bytes().windows(2) {
            if let Ok(dinuc) = std::str::from_utf8(window) {
                if let Some(pos) = dinucleotides.iter().position(|&x| x.eq_ignore_ascii_case(dinuc)) {
                    dinuc_counts[pos] += 1;
                }
            }
        }
        
        let total_dinucs = dinuc_counts.iter().sum::<u32>() as f64;
        if total_dinucs > 0.0 {
            for count in dinuc_counts {
                features.push(count as f64 / total_dinucs);
            }
        } else {
            features.extend(vec![0.0; 16]);
        }
        
        Ok(features)
    }
    
    /// Codon usage bias features
    fn extract_codon_usage_features(&self, sequence: &str) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        // Extract all codons
        let mut codon_counts = AHashMap::new();
        for window in sequence.as_bytes().windows(3) {
            if window.len() == 3 {
                if let Ok(codon) = std::str::from_utf8(window) {
                    let codon_upper = codon.to_uppercase();
                    *codon_counts.entry(codon_upper).or_insert(0) += 1;
                }
            }
        }
        
        let total_codons = codon_counts.values().sum::<u32>() as f64;
        
        if total_codons > 0.0 {
            // Calculate codon usage bias for each amino acid
            for amino_acid in &["F", "L", "S", "Y", "C", "W", "P", "H", "Q", "R", "I", "M", "T", "N", "K", "V", "A", "D", "E", "G"] {
                let aa_codons = self.codon_tables.get_codons_for_amino_acid(amino_acid);
                let aa_total: u32 = aa_codons.iter()
                    .map(|codon| codon_counts.get(*codon).copied().unwrap_or(0))
                    .sum();
                
                if aa_total > 0 {
                    // Calculate relative synonymous codon usage (RSCU)
                    let expected_freq = 1.0 / aa_codons.len() as f64;
                    for codon in aa_codons {
                        let observed_freq = codon_counts.get(codon).copied().unwrap_or(0) as f64 / aa_total as f64;
                        let rscu = if expected_freq > 0.0 { observed_freq / expected_freq } else { 0.0 };
                        features.push(rscu);
                    }
                } else {
                    // No codons for this amino acid
                    features.extend(vec![0.0; aa_codons.len()]);
                }
            }
        } else {
            // No codons found
            features.extend(vec![0.0; 61]); // 61 sense codons
        }
        
        Ok(features)
    }
    
    /// Pattern recognition features
    fn extract_pattern_features(&self, sequence: &str) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        // Tandem repeats
        features.push(self.pattern_recognizers.detect_tandem_repeats(sequence));
        
        // Inverted repeats (palindromes)
        features.push(self.pattern_recognizers.detect_inverted_repeats(sequence));
        
        // Direct repeats
        features.push(self.pattern_recognizers.detect_direct_repeats(sequence));
        
        // CpG islands
        features.push(self.pattern_recognizers.detect_cpg_islands(sequence));
        
        // Low complexity regions
        features.push(self.pattern_recognizers.detect_low_complexity_regions(sequence));
        
        // Homopolymer runs
        features.extend(self.pattern_recognizers.analyze_homopolymer_runs(sequence));
        
        // Periodic patterns
        features.extend(self.pattern_recognizers.detect_periodic_patterns(sequence));
        
        Ok(features)
    }
    
    /// Sequence complexity features
    fn extract_complexity_features(&self, sequence: &str) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        // Shannon entropy
        features.push(calculate_sequence_complexity(sequence));
        
        // Linguistic complexity (number of unique substrings)
        features.push(self.calculate_linguistic_complexity(sequence));
        
        // Compressibility (estimate using simple compression)
        features.push(self.estimate_compressibility(sequence));
        
        // k-mer diversity for various k
        for k in &[3, 4, 5, 6] {
            features.push(self.calculate_kmer_diversity(sequence, *k));
        }
        
        // Fractal dimension approximation
        features.push(self.estimate_fractal_dimension(sequence));
        
        Ok(features)
    }
    
    /// k-mer based features using various signatures
    fn extract_kmer_features(&self, sequence: &str) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        for &k in &self.config.kmer_sizes {
            // Extract k-mers
            let kmers = self.extract_kmers(sequence, k);
            
            // k-mer frequency spectrum
            let freq_spectrum = self.calculate_kmer_frequency_spectrum(&kmers);
            features.extend(freq_spectrum);
            
            // k-mer signature (hash-based)
            let signature = self.kmer_signatures.calculate_signature(&kmers, k)?;
            features.extend(signature);
        }
        
        // Limit to configured dimension
        features.truncate(self.config.kmer_feature_dim);
        features.resize(self.config.kmer_feature_dim, 0.0);
        
        Ok(features)
    }
    
    /// Graph topology features
    fn extract_topology_features(&self, graph: &AssemblyGraph) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        // Basic graph statistics
        features.push(graph.graph_fragment.nodes.len() as f64);
        features.push(graph.graph_fragment.edges.len() as f64);
        
        // Degree distribution
        let degrees = self.calculate_degree_distribution(graph);
        features.push(degrees.mean);
        features.push(degrees.variance);
        features.push(degrees.max as f64);
        features.push(degrees.min as f64);
        
        // Connectivity
        features.push(self.calculate_average_path_length(graph));
        features.push(self.calculate_diameter(graph));
        features.push(self.count_connected_components(graph) as f64);
        
        // Graph density
        let n = graph.graph_fragment.nodes.len() as f64;
        let m = graph.graph_fragment.edges.len() as f64;
        let density = if n > 1.0 { 2.0 * m / (n * (n - 1.0)) } else { 0.0 };
        features.push(density);
        
        Ok(features)
    }
    
    /// Centrality measures for nodes
    fn extract_centrality_features(&self, graph: &AssemblyGraph) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        // Calculate betweenness centrality statistics
        let betweenness_values = self.calculate_betweenness_centrality(graph)?;
        features.push(betweenness_values.iter().sum::<f64>() / betweenness_values.len() as f64); // Mean
        features.push(self.calculate_variance(&betweenness_values)); // Variance
        
        // Calculate closeness centrality statistics
        let closeness_values = self.calculate_closeness_centrality(graph)?;
        features.push(closeness_values.iter().sum::<f64>() / closeness_values.len() as f64); // Mean
        features.push(self.calculate_variance(&closeness_values)); // Variance
        
        // Calculate eigenvector centrality statistics
        let eigenvector_values = self.calculate_eigenvector_centrality(graph)?;
        features.push(eigenvector_values.iter().sum::<f64>() / eigenvector_values.len() as f64); // Mean
        features.push(self.calculate_variance(&eigenvector_values)); // Variance
        
        Ok(features)
    }
    
    /// Clustering features
    fn extract_clustering_features(&self, graph: &AssemblyGraph) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        // Global clustering coefficient
        features.push(self.calculate_global_clustering_coefficient(graph));
        
        // Average local clustering coefficient
        features.push(self.calculate_average_local_clustering_coefficient(graph));
        
        // Transitivity
        features.push(self.calculate_transitivity(graph));
        
        // Modularity (simplified)
        features.push(self.estimate_modularity(graph));
        
        Ok(features)
    }
    
    /// Local topology features for a specific node
    fn extract_local_topology_features(&self, node: &GraphNode, graph: &AssemblyGraph) -> Result<Vec<f64>> {
        let mut features = Vec::new();
        
        // Calculate node degree
        let (in_degree, out_degree) = self.calculate_node_degrees(node.kmer.hash, graph);
        features.push(in_degree as f64);
        features.push(out_degree as f64);
        features.push((in_degree + out_degree) as f64); // Total degree
        
        // Local clustering coefficient
        features.push(self.calculate_local_clustering_coefficient(node.kmer.hash, graph));
        
        // Number of triangles involving this node
        features.push(self.count_triangles_at_node(node.kmer.hash, graph) as f64);
        
        // Average neighbor degree
        features.push(self.calculate_average_neighbor_degree(node.kmer.hash, graph));
        
        Ok(features)
    }
    
    // Helper methods for feature calculation
    
    fn calculate_linguistic_complexity(&self, sequence: &str) -> f64 {
        let mut substrings = AHashSet::new();
        let k = 6; // Use 6-mers for linguistic complexity
        
        for window in sequence.as_bytes().windows(k) {
            if let Ok(substr) = std::str::from_utf8(window) {
                substrings.insert(substr.to_string());
            }
        }
        
        if sequence.len() >= k {
            substrings.len() as f64 / (sequence.len() - k + 1) as f64
        } else {
            0.0
        }
    }
    
    fn estimate_compressibility(&self, sequence: &str) -> f64 {
        // Simple run-length encoding estimate
        let mut compressed_length = 0;
        let mut chars = sequence.chars().peekable();
        
        while let Some(current_char) = chars.next() {
            let mut count = 1;
            while let Some(&next_char) = chars.peek() {
                if next_char == current_char {
                    count += 1;
                    chars.next();
                } else {
                    break;
                }
            }
            compressed_length += if count > 1 { 2 } else { 1 }; // Character + count or just character
        }
        
        if sequence.len() > 0 {
            compressed_length as f64 / sequence.len() as f64
        } else {
            1.0
        }
    }
    
    fn calculate_kmer_diversity(&self, sequence: &str, k: usize) -> f64 {
        let mut unique_kmers = AHashSet::new();
        let mut total_kmers = 0;
        
        for window in sequence.as_bytes().windows(k) {
            if let Ok(kmer) = std::str::from_utf8(window) {
                unique_kmers.insert(kmer.to_string());
                total_kmers += 1;
            }
        }
        
        if total_kmers > 0 {
            unique_kmers.len() as f64 / total_kmers as f64
        } else {
            0.0
        }
    }
    
    fn estimate_fractal_dimension(&self, sequence: &str) -> f64 {
        // Box-counting method approximation
        let mut dimension = 0.0;
        let max_box_size = sequence.len().min(100);
        
        for box_size in (1..=max_box_size).step_by(5) {
            let boxes_needed = (sequence.len() + box_size - 1) / box_size;
            if boxes_needed > 1 {
                dimension = (boxes_needed as f64).log2() / (1.0 / box_size as f64).log2();
            }
        }
        
        dimension.min(2.0) // Cap at 2D
    }
    
    fn extract_kmers(&self, sequence: &str, k: usize) -> Vec<String> {
        sequence.as_bytes()
            .windows(k)
            .filter_map(|window| std::str::from_utf8(window).ok())
            .map(|s| s.to_string())
            .collect()
    }
    
    fn calculate_kmer_frequency_spectrum(&self, kmers: &[String]) -> Vec<f64> {
        let mut counts = AHashMap::new();
        for kmer in kmers {
            *counts.entry(kmer.clone()).or_insert(0) += 1;
        }
        
        // Create frequency spectrum (frequency of frequencies)
        let mut spectrum = vec![0.0; 10]; // Track frequencies 1-10+
        for &count in counts.values() {
            let bin = (count - 1).min(9) as usize;
            spectrum[bin] += 1.0;
        }
        
        // Normalize
        let total: f64 = spectrum.iter().sum();
        if total > 0.0 {
            for freq in &mut spectrum {
                *freq /= total;
            }
        }
        
        spectrum
    }
    
    fn calculate_degree_distribution(&self, graph: &AssemblyGraph) -> DegreeStats {
        let mut degrees = Vec::new();
        
        for &node_hash in graph.graph_fragment.nodes.keys() {
            let (in_deg, out_deg) = self.calculate_node_degrees(node_hash, graph);
            degrees.push(in_deg + out_deg);
        }
        
        if degrees.is_empty() {
            return DegreeStats { mean: 0.0, variance: 0.0, min: 0, max: 0 };
        }
        
        let mean = degrees.iter().sum::<usize>() as f64 / degrees.len() as f64;
        let variance = degrees.iter()
            .map(|&d| (d as f64 - mean).powi(2))
            .sum::<f64>() / degrees.len() as f64;
        
        DegreeStats {
            mean,
            variance,
            min: *degrees.iter().min().unwrap(),
            max: *degrees.iter().max().unwrap(),
        }
    }
    
    fn calculate_node_degrees(&self, node_hash: u64, graph: &AssemblyGraph) -> (usize, usize) {
        let mut in_degree = 0;
        let mut out_degree = 0;
        
        for edge in &graph.graph_fragment.edges {
            if edge.from_hash == node_hash {
                out_degree += 1;
            }
            if edge.to_hash == node_hash {
                in_degree += 1;
            }
        }
        
        (in_degree, out_degree)
    }
    
    fn calculate_local_clustering_coefficient(&self, node_hash: u64, graph: &AssemblyGraph) -> f64 {
        // Get neighbors
        let mut neighbors = AHashSet::new();
        for edge in &graph.graph_fragment.edges {
            if edge.from_hash == node_hash {
                neighbors.insert(edge.to_hash);
            }
            if edge.to_hash == node_hash {
                neighbors.insert(edge.from_hash);
            }
        }
        
        if neighbors.len() < 2 {
            return 0.0;
        }
        
        // Count edges between neighbors
        let mut neighbor_edges = 0;
        for edge in &graph.graph_fragment.edges {
            if neighbors.contains(&edge.from_hash) && neighbors.contains(&edge.to_hash) {
                neighbor_edges += 1;
            }
        }
        
        let possible_edges = neighbors.len() * (neighbors.len() - 1) / 2;
        if possible_edges > 0 {
            neighbor_edges as f64 / possible_edges as f64
        } else {
            0.0
        }
    }
    
    fn combine_features(&self, seq_features: &Array1<f64>, graph_features: &Array1<f64>, kmer_features: &Array1<f64>) -> Result<Array1<f64>> {
        let total_dim = seq_features.len() + graph_features.len() + kmer_features.len();
        let mut combined = Vec::with_capacity(total_dim);
        
        combined.extend(seq_features.iter());
        combined.extend(graph_features.iter());
        combined.extend(kmer_features.iter());
        
        Ok(Array1::from_vec(combined))
    }
    
    fn normalize_feature_dimensions(&self, features: &mut Vec<f64>, target_dim: usize) {
        if features.len() < target_dim {
            features.resize(target_dim, 0.0);
        } else if features.len() > target_dim {
            features.truncate(target_dim);
        }
    }
    
    // Feature name getters for interpretability
    
    fn get_composition_feature_names(&self) -> Vec<String> {
        let mut names = Vec::new();
        
        // Single nucleotides
        for nuc in &["A", "C", "G", "T"] {
            names.push(format!("freq_{}", nuc));
        }
        
        // Dinucleotides
        for nuc1 in &["A", "C", "G", "T"] {
            for nuc2 in &["A", "C", "G", "T"] {
                names.push(format!("dinuc_{}_{}", nuc1, nuc2));
            }
        }
        
        names
    }
    
    fn get_codon_feature_names(&self) -> Vec<String> {
        // Simplified - would include all codon names
        (0..61).map(|i| format!("codon_{}", i)).collect()
    }
    
    fn get_pattern_feature_names(&self) -> Vec<String> {
        vec![
            "tandem_repeats".to_string(),
            "inverted_repeats".to_string(),
            "direct_repeats".to_string(),
            "cpg_islands".to_string(),
            "low_complexity".to_string(),
            "homopolymer_a".to_string(),
            "homopolymer_c".to_string(),
            "homopolymer_g".to_string(),
            "homopolymer_t".to_string(),
            "periodic_2".to_string(),
            "periodic_3".to_string(),
            "periodic_4".to_string(),
        ]
    }
    
    fn get_complexity_feature_names(&self) -> Vec<String> {
        let mut names = vec![
            "shannon_entropy".to_string(),
            "linguistic_complexity".to_string(),
            "compressibility".to_string(),
        ];
        
        for k in &[3, 4, 5, 6] {
            names.push(format!("kmer_diversity_{}", k));
        }
        
        names.push("fractal_dimension".to_string());
        names
    }
    
    fn get_kmer_feature_names(&self) -> Vec<String> {
        let mut names = Vec::new();
        
        for &k in &self.config.kmer_sizes {
            for i in 0..10 {
                names.push(format!("kmer_{}_freq_{}", k, i));
            }
            for i in 0..16 {
                names.push(format!("kmer_{}_sig_{}", k, i));
            }
        }
        
        names
    }
    
    // Placeholder implementations for complex calculations
    
    fn calculate_average_path_length(&self, _graph: &AssemblyGraph) -> f64 {
        // Simplified implementation
        5.0 // Would calculate actual average shortest path length
    }
    
    fn calculate_diameter(&self, _graph: &AssemblyGraph) -> f64 {
        // Simplified implementation
        10.0 // Would calculate actual graph diameter
    }
    
    fn count_connected_components(&self, _graph: &AssemblyGraph) -> usize {
        // Simplified implementation
        1 // Would calculate actual connected components
    }
    
    fn calculate_betweenness_centrality(&self, graph: &AssemblyGraph) -> Result<Vec<f64>> {
        // Simplified implementation
        Ok(vec![0.5; graph.graph_fragment.nodes.len()])
    }
    
    fn calculate_closeness_centrality(&self, graph: &AssemblyGraph) -> Result<Vec<f64>> {
        // Simplified implementation
        Ok(vec![0.5; graph.graph_fragment.nodes.len()])
    }
    
    fn calculate_eigenvector_centrality(&self, graph: &AssemblyGraph) -> Result<Vec<f64>> {
        // Simplified implementation
        Ok(vec![0.5; graph.graph_fragment.nodes.len()])
    }
    
    fn calculate_global_clustering_coefficient(&self, _graph: &AssemblyGraph) -> f64 {
        0.3 // Simplified
    }
    
    fn calculate_average_local_clustering_coefficient(&self, _graph: &AssemblyGraph) -> f64 {
        0.4 // Simplified
    }
    
    fn calculate_transitivity(&self, _graph: &AssemblyGraph) -> f64 {
        0.35 // Simplified
    }
    
    fn estimate_modularity(&self, _graph: &AssemblyGraph) -> f64 {
        0.6 // Simplified
    }
    
    fn count_triangles_at_node(&self, _node_hash: u64, _graph: &AssemblyGraph) -> usize {
        2 // Simplified
    }
    
    fn calculate_average_neighbor_degree(&self, _node_hash: u64, _graph: &AssemblyGraph) -> f64 {
        3.5 // Simplified
    }
    
    fn calculate_variance(&self, values: &[f64]) -> f64 {
        if values.is_empty() {
            return 0.0;
        }
        
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        values.iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f64>() / values.len() as f64
    }
}

// Supporting data structures

#[derive(Debug, Clone)]
struct DegreeStats {
    mean: f64,
    variance: f64,
    min: usize,
    max: usize,
}

/// Codon usage tables for different organisms
pub struct CodonUsageTables {
    standard_table: AHashMap<String, String>, // Codon -> Amino acid
    amino_acid_codons: AHashMap<String, Vec<&'static str>>, // Amino acid -> Codons
}

impl CodonUsageTables {
    fn load_standard_tables() -> Result<Self> {
        let mut standard_table = AHashMap::new();
        let mut amino_acid_codons = AHashMap::new();
        
        // Standard genetic code
        let codon_map = [
            ("TTT", "F"), ("TTC", "F"), ("TTA", "L"), ("TTG", "L"),
            ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"),
            ("TAT", "Y"), ("TAC", "Y"), ("TAA", "*"), ("TAG", "*"),
            ("TGT", "C"), ("TGC", "C"), ("TGA", "*"), ("TGG", "W"),
            // ... (would include all 64 codons)
        ];
        
        for (codon, aa) in &codon_map {
            standard_table.insert(codon.to_string(), aa.to_string());
            amino_acid_codons.entry(aa.to_string())
                .or_insert_with(Vec::new)
                .push(codon);
        }
        
        Ok(Self {
            standard_table,
            amino_acid_codons,
        })
    }
    
    fn get_codons_for_amino_acid(&self, amino_acid: &str) -> Vec<&str> {
        self.amino_acid_codons.get(amino_acid)
            .map(|codons| codons.iter().copied().collect())
            .unwrap_or_default()
    }
}

/// k-mer signature calculator
pub struct KmerSignatures {
    hash_functions: Vec<RandomHashFunction>,
}

impl KmerSignatures {
    fn new(kmer_sizes: &[usize]) -> Result<Self> {
        let mut hash_functions = Vec::new();
        
        for _ in 0..16 { // 16 hash functions for signatures
            hash_functions.push(RandomHashFunction::new());
        }
        
        Ok(Self { hash_functions })
    }
    
    fn calculate_signature(&self, kmers: &[String], _k: usize) -> Result<Vec<f64>> {
        let mut signature = Vec::new();
        
        for hash_func in &self.hash_functions {
            let mut min_hash = u64::MAX;
            for kmer in kmers {
                let hash = hash_func.hash(kmer);
                min_hash = min_hash.min(hash);
            }
            signature.push((min_hash % 10000) as f64 / 10000.0); // Normalize
        }
        
        Ok(signature)
    }
}

struct RandomHashFunction {
    seed: u64,
}

impl RandomHashFunction {
    fn new() -> Self {
        Self {
            seed: fastrand::u64(..),
        }
    }
    
    fn hash(&self, input: &str) -> u64 {
        let mut hash = self.seed;
        for byte in input.bytes() {
            hash = hash.wrapping_mul(31).wrapping_add(byte as u64);
        }
        hash
    }
}

/// Pattern recognizers for biological sequences
pub struct PatternRecognizers;

impl PatternRecognizers {
    fn new() -> Result<Self> {
        Ok(Self)
    }
    
    fn detect_tandem_repeats(&self, sequence: &str) -> f64 {
        // Simplified tandem repeat detection
        let mut max_repeat_length = 0;
        
        for period in 2..=10 {
            for start in 0..sequence.len().saturating_sub(period * 2) {
                let pattern = &sequence[start..start + period];
                let mut repeat_count = 1;
                
                let mut pos = start + period;
                while pos + period <= sequence.len() {
                    if &sequence[pos..pos + period] == pattern {
                        repeat_count += 1;
                        pos += period;
                    } else {
                        break;
                    }
                }
                
                if repeat_count >= 3 {
                    max_repeat_length = max_repeat_length.max(repeat_count * period);
                }
            }
        }
        
        max_repeat_length as f64 / sequence.len() as f64
    }
    
    fn detect_inverted_repeats(&self, sequence: &str) -> f64 {
        // Simplified inverted repeat detection
        let mut max_palindrome = 0;
        
        for center in 0..sequence.len() {
            // Odd-length palindromes
            let mut length = 0;
            while center >= length && center + length < sequence.len() {
                if sequence.chars().nth(center - length) == 
                   self.complement(sequence.chars().nth(center + length).unwrap_or('N')) {
                    length += 1;
                } else {
                    break;
                }
            }
            max_palindrome = max_palindrome.max(length * 2 - 1);
            
            // Even-length palindromes
            if center + 1 < sequence.len() {
                let mut length = 0;
                while center >= length && center + 1 + length < sequence.len() {
                    if sequence.chars().nth(center - length) == 
                       self.complement(sequence.chars().nth(center + 1 + length).unwrap_or('N')) {
                        length += 1;
                    } else {
                        break;
                    }
                }
                max_palindrome = max_palindrome.max(length * 2);
            }
        }
        
        max_palindrome as f64 / sequence.len() as f64
    }
    
    fn detect_direct_repeats(&self, _sequence: &str) -> f64 {
        0.1 // Simplified
    }
    
    fn detect_cpg_islands(&self, sequence: &str) -> f64 {
        // Simplified CpG island detection
        let mut cpg_count = 0;
        let mut total_dinucs = 0;
        
        for window in sequence.as_bytes().windows(2) {
            if let Ok(dinuc) = std::str::from_utf8(window) {
                total_dinucs += 1;
                if dinuc.to_uppercase() == "CG" {
                    cpg_count += 1;
                }
            }
        }
        
        if total_dinucs > 0 {
            cpg_count as f64 / total_dinucs as f64
        } else {
            0.0
        }
    }
    
    fn detect_low_complexity_regions(&self, sequence: &str) -> f64 {
        // Use entropy-based detection
        1.0 - calculate_sequence_complexity(sequence)
    }
    
    fn analyze_homopolymer_runs(&self, sequence: &str) -> Vec<f64> {
        let mut features = vec![0.0; 4]; // A, C, G, T
        let bases = ['A', 'C', 'G', 'T'];
        
        for (i, &base) in bases.iter().enumerate() {
            let mut max_run = 0;
            let mut current_run = 0;
            
            for ch in sequence.chars() {
                if ch.to_ascii_uppercase() == base {
                    current_run += 1;
                    max_run = max_run.max(current_run);
                } else {
                    current_run = 0;
                }
            }
            
            features[i] = max_run as f64 / sequence.len() as f64;
        }
        
        features
    }
    
    fn detect_periodic_patterns(&self, sequence: &str) -> Vec<f64> {
        let mut features = Vec::new();
        
        for period in 2..=4 {
            let mut period_score = 0.0;
            let mut comparisons = 0;
            
            for i in 0..sequence.len().saturating_sub(period) {
                if let (Some(a), Some(b)) = (sequence.chars().nth(i), sequence.chars().nth(i + period)) {
                    if a == b {
                        period_score += 1.0;
                    }
                    comparisons += 1;
                }
            }
            
            if comparisons > 0 {
                features.push(period_score / comparisons as f64);
            } else {
                features.push(0.0);
            }
        }
        
        features
    }
    
    fn complement(&self, base: char) -> char {
        match base.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_sequence_feature_extraction() {
        let config = FeatureConfig::default();
        let extractor = AdvancedFeatureExtractor::new(config).unwrap();
        
        let sequence = "ATCGATCGATCGATCGATCGATCGATCG";
        let features = extractor.extract_sequence_features(sequence).unwrap();
        
        assert_eq!(features.sequence_features.len(), extractor.config.sequence_feature_dim);
        assert!(features.metadata.gc_content > 0.0);
        assert!(features.metadata.complexity_score > 0.0);
    }
    
    #[test]
    fn test_composition_features() {
        let config = FeatureConfig::default();
        let extractor = AdvancedFeatureExtractor::new(config).unwrap();
        
        let sequence = "AAACCCGGGTTT";
        let features = extractor.extract_composition_features(sequence).unwrap();
        
        assert!(!features.is_empty());
        // Should have equal frequencies for A, C, G, T
        assert!((features[0] - 0.25).abs() < 0.01); // A frequency
        assert!((features[1] - 0.25).abs() < 0.01); // C frequency
    }
    
    #[test]
    fn test_pattern_recognition() {
        let recognizers = PatternRecognizers::new().unwrap();
        
        // Test tandem repeats
        let tandem_seq = "ATCGATCGATCGATCG";
        let tandem_score = recognizers.detect_tandem_repeats(tandem_seq);
        assert!(tandem_score > 0.0);
        
        // Test palindromes
        let palindrome_seq = "ATCGATCGAT";
        let palindrome_score = recognizers.detect_inverted_repeats(palindrome_seq);
        assert!(palindrome_score > 0.0);
    }
    
    #[test]
    fn test_feature_normalization() {
        let config = FeatureConfig::default();
        let extractor = AdvancedFeatureExtractor::new(config).unwrap();
        
        let mut features = vec![1.0, 2.0, 3.0];
        extractor.normalize_feature_dimensions(&mut features, 5);
        
        assert_eq!(features.len(), 5);
        assert_eq!(features[3], 0.0); // Padded with zeros
        assert_eq!(features[4], 0.0);
    }
}