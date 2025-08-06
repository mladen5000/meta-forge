use std::collections::{HashMap, HashSet};
use anyhow::{Result, anyhow};
use serde::{Serialize, Deserialize};
use crate::core_data_structures::{GraphFragment, GraphEdge};

/// GNN-based system for resolving repeats in minimiser graphs
pub struct RepeatResolverGNN {
    /// The graph neural network model
    model: GraphNeuralNetwork,
    /// Node feature extractor
    feature_extractor: NodeFeatureExtractor,
    /// Edge scoring threshold
    confidence_threshold: f32,
    /// ONNX runtime for inference (CPU-only)
    onnx_session: Option<OnnxSession>,
}

/// Simplified ONNX session wrapper for CPU inference
pub struct OnnxSession {
    model_path: String,
    input_names: Vec<String>,
    output_names: Vec<String>,
}

/// Graph Neural Network implementation
#[derive(Clone, Serialize, Deserialize)]
pub struct GraphNeuralNetwork {
    /// Message passing layers
    message_layers: Vec<MessagePassingLayer>,
    /// Final edge classifier
    edge_classifier: EdgeClassifier,
    /// Model hyperparameters
    config: GNNConfig,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct MessagePassingLayer {
    /// Node update weights: [input_dim, hidden_dim]
    node_weights: Vec<Vec<f32>>,
    /// Message aggregation weights
    message_weights: Vec<Vec<f32>>,
    /// Bias vectors
    bias: Vec<f32>,
    /// Layer normalization parameters
    layer_norm_scale: Vec<f32>,
    layer_norm_bias: Vec<f32>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct EdgeClassifier {
    /// Final classification weights
    weights: Vec<Vec<f32>>,
    bias: Vec<f32>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct GNNConfig {
    node_feature_dim: usize,
    hidden_dim: usize,
    num_layers: usize,
    dropout_rate: f32,
    learning_rate: f32,
}

/// Node feature extraction from minimiser graph
pub struct NodeFeatureExtractor {
    feature_dim: usize,
}

/// Graph representation for GNN input
#[derive(Clone)]
pub struct MinimiserGraphGNN {
    /// Node features: node_id -> feature_vector
    pub node_features: HashMap<u64, Vec<f32>>,
    /// Edge list with features: (from, to) -> edge_features
    pub edges: HashMap<(u64, u64), Vec<f32>>,
    /// Adjacency list for message passing
    pub adjacency: HashMap<u64, Vec<u64>>,
}

/// Edge prediction result
#[derive(Clone, Debug)]
pub struct EdgePrediction {
    pub from_node: u64,
    pub to_node: u64,
    pub confidence_score: f32,
    pub is_correct: bool,
    pub repeat_type: RepeatType,
}

#[derive(Clone, Debug, PartialEq)]
pub enum RepeatType {
    TandemRepeat,
    InterspersedRepeat,
    InvertedRepeat,
    NotRepeat,
}

impl RepeatResolverGNN {
    pub fn new(model_path: Option<&str>) -> Result<Self> {
        let config = GNNConfig {
            node_feature_dim: 64,
            hidden_dim: 128,
            num_layers: 3,
            dropout_rate: 0.1,
            learning_rate: 0.001,
        };

        let model = GraphNeuralNetwork::new(config.clone());
        let feature_extractor = NodeFeatureExtractor::new(config.node_feature_dim);

        let onnx_session = if let Some(path) = model_path {
            Some(OnnxSession::load(path)?)
        } else {
            None
        };

        Ok(Self {
            model,
            feature_extractor,
            confidence_threshold: 0.7,
            onnx_session,
        })
    }

    /// Convert minimiser graph to GNN input format
    pub fn prepare_graph_input(&self, graph: &GraphFragment) -> Result<MinimiserGraphGNN> {
        let mut gnn_graph = MinimiserGraphGNN {
            node_features: HashMap::new(),
            edges: HashMap::new(),
            adjacency: HashMap::new(),
        };

        // Extract node features from the graph nodes
        for (&node_id, _node) in &graph.nodes {
            let features = self.feature_extractor.extract_node_features(
                node_id, 
                &graph.edges, 
                graph
            )?;
            gnn_graph.node_features.insert(node_id, features);
            
            // Build adjacency list from edges
            let mut adj_list = Vec::new();
            for edge in &graph.edges {
                if edge.from_hash == node_id {
                    adj_list.push(edge.to_hash);
                } else if edge.to_hash == node_id {
                    adj_list.push(edge.from_hash);
                }
            }
            gnn_graph.adjacency.insert(node_id, adj_list);
        }

        // Extract edge features
        for edge in &graph.edges {
            let edge_features = self.feature_extractor.extract_edge_features(
                edge.from_hash, 
                edge.to_hash, 
                edge.weight, 
                graph
            )?;
            gnn_graph.edges.insert((edge.from_hash, edge.to_hash), edge_features);
        }

        Ok(gnn_graph)
    }

    /// Predict edge correctness using GNN
    pub fn predict_edge_correctness(&mut self, graph: &MinimiserGraphGNN) -> Result<Vec<EdgePrediction>> {
        if let Some(ref onnx_session) = self.onnx_session {
            // Use ONNX runtime for inference
            self.predict_with_onnx(graph, onnx_session)
        } else {
            // Use built-in PyTorch-style model
            self.predict_with_builtin(graph)
        }
    }

    fn predict_with_builtin(&self, graph: &MinimiserGraphGNN) -> Result<Vec<EdgePrediction>> {
        let mut predictions = Vec::new();
        
        // Forward pass through GNN
        let node_embeddings = self.model.forward_pass(graph)?;
        
        // Predict edge correctness
        for (&(from_node, to_node), _) in &graph.edges {
            let from_embedding = node_embeddings.get(&from_node)
                .ok_or_else(|| anyhow!("Missing embedding for node {}", from_node))?;
            let to_embedding = node_embeddings.get(&to_node)
                .ok_or_else(|| anyhow!("Missing embedding for node {}", to_node))?;
            
            let edge_score = self.model.predict_edge(from_embedding, to_embedding)?;
            let is_correct = edge_score > self.confidence_threshold;
            
            // Classify repeat type based on graph structure
            let repeat_type = self.classify_repeat_type(from_node, to_node, graph)?;
            
            predictions.push(EdgePrediction {
                from_node,
                to_node,
                confidence_score: edge_score,
                is_correct,
                repeat_type,
            });
        }

        Ok(predictions)
    }

    fn predict_with_onnx(&self, graph: &MinimiserGraphGNN, session: &OnnxSession) -> Result<Vec<EdgePrediction>> {
        // Convert graph to ONNX input tensors
        let (node_tensor, edge_tensor, adjacency_tensor) = self.graph_to_tensors(graph)?;
        
        // Run ONNX inference
        let outputs = session.run(&[
            ("node_features", node_tensor),
            ("edge_features", edge_tensor),
            ("adjacency", adjacency_tensor),
        ])?;
        
        // Parse output tensors back to predictions
        let edge_scores = outputs.get("edge_scores")
            .ok_or_else(|| anyhow!("Missing edge_scores output"))?;
        
        self.tensors_to_predictions(graph, edge_scores)
    }

    fn classify_repeat_type(&self, from_node: u64, to_node: u64, graph: &MinimiserGraphGNN) -> Result<RepeatType> {
        // Analyze local graph structure to classify repeat type
        
        let from_neighbors = graph.adjacency.get(&from_node).map(|v| v.len()).unwrap_or(0);
        let to_neighbors = graph.adjacency.get(&to_node).map(|v| v.len()).unwrap_or(0);
        
        // Simple heuristics - in practice, use learned features
        if from_neighbors > 5 && to_neighbors > 5 {
            // High-degree nodes suggest repeats
            if self.check_tandem_pattern(from_node, to_node, graph)? {
                Ok(RepeatType::TandemRepeat)
            } else if self.check_inverted_pattern(from_node, to_node, graph)? {
                Ok(RepeatType::InvertedRepeat)
            } else {
                Ok(RepeatType::InterspersedRepeat)
            }
        } else {
            Ok(RepeatType::NotRepeat)
        }
    }

    fn check_tandem_pattern(&self, from_node: u64, to_node: u64, graph: &MinimiserGraphGNN) -> Result<bool> {
        // Check if nodes form part of tandem repeat pattern
        // Look for short cycles indicating tandem repeats
        if let Some(from_neighbors) = graph.adjacency.get(&from_node) {
            if let Some(to_neighbors) = graph.adjacency.get(&to_node) {
                // Check for common neighbors (triangle pattern)
                let common_neighbors: HashSet<_> = from_neighbors.iter()
                    .filter(|&&n| to_neighbors.contains(&n))
                    .collect();
                
                // Tandem repeats often create triangular patterns
                Ok(common_neighbors.len() > 2)
            } else {
                Ok(false)
            }
        } else {
            Ok(false)
        }
    }

    fn check_inverted_pattern(&self, from_node: u64, to_node: u64, graph: &MinimiserGraphGNN) -> Result<bool> {
        // Check for inverted repeat patterns
        // Look for palindromic-like structures in the graph
        
        // Simple heuristic: check if reverse complement nodes exist
        // (This would require actual k-mer sequence analysis in practice)
        let from_neighbors = graph.adjacency.get(&from_node).map(|v| v.len()).unwrap_or(0);
        let to_neighbors = graph.adjacency.get(&to_node).map(|v| v.len()).unwrap_or(0);
        
        // Inverted repeats often have symmetric degree patterns
        Ok((from_neighbors == to_neighbors) && from_neighbors > 3)
    }

    /// Apply GNN predictions to resolve repeats in assembly
    pub fn resolve_repeats(&self, 
        graph: &mut GraphFragment, 
        predictions: &[EdgePrediction]
    ) -> Result<ResolvedGraph> {
        let mut resolved = ResolvedGraph {
            original_edges: graph.edges.clone(),
            resolved_edges: graph.edges.clone(),
            removed_edges: Vec::new(),
            repeat_regions: Vec::new(),
        };

        // Group predictions by confidence and repeat type
        let mut high_confidence_correct = Vec::new();
        let mut low_confidence_edges = Vec::new();
        let mut repeat_edges = Vec::new();

        for pred in predictions {
            if pred.confidence_score > self.confidence_threshold && pred.is_correct {
                high_confidence_correct.push(pred);
            } else if pred.repeat_type != RepeatType::NotRepeat {
                repeat_edges.push(pred);
            } else {
                low_confidence_edges.push(pred);
            }
        }

        // Remove low-confidence edges
        for pred in &low_confidence_edges {
            if !pred.is_correct {
                // Remove edges from both graph and resolved edges
                graph.edges.retain(|edge| {
                    !(edge.from_hash == pred.from_node && edge.to_hash == pred.to_node)
                });
                resolved.resolved_edges.retain(|edge| {
                    !(edge.from_hash == pred.from_node && edge.to_hash == pred.to_node)
                });
                resolved.removed_edges.push((pred.from_node, pred.to_node));
            }
        }

        // Handle repeat edges specially
        for pred in &repeat_edges {
            let repeat_region = self.analyze_repeat_region(pred, graph)?;
            resolved.repeat_regions.push(repeat_region);
        }

        // Edges have been updated in-place during removal

        Ok(resolved)
    }

    fn analyze_repeat_region(&self, pred: &EdgePrediction, graph: &GraphFragment) -> Result<RepeatRegion> {
        Ok(RepeatRegion {
            repeat_type: pred.repeat_type.clone(),
            nodes: vec![pred.from_node, pred.to_node],
            confidence: pred.confidence_score,
            resolution_strategy: match pred.repeat_type {
                RepeatType::TandemRepeat => ResolutionStrategy::UnrollTandem,
                RepeatType::InterspersedRepeat => ResolutionStrategy::SkipRepeat, 
                RepeatType::InvertedRepeat => ResolutionStrategy::ResolveByPairing,
                RepeatType::NotRepeat => ResolutionStrategy::KeepOriginal,
            },
        })
    }

    fn graph_to_tensors(&self, graph: &MinimiserGraphGNN) -> Result<(Vec<f32>, Vec<f32>, Vec<f32>)> {
        // Convert graph to tensor format for ONNX
        let num_nodes = graph.node_features.len();
        let num_edges = graph.edges.len();
        let feature_dim = self.model.config.node_feature_dim;

        // Node features tensor [num_nodes, feature_dim]
        let mut node_tensor = vec![0.0; num_nodes * feature_dim];
        for (i, (_, features)) in graph.node_features.iter().enumerate() {
            for (j, &feature) in features.iter().enumerate() {
                node_tensor[i * feature_dim + j] = feature;
            }
        }

        // Edge features tensor [num_edges, edge_feature_dim]
        let edge_feature_dim = 8; // Simplified edge features
        let mut edge_tensor = vec![0.0; num_edges * edge_feature_dim];
        for (i, (_, features)) in graph.edges.iter().enumerate() {
            for (j, &feature) in features.iter().take(edge_feature_dim).enumerate() {
                edge_tensor[i * edge_feature_dim + j] = feature;
            }
        }

        // Adjacency tensor as edge list [2, num_edges]
        let mut adjacency_tensor = vec![0f32; 2 * num_edges];
        let node_to_idx: HashMap<u64, usize> = graph.node_features.keys()
            .enumerate()
            .map(|(i, &k)| (k, i))
            .collect();

        for (i, &(from, to)) in graph.edges.keys().enumerate() {
            adjacency_tensor[i] = *node_to_idx.get(&from).unwrap_or(&0) as f32;
            adjacency_tensor[num_edges + i] = *node_to_idx.get(&to).unwrap_or(&0) as f32;
        }

        Ok((node_tensor, edge_tensor, adjacency_tensor))
    }

    fn tensors_to_predictions(&self, graph: &MinimiserGraphGNN, edge_scores: &[f32]) -> Result<Vec<EdgePrediction>> {
        let mut predictions = Vec::new();
        
        for (i, &(from_node, to_node)) in graph.edges.keys().enumerate() {
            if i < edge_scores.len() {
                let score = edge_scores[i];
                let is_correct = score > self.confidence_threshold;
                let repeat_type = self.classify_repeat_type(from_node, to_node, graph)?;
                
                predictions.push(EdgePrediction {
                    from_node,
                    to_node,
                    confidence_score: score,
                    is_correct,
                    repeat_type,
                });
            }
        }
        
        Ok(predictions)
    }
}

impl GraphNeuralNetwork {
    fn new(config: GNNConfig) -> Self {
        let mut message_layers = Vec::new();
        
        for layer_idx in 0..config.num_layers {
            let input_dim = if layer_idx == 0 { 
                config.node_feature_dim 
            } else { 
                config.hidden_dim 
            };
            
            message_layers.push(MessagePassingLayer::new(input_dim, config.hidden_dim));
        }

        let edge_classifier = EdgeClassifier::new(config.hidden_dim * 2); // Concatenated node embeddings

        Self {
            message_layers,
            edge_classifier,
            config,
        }
    }

    fn forward_pass(&self, graph: &MinimiserGraphGNN) -> Result<HashMap<u64, Vec<f32>>> {
        let mut node_embeddings: HashMap<u64, Vec<f32>> = graph.node_features.clone();
        
        // Message passing layers
        for layer in &self.message_layers {
            node_embeddings = layer.forward(&node_embeddings, &graph.adjacency)?;
        }

        Ok(node_embeddings)
    }

    fn predict_edge(&self, from_embedding: &[f32], to_embedding: &[f32]) -> Result<f32> {
        // Concatenate node embeddings
        let mut edge_features = Vec::with_capacity(from_embedding.len() + to_embedding.len());
        edge_features.extend_from_slice(from_embedding);
        edge_features.extend_from_slice(to_embedding);
        
        self.edge_classifier.predict(&edge_features)
    }
}

impl MessagePassingLayer {
    fn new(input_dim: usize, hidden_dim: usize) -> Self {
        let mut rng = fastrand::Rng::new();
        
        // Initialize weights with Xavier/Glorot initialization
        let node_weights = (0..hidden_dim)
            .map(|_| {
                (0..input_dim)
                    .map(|_| {
                        let limit = (6.0 / (input_dim + hidden_dim) as f32).sqrt();
                        rng.f32() * 2.0 * limit - limit
                    })
                    .collect()
            })
            .collect();

        let message_weights = (0..hidden_dim)
            .map(|_| {
                (0..hidden_dim)
                    .map(|_| {
                        let limit = (6.0 / (hidden_dim + hidden_dim) as f32).sqrt();
                        rng.f32() * 2.0 * limit - limit
                    })
                    .collect()
            })
            .collect();

        Self {
            node_weights,
            message_weights,
            bias: vec![0.0; hidden_dim],
            layer_norm_scale: vec![1.0; hidden_dim],
            layer_norm_bias: vec![0.0; hidden_dim],
        }
    }

    fn forward(
        &self,
        node_embeddings: &HashMap<u64, Vec<f32>>,
        adjacency: &HashMap<u64, Vec<u64>>,
    ) -> Result<HashMap<u64, Vec<f32>>> {
        let mut new_embeddings = HashMap::new();

        for (&node_id, current_embedding) in node_embeddings {
            // Aggregate messages from neighbors
            let mut aggregated_message = vec![0.0; self.message_weights.len()];
            
            if let Some(neighbors) = adjacency.get(&node_id) {
                for &neighbor_id in neighbors {
                    if let Some(neighbor_embedding) = node_embeddings.get(&neighbor_id) {
                        // Compute message from neighbor
                        let message = self.compute_message(neighbor_embedding)?;
                        
                        // Add to aggregated message
                        for (i, &msg_val) in message.iter().enumerate() {
                            aggregated_message[i] += msg_val;
                        }
                    }
                }
                
                // Average the messages
                let neighbor_count = neighbors.len() as f32;
                if neighbor_count > 0.0 {
                    for val in &mut aggregated_message {
                        *val /= neighbor_count;
                    }
                }
            }

            // Update node embedding
            let updated_embedding = self.update_node(current_embedding, &aggregated_message)?;
            new_embeddings.insert(node_id, updated_embedding);
        }

        Ok(new_embeddings)
    }

    fn compute_message(&self, embedding: &[f32]) -> Result<Vec<f32>> {
        let mut message = vec![0.0; self.message_weights.len()];
        
        for (i, neuron_weights) in self.message_weights.iter().enumerate() {
            let weighted_sum: f32 = embedding.iter()
                .zip(neuron_weights.iter())
                .map(|(x, w)| x * w)
                .sum();
            message[i] = Self::relu(weighted_sum);
        }
        
        Ok(message)
    }

    fn update_node(&self, current_embedding: &[f32], message: &[f32]) -> Result<Vec<f32>> {
        let mut updated = vec![0.0; self.node_weights.len()];
        
        // Transform current embedding
        for (i, neuron_weights) in self.node_weights.iter().enumerate() {
            let weighted_sum: f32 = current_embedding.iter()
                .zip(neuron_weights.iter())
                .map(|(x, w)| x * w)
                .sum();
            
            // Add message and bias
            let combined = weighted_sum + message.get(i).unwrap_or(&0.0) + self.bias[i];
            updated[i] = Self::relu(combined);
        }

        // Apply layer normalization
        let updated = self.layer_normalize(&updated)?;
        
        Ok(updated)
    }

    fn layer_normalize(&self, input: &[f32]) -> Result<Vec<f32>> {
        // Compute mean and variance
        let mean = input.iter().sum::<f32>() / input.len() as f32;
        let variance = input.iter()
            .map(|x| (x - mean).powi(2))
            .sum::<f32>() / input.len() as f32;
        
        let std_dev = (variance + 1e-5).sqrt(); // Add epsilon for numerical stability
        
        // Normalize and scale
        let normalized: Vec<f32> = input.iter()
            .zip(self.layer_norm_scale.iter())
            .zip(self.layer_norm_bias.iter())
            .map(|((x, scale), bias)| {
                ((x - mean) / std_dev) * scale + bias
            })
            .collect();
        
        Ok(normalized)
    }

    fn relu(x: f32) -> f32 {
        x.max(0.0)
    }
}

impl EdgeClassifier {
    fn new(input_dim: usize) -> Self {
        let mut rng = fastrand::Rng::new();
        
        // Single layer classifier for edge prediction
        let weights = vec![
            (0..input_dim)
                .map(|_| rng.f32() * 2.0 - 1.0)
                .collect()
        ];
        
        Self {
            weights,
            bias: vec![0.0],
        }
    }

    fn predict(&self, edge_features: &[f32]) -> Result<f32> {
        if edge_features.len() != self.weights[0].len() {
            return Err(anyhow!("Input dimension mismatch"));
        }

        let weighted_sum: f32 = edge_features.iter()
            .zip(self.weights[0].iter())
            .map(|(x, w)| x * w)
            .sum();
        
        let output = weighted_sum + self.bias[0];
        
        // Apply sigmoid activation for binary classification
        Ok(Self::sigmoid(output))
    }

    fn sigmoid(x: f32) -> f32 {
        1.0 / (1.0 + (-x).exp())
    }
}

impl NodeFeatureExtractor {
    fn new(feature_dim: usize) -> Self {
        Self { feature_dim }
    }

    fn extract_node_features(
        &self,
        node_id: u64,
        edges: &[GraphEdge],
        graph: &GraphFragment,
    ) -> Result<Vec<f32>> {
        let mut features = vec![0.0; self.feature_dim];
        
        // Count neighbors and collect weights from edges
        let mut neighbors = HashMap::new();
        for edge in edges {
            if edge.from_hash == node_id {
                neighbors.insert(edge.to_hash, edge.weight);
            } else if edge.to_hash == node_id {
                neighbors.insert(edge.from_hash, edge.weight);
            }
        }
        
        // Feature 1-4: Basic graph properties
        features[0] = neighbors.len() as f32; // Degree
        features[1] = neighbors.values().sum::<u32>() as f32; // Total weight
        features[2] = neighbors.values().map(|&w| w as f32).sum::<f32>() / neighbors.len() as f32; // Avg weight
        features[3] = (node_id % 1000) as f32 / 1000.0; // Hash-based feature
        
        // Feature 5-8: Local clustering coefficient
        let clustering = self.compute_clustering_coefficient(node_id, graph)?;
        features[4] = clustering;
        features[5] = if clustering > 0.5 { 1.0 } else { 0.0 }; // High clustering indicator
        
        // Feature 6-10: K-mer composition features (simplified)
        let composition = self.extract_composition_features(node_id)?;
        for (i, &comp) in composition.iter().take(5).enumerate() {
            if 6 + i < features.len() {
                features[6 + i] = comp;
            }
        }
        
        // Feature 11-15: Repeat indicators
        features[10] = if neighbors.len() > 10 { 1.0 } else { 0.0 }; // High degree = potential repeat
        features[11] = self.count_palindromic_neighbors(&neighbors) as f32;
        
        // Pad remaining features with normalized hash values
        for i in 12..self.feature_dim {
            features[i] = ((node_id.wrapping_mul(i as u64 + 1)) % 1000) as f32 / 1000.0;
        }
        
        Ok(features)
    }

    fn extract_edge_features(
        &self,
        from_node: u64,
        to_node: u64,
        weight: u32,
        graph: &GraphFragment,
    ) -> Result<Vec<f32>> {
        let mut features = vec![0.0; 8]; // Simplified edge features
        
        features[0] = weight as f32;
        features[1] = (weight as f32).ln(); // Log weight
        
        // Degree features
        let from_degree = graph.edges.iter()
            .filter(|e| e.from_hash == from_node || e.to_hash == from_node)
            .count() as f32;
        let to_degree = graph.edges.iter()
            .filter(|e| e.from_hash == to_node || e.to_hash == to_node)
            .count() as f32;
        features[2] = from_degree;
        features[3] = to_degree;
        features[4] = (from_degree - to_degree).abs();
        
        // Structural features
        features[5] = if self.nodes_share_neighbors(from_node, to_node, graph)? { 1.0 } else { 0.0 };
        features[6] = ((from_node ^ to_node) % 1000) as f32 / 1000.0; // Interaction hash
        features[7] = if weight > 5 { 1.0 } else { 0.0 }; // High weight indicator
        
        Ok(features)
    }

    fn compute_clustering_coefficient(&self, node_id: u64, graph: &GraphFragment) -> Result<f32> {
        // Get neighbors from edges
        let neighbor_ids: Vec<u64> = graph.edges.iter()
            .filter_map(|e| {
                if e.from_hash == node_id {
                    Some(e.to_hash)
                } else if e.to_hash == node_id {
                    Some(e.from_hash)
                } else {
                    None
                }
            })
            .collect();
            
        if neighbor_ids.len() < 2 {
            return Ok(0.0);
        }
        
        let mut triangles = 0;
        let possible_triangles = neighbor_ids.len() * (neighbor_ids.len() - 1) / 2;
            
        // Count triangles
        for i in 0..neighbor_ids.len() {
            for j in i+1..neighbor_ids.len() {
                let node_a = neighbor_ids[i];
                let node_b = neighbor_ids[j];
                
                // Check if node_a and node_b are connected
                let are_connected = graph.edges.iter().any(|e| {
                    (e.from_hash == node_a && e.to_hash == node_b) ||
                    (e.from_hash == node_b && e.to_hash == node_a)
                });
                
                if are_connected {
                    triangles += 1;
                }
            }
        }
        
        Ok(triangles as f32 / possible_triangles as f32)
    }

    fn extract_composition_features(&self, node_id: u64) -> Result<Vec<f32>> {
        // Simplified composition based on hash (in practice, use actual k-mer)
        let mut composition = vec![0.0; 4];
        
        for i in 0..4 {
            let shift = i * 16;
            composition[i] = ((node_id >> shift) & 0xFFFF) as f32 / 65535.0;
        }
        
        Ok(composition)
    }

    fn count_palindromic_neighbors(&self, neighbors: &HashMap<u64, u32>) -> usize {
        // Simplified palindrome detection (would use actual sequences in practice)
        neighbors.keys()
            .filter(|&&neighbor| {
                // Simple heuristic: even hash values as "palindromic"
                neighbor % 2 == 0
            })
            .count()
    }

    fn nodes_share_neighbors(&self, node_a: u64, node_b: u64, graph: &GraphFragment) -> Result<bool> {
        // Get neighbors for both nodes
        let neighbors_a: HashSet<u64> = graph.edges.iter()
            .filter_map(|e| {
                if e.from_hash == node_a { Some(e.to_hash) }
                else if e.to_hash == node_a { Some(e.from_hash) }
                else { None }
            })
            .collect();
            
        let neighbors_b: HashSet<u64> = graph.edges.iter()
            .filter_map(|e| {
                if e.from_hash == node_b { Some(e.to_hash) }
                else if e.to_hash == node_b { Some(e.from_hash) }
                else { None }
            })
            .collect();
            
        let common_neighbors: HashSet<_> = neighbors_a.intersection(&neighbors_b).collect();
        Ok(common_neighbors.len() > 0)
    }
}

impl OnnxSession {
    fn load(model_path: &str) -> Result<Self> {
        // In a real implementation, you'd load the ONNX model
        // For now, return a placeholder
        Ok(Self {
            model_path: model_path.to_string(),
            input_names: vec!["node_features".to_string(), "edge_features".to_string(), "adjacency".to_string()],
            output_names: vec!["edge_scores".to_string()],
        })
    }

    fn run(&self, inputs: &[(&str, Vec<f32>)]) -> Result<HashMap<String, Vec<f32>>> {
        // Placeholder ONNX inference
        // In practice, use ort crate or similar for ONNX runtime
        let mut outputs = HashMap::new();
        
        // Mock output - in reality this comes from the trained model
        let num_edges = inputs.iter()
            .find(|(name, _)| *name == "edge_features")
            .map(|(_, tensor)| tensor.len() / 8) // 8 edge features per edge
            .unwrap_or(0);
        
        let mock_scores: Vec<f32> = (0..num_edges)
            .map(|i| fastrand::f32() * 0.5 + 0.25) // Random scores between 0.25-0.75
            .collect();
        
        outputs.insert("edge_scores".to_string(), mock_scores);
        Ok(outputs)
    }
}

/// Results of repeat resolution
#[derive(Clone)]
pub struct ResolvedGraph {
    pub original_edges: Vec<GraphEdge>,
    pub resolved_edges: Vec<GraphEdge>,
    pub removed_edges: Vec<(u64, u64)>,
    pub repeat_regions: Vec<RepeatRegion>,
}

#[derive(Clone)]
pub struct RepeatRegion {
    pub repeat_type: RepeatType,
    pub nodes: Vec<u64>,
    pub confidence: f32,
    pub resolution_strategy: ResolutionStrategy,
}

#[derive(Clone, Debug)]
pub enum ResolutionStrategy {
    UnrollTandem,
    SkipRepeat,
    ResolveByPairing,
    KeepOriginal,
}

/// Integration with existing pipeline
pub fn resolve_repeats_in_pipeline(
    graph: &mut GraphFragment,
    model_path: Option<&str>,
) -> Result<ResolvedGraph> {
    let mut resolver = RepeatResolverGNN::new(model_path)?;
    
    // Convert graph to GNN format
    let gnn_graph = resolver.prepare_graph_input(graph)?;
    
    // Get edge predictions
    let predictions = resolver.predict_edge_correctness(&gnn_graph)?;
    
    // Apply resolutions
    resolver.resolve_repeats(graph, &predictions)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feature_extraction() {
        let extractor = NodeFeatureExtractor::new(32);
        let mut neighbors = HashMap::new();
        neighbors.insert(123, 5);
        neighbors.insert(456, 3);
        
        let mut graph = GraphFragment::default();
        graph.edges.insert(789, neighbors.clone());
        
        let features = extractor.extract_node_features(789, &neighbors, &graph).unwrap();
        assert_eq!(features.len(), 32);
        assert_eq!(features[0], 2.0); // Degree = 2
        assert_eq!(features[1], 8.0); // Total weight = 5 + 3
    }

    #[test]
    fn test_repeat_classification() {
        let resolver = RepeatResolverGNN::new(None).unwrap();
        let mut graph = MinimiserGraphGNN {
            node_features: HashMap::new(),
            edges: HashMap::new(),
            adjacency: HashMap::new(),
        };
        
        // Create high-degree nodes (typical of repeats)
        graph.adjacency.insert(1, vec![2, 3, 4, 5, 6, 7]);
        graph.adjacency.insert(2, vec![1, 3, 4, 5, 6, 7]);
        
        let repeat_type = resolver.classify_repeat_type(1, 2, &graph).unwrap();
        assert_ne!(repeat_type, RepeatType::NotRepeat);
    }
}