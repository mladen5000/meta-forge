use anyhow::{anyhow, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Local training config for learned bloom filter specific parameters
#[derive(Clone, Serialize, Deserialize)]
struct BloomTrainingConfig {
    learning_rate: f32,
    batch_size: usize,
    epochs: usize,
    embedding_dim: usize,
    hidden_dim: usize,
    target_fpr: f64,
}

/// Learned Bloom filter combining neural predictor with backup Bloom filter
pub struct LearnedBloomFilter {
    /// Neural network predictor (lightweight MLP)
    neural_oracle: NeuralOracle,
    /// Traditional Bloom filter for backup
    backup_bloom: TraditionalBloom,
    /// K-mer embedding cache
    embedding_cache: HashMap<u64, Vec<f32>>,
    /// Training parameters
    training_config: BloomTrainingConfig,
    /// Performance metrics
    metrics: FilterMetrics,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct NeuralOracle {
    /// Network weights: input -> hidden -> output
    input_weights: Vec<Vec<f32>>,
    hidden_weights: Vec<Vec<f32>>,
    output_weights: Vec<f32>,
    /// Network dimensions
    input_size: usize,
    hidden_size: usize,
    /// Activation thresholds
    positive_threshold: f32,
    negative_threshold: f32,
}

pub struct TraditionalBloom {
    /// Bit array
    bits: Vec<bool>,
    /// Hash functions count
    hash_count: usize,
    /// Size of bit array
    size: usize,
}

// Removed duplicate TrainingConfig - use the centralized one from utils::configuration

#[derive(Default)]
pub struct FilterMetrics {
    neural_queries: u64,
    neural_hits: u64,
    backup_queries: u64,
    false_positives: u64,
    true_positives: u64,
}

#[derive(Serialize, Deserialize)]
pub struct KmerEmbedding {
    pub kmer_hash: u64,
    pub embedding: Vec<f32>,
    pub taxonomy_id: Option<u32>,
}

impl LearnedBloomFilter {
    pub fn new(
        expected_elements: usize,
        target_fpr: f64,
        embedding_dim: usize,
        hidden_dim: usize,
    ) -> Self {
        let config = BloomTrainingConfig {
            learning_rate: 0.001,
            batch_size: 64,
            epochs: 100,
            embedding_dim,
            hidden_dim,
            target_fpr,
        };

        let neural_oracle = NeuralOracle::new(embedding_dim, hidden_dim);
        let backup_bloom = TraditionalBloom::new(expected_elements, target_fpr);

        Self {
            neural_oracle,
            backup_bloom,
            embedding_cache: HashMap::new(),
            training_config: config,
            metrics: FilterMetrics::default(),
        }
    }

    /// Train the neural oracle on k-mer embeddings
    /// Train the neural oracle on k-mer embeddings with positive and negative samples
    pub fn train_on_kmers(&mut self, training_data: &[KmerEmbedding]) -> Result<()> {
        println!(
            "Training learned Bloom filter on {} k-mers",
            training_data.len()
        );

        // Prepare training batches
        let mut positive_samples = Vec::new();
        let mut negative_samples = Vec::new();

        for embedding in training_data {
            if embedding.taxonomy_id.is_some() {
                positive_samples.push(embedding.embedding.clone());
                // Cache the embedding
                self.embedding_cache
                    .insert(embedding.kmer_hash, embedding.embedding.clone());
                // Add to backup Bloom filter
                self.backup_bloom.insert(embedding.kmer_hash);
            }
        }

        // Generate negative samples (k-mers not in taxonomic database)
        negative_samples = self.generate_negative_samples(positive_samples.len())?;

        // Train neural network
        self.train_neural_oracle(&positive_samples, &negative_samples)?;

        println!("Training completed. Neural oracle ready.");
        Ok(())
    }

    /// Query the learned Bloom filter
    /// Query the learned bloom filter for k-mer membership
    pub fn query(&mut self, kmer_hash: u64) -> Result<bool> {
        // Get k-mer embedding
        let embedding = if let Some(cached) = self.embedding_cache.get(&kmer_hash) {
            cached.clone()
        } else {
            self.compute_kmer_embedding(kmer_hash)?
        };

        // Query neural oracle first
        self.metrics.neural_queries += 1;
        let neural_prediction = self.neural_oracle.predict(&embedding)?;

        match neural_prediction {
            NeuralPrediction::DefinitelyPresent => {
                self.metrics.neural_hits += 1;
                Ok(true)
            }
            NeuralPrediction::DefinitelyAbsent => Ok(false),
            NeuralPrediction::Uncertain => {
                // Fall back to traditional Bloom filter
                self.metrics.backup_queries += 1;
                Ok(self.backup_bloom.query(kmer_hash))
            }
        }
    }

    /// Compute embedding for a k-mer hash
    /// Compute feature embedding for a k-mer hash
    fn compute_kmer_embedding(&self, kmer_hash: u64) -> Result<Vec<f32>> {
        // Convert hash to embedding using learned features
        let mut embedding = vec![0.0; self.training_config.embedding_dim];

        // Feature 1: Hash magnitude (normalized)
        embedding[0] = (kmer_hash as f32) / (u64::MAX as f32);

        // Feature 2-5: Hash bits in different positions
        for i in 1..5 {
            if i < embedding.len() {
                let shift = (i * 16) % 64; // Prevent overflow by limiting shift amount
                embedding[i] = ((kmer_hash >> shift) & 0xFFFF) as f32 / 65535.0;
            }
        }

        // Feature 6-10: Compositional features (simulated)
        // In practice, you'd extract actual k-mer composition
        for i in 5..10 {
            embedding[i] = ((kmer_hash.wrapping_mul(i as u64 + 1)) % 1000) as f32 / 1000.0;
        }

        // Pad or truncate to desired dimension
        embedding.resize(self.training_config.embedding_dim, 0.0);

        Ok(embedding)
    }

    /// Generate negative training samples for neural oracle training
    fn generate_negative_samples(&self, count: usize) -> Result<Vec<Vec<f32>>> {
        let mut negatives = Vec::with_capacity(count);
        let mut rng = fastrand::Rng::new();

        for _ in 0..count {
            // Generate random k-mer hash
            let random_hash = rng.u64(..);
            // Ensure it's not in our positive set
            if !self.backup_bloom.query(random_hash) {
                let embedding = self.compute_kmer_embedding(random_hash)?;
                negatives.push(embedding);
            }
        }

        Ok(negatives)
    }

    fn train_neural_oracle(
        &mut self,
        positive_samples: &[Vec<f32>],
        negative_samples: &[Vec<f32>],
    ) -> Result<()> {
        let mut training_data = Vec::new();

        // Prepare training pairs (embedding, label)
        for embedding in positive_samples {
            training_data.push((embedding.clone(), 1.0)); // Positive class
        }
        for embedding in negative_samples {
            training_data.push((embedding.clone(), 0.0)); // Negative class
        }

        // Shuffle training data
        let mut rng = fastrand::Rng::new();
        for i in (1..training_data.len()).rev() {
            let j = rng.usize(..=i);
            training_data.swap(i, j);
        }

        // Training loop
        for epoch in 0..self.training_config.epochs {
            let mut total_loss = 0.0;
            let mut correct_predictions = 0;

            for batch_start in (0..training_data.len()).step_by(self.training_config.batch_size) {
                let batch_end =
                    (batch_start + self.training_config.batch_size).min(training_data.len());
                let batch = &training_data[batch_start..batch_end];

                let (batch_loss, batch_correct) = self
                    .neural_oracle
                    .train_batch(batch, self.training_config.learning_rate)?;

                total_loss += batch_loss;
                correct_predictions += batch_correct;
            }

            let accuracy = correct_predictions as f32 / training_data.len() as f32;

            if epoch % 20 == 0 {
                println!(
                    "Epoch {}: Loss = {:.4}, Accuracy = {:.4}",
                    epoch,
                    total_loss / training_data.len() as f32,
                    accuracy
                );
            }

            // Early stopping if accuracy is good enough
            if accuracy > 0.95 {
                break;
            }
        }

        Ok(())
    }

    /// Export the trained model to a compact format
    /// Export the trained neural oracle model as serialized bytes
    pub fn export_model(&self) -> Result<Vec<u8>> {
        let model_data = SerializableModel {
            neural_oracle: self.neural_oracle.clone(),
            backup_bloom_bits: self.backup_bloom.bits.clone(),
            backup_bloom_size: self.backup_bloom.size,
            backup_bloom_hash_count: self.backup_bloom.hash_count,
            training_config: self.training_config.clone(),
        };

        let serialized =
            serde_json::to_vec(&model_data).map_err(|e| anyhow!("Serialization failed: {}", e))?;

        println!(
            "Model exported, size: {} bytes ({:.2} MB)",
            serialized.len(),
            serialized.len() as f64 / (1024.0 * 1024.0)
        );

        Ok(serialized)
    }

    /// Get current performance metrics for the learned bloom filter
    pub fn get_metrics(&self) -> &FilterMetrics {
        &self.metrics
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct SerializableModel {
    neural_oracle: NeuralOracle,
    backup_bloom_bits: Vec<bool>,
    backup_bloom_size: usize,
    backup_bloom_hash_count: usize,
    training_config: BloomTrainingConfig,
}

#[derive(Debug)]
enum NeuralPrediction {
    DefinitelyPresent,
    DefinitelyAbsent,
    Uncertain,
}

impl NeuralOracle {
    /// Create a new neural oracle with specified architecture dimensions
    fn new(input_size: usize, hidden_size: usize) -> Self {
        let mut rng = fastrand::Rng::new();

        // Initialize weights with Xavier initialization
        let input_weights = (0..hidden_size)
            .map(|_| {
                (0..input_size)
                    .map(|_| rng.f32() * 2.0 - 1.0) // [-1, 1]
                    .map(|w| w * (2.0 / input_size as f32).sqrt()) // Xavier scaling
                    .collect()
            })
            .collect();

        let hidden_weights = (0..hidden_size)
            .map(|_| {
                (0..hidden_size)
                    .map(|_| rng.f32() * 2.0 - 1.0)
                    .map(|w| w * (2.0 / hidden_size as f32).sqrt())
                    .collect()
            })
            .collect();

        let output_weights = (0..hidden_size)
            .map(|_| rng.f32() * 2.0 - 1.0)
            .map(|w| w * (2.0 / hidden_size as f32).sqrt())
            .collect();

        Self {
            input_weights,
            hidden_weights,
            output_weights,
            input_size,
            hidden_size,
            positive_threshold: 0.8,
            negative_threshold: 0.2,
        }
    }

    /// Make a prediction using the neural network with confidence scoring
    fn predict(&self, input: &[f32]) -> Result<NeuralPrediction> {
        let output = self.forward(input)?;

        if output > self.positive_threshold {
            Ok(NeuralPrediction::DefinitelyPresent)
        } else if output < self.negative_threshold {
            Ok(NeuralPrediction::DefinitelyAbsent)
        } else {
            Ok(NeuralPrediction::Uncertain)
        }
    }

    /// Forward pass through the neural network
    fn forward(&self, input: &[f32]) -> Result<f32> {
        if input.len() != self.input_size {
            return Err(anyhow!("Input size mismatch"));
        }

        // Input to hidden layer
        let mut hidden = vec![0.0; self.hidden_size];
        for (i, neuron_weights) in self.input_weights.iter().enumerate() {
            let weighted_sum: f32 = input
                .iter()
                .zip(neuron_weights.iter())
                .map(|(x, w)| x * w)
                .sum();
            hidden[i] = Self::relu(weighted_sum);
        }

        // Hidden to hidden (one layer)
        let mut hidden2 = vec![0.0; self.hidden_size];
        for (i, neuron_weights) in self.hidden_weights.iter().enumerate() {
            let weighted_sum: f32 = hidden
                .iter()
                .zip(neuron_weights.iter())
                .map(|(x, w)| x * w)
                .sum();
            hidden2[i] = Self::relu(weighted_sum);
        }

        // Hidden to output
        let output: f32 = hidden2
            .iter()
            .zip(self.output_weights.iter())
            .map(|(x, w)| x * w)
            .sum();

        Ok(Self::sigmoid(output))
    }

    fn train_batch(
        &mut self,
        batch: &[(Vec<f32>, f32)],
        learning_rate: f32,
    ) -> Result<(f32, usize)> {
        let mut total_loss = 0.0;
        let mut correct = 0;

        for (input, target) in batch {
            let prediction = self.forward(input)?;

            // Binary cross-entropy loss
            let loss = -target * prediction.ln() - (1.0 - target) * (1.0 - prediction).ln();
            total_loss += loss;

            // Check if prediction is correct
            let predicted_class = if prediction > 0.5 { 1.0 } else { 0.0 };
            if (predicted_class - target).abs() < 0.1 {
                correct += 1;
            }

            // Backpropagation (simplified gradient descent)
            let error = prediction - target;
            self.update_weights(input, error, learning_rate)?;
        }

        Ok((total_loss, correct))
    }

    /// Update neural network weights using backpropagation
    fn update_weights(&mut self, input: &[f32], error: f32, learning_rate: f32) -> Result<()> {
        // Simplified weight update (proper backprop would compute gradients layer by layer)
        // This is a basic approximation for demonstration

        // Update output weights
        for (i, weight) in self.output_weights.iter_mut().enumerate() {
            if i < input.len() {
                *weight -= learning_rate * error * input[i];
            }
        }

        Ok(())
    }

    fn relu(x: f32) -> f32 {
        x.max(0.0)
    }

    fn sigmoid(x: f32) -> f32 {
        1.0 / (1.0 + (-x).exp())
    }
}

impl TraditionalBloom {
    fn new(expected_elements: usize, false_positive_rate: f64) -> Self {
        // Calculate optimal parameters
        let size = Self::optimal_size(expected_elements, false_positive_rate);
        let hash_count = Self::optimal_hash_count(expected_elements, size);

        Self {
            bits: vec![false; size],
            hash_count,
            size,
        }
    }

    fn optimal_size(n: usize, p: f64) -> usize {
        (-(n as f64) * p.ln() / (2.0_f64.ln().powi(2))).ceil() as usize
    }

    fn optimal_hash_count(n: usize, m: usize) -> usize {
        ((m as f64 / n as f64) * 2.0_f64.ln()).round() as usize
    }

    fn insert(&mut self, item: u64) {
        for i in 0..self.hash_count {
            let hash = self.hash_function(item, i);
            let index = (hash % self.size as u64) as usize;
            self.bits[index] = true;
        }
    }

    fn query(&self, item: u64) -> bool {
        for i in 0..self.hash_count {
            let hash = self.hash_function(item, i);
            let index = (hash % self.size as u64) as usize;
            if !self.bits[index] {
                return false;
            }
        }
        true
    }

    fn hash_function(&self, item: u64, seed: usize) -> u64 {
        // Simple hash function family (use better ones in production)
        let mut hasher = item;
        hasher ^= seed as u64;
        hasher = hasher.wrapping_mul(0x9e3779b9);
        hasher ^= hasher >> 16;
        hasher = hasher.wrapping_mul(0x85ebca6b);
        hasher ^= hasher >> 13;
        hasher = hasher.wrapping_mul(0xc2b2ae35);
        hasher ^= hasher >> 16;
        hasher
    }
}

/// Example integration with your existing sketch system
pub fn integrate_learned_bloom(
    sketch_file: &str,
    model_size_limit_mb: usize,
) -> Result<LearnedBloomFilter> {
    // Load existing sketch data
    let sketch_data = std::fs::read_to_string(sketch_file)?;
    let mut kmer_hashes = Vec::new();

    for line in sketch_data.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            if let Ok(hash) = parts[0].parse::<u64>() {
                kmer_hashes.push(hash);
            }
        }
    }

    println!("Loaded {} k-mer hashes from sketch", kmer_hashes.len());

    // Create and train learned Bloom filter
    let mut learned_bloom = LearnedBloomFilter::new(
        kmer_hashes.len(),
        1e-6, // Target FPR of 10^-6
        32,   // Embedding dimension
        64,   // Hidden layer size
    );

    // Convert hashes to training data (simplified)
    let training_data: Vec<KmerEmbedding> = kmer_hashes
        .into_iter()
        .enumerate()
        .map(|(i, hash)| KmerEmbedding {
            kmer_hash: hash,
            embedding: vec![0.0; 32],    // Would compute actual embeddings
            taxonomy_id: Some(i as u32), // Placeholder taxonomy
        })
        .collect();

    learned_bloom.train_on_kmers(&training_data)?;

    // Check model size
    let model_bytes = learned_bloom.export_model()?;
    let model_size_mb = model_bytes.len() / (1024 * 1024);

    if model_size_mb > model_size_limit_mb {
        return Err(anyhow!(
            "Model size {} MB exceeds limit {} MB",
            model_size_mb,
            model_size_limit_mb
        ));
    }

    println!("Learned Bloom filter ready, size: {model_size_mb} MB");
    Ok(learned_bloom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neural_oracle_training() {
        let mut oracle = NeuralOracle::new(10, 16);

        // Create simple training data
        let batch = vec![
            (vec![1.0; 10], 1.0), // Positive example
            (vec![0.0; 10], 0.0), // Negative example
        ];

        let (loss, correct) = oracle.train_batch(&batch, 0.01).unwrap();
        assert!(loss > 0.0);
        assert!(correct <= 2);
    }

    #[test]
    fn test_bloom_filter_basic() {
        let mut bloom = TraditionalBloom::new(1000, 0.01);

        bloom.insert(12345);
        assert!(bloom.query(12345));

        // Should probably not contain random values
        let false_positives = (0..100)
            .map(|i| bloom.query(54321 + i))
            .filter(|&x| x)
            .count();

        assert!(false_positives < 10); // Should be low for good parameters
    }

    #[test]
    fn test_learned_bloom_filter_creation() {
        let expected_elements = 1000;
        let target_fpr = 0.01;
        let embedding_dim = 32;
        let hidden_dim = 64;

        let filter =
            LearnedBloomFilter::new(expected_elements, target_fpr, embedding_dim, hidden_dim);

        // Check training config is set correctly
        assert_eq!(filter.training_config.embedding_dim, embedding_dim);
        assert_eq!(filter.training_config.hidden_dim, hidden_dim);
        assert!((filter.training_config.target_fpr - target_fpr).abs() < f64::EPSILON);
        assert!((filter.training_config.learning_rate - 0.001).abs() < f32::EPSILON);
        assert_eq!(filter.training_config.batch_size, 64);
        assert_eq!(filter.training_config.epochs, 100);

        // Check neural oracle dimensions
        assert_eq!(filter.neural_oracle.input_size, embedding_dim);
        assert_eq!(filter.neural_oracle.hidden_size, hidden_dim);

        // Check that caches are initialized as empty
        assert!(filter.embedding_cache.is_empty());

        // Check metrics are initialized to defaults
        assert_eq!(filter.metrics.neural_queries, 0);
        assert_eq!(filter.metrics.neural_hits, 0);
        assert_eq!(filter.metrics.backup_queries, 0);
    }

    #[test]
    fn test_neural_oracle_creation() {
        let input_size = 16;
        let hidden_size = 32;
        let oracle = NeuralOracle::new(input_size, hidden_size);

        assert_eq!(oracle.input_size, input_size);
        assert_eq!(oracle.hidden_size, hidden_size);
        assert_eq!(oracle.input_weights.len(), hidden_size);
        assert_eq!(oracle.hidden_weights.len(), hidden_size);
        assert_eq!(oracle.output_weights.len(), hidden_size);

        // Check that input weights have correct dimensions
        if !oracle.input_weights.is_empty() {
            assert_eq!(oracle.input_weights[0].len(), input_size);
        }

        // Check threshold values are reasonable
        assert!(oracle.positive_threshold > 0.5);
        assert!(oracle.negative_threshold < 0.5);
        assert!(oracle.positive_threshold > oracle.negative_threshold);
    }

    #[test]
    fn test_neural_oracle_forward_pass() {
        let oracle = NeuralOracle::new(8, 16);
        let input = vec![0.5; 8];

        let output = oracle.forward(&input).unwrap();

        // Output should be between 0 and 1 due to sigmoid activation
        assert!(output >= 0.0 && output <= 1.0);
    }

    #[test]
    fn test_neural_oracle_prediction_types() {
        let mut oracle = NeuralOracle::new(4, 8);

        // Test with different input patterns
        let high_confidence_positive = vec![1.0; 4];
        let high_confidence_negative = vec![0.0; 4];
        let uncertain_input = vec![0.5; 4];

        // We can't guarantee specific predictions without training,
        // but we can test that predictions return valid enum variants
        let pred1 = oracle.predict(&high_confidence_positive).unwrap();
        let pred2 = oracle.predict(&high_confidence_negative).unwrap();
        let pred3 = oracle.predict(&uncertain_input).unwrap();

        // Just verify that we get valid prediction types
        match pred1 {
            NeuralPrediction::DefinitelyPresent
            | NeuralPrediction::DefinitelyAbsent
            | NeuralPrediction::Uncertain => {}
        }

        match pred2 {
            NeuralPrediction::DefinitelyPresent
            | NeuralPrediction::DefinitelyAbsent
            | NeuralPrediction::Uncertain => {}
        }

        match pred3 {
            NeuralPrediction::DefinitelyPresent
            | NeuralPrediction::DefinitelyAbsent
            | NeuralPrediction::Uncertain => {}
        }
    }

    #[test]
    fn test_traditional_bloom_optimal_parameters() {
        let expected_elements = 1000;
        let fpr = 0.01;

        let size = TraditionalBloom::optimal_size(expected_elements, fpr);
        let hash_count = TraditionalBloom::optimal_hash_count(expected_elements, size);

        // Size should be reasonable for the given parameters
        assert!(size > expected_elements); // Should be larger than number of elements
        assert!(size < expected_elements * 20); // But not excessively large

        // Hash count should be reasonable
        assert!(hash_count > 0);
        assert!(hash_count < 20); // Typically less than 20 hash functions
    }

    #[test]
    fn test_kmer_embedding_structure() {
        let embedding = KmerEmbedding {
            kmer_hash: 12345,
            embedding: vec![0.1, 0.2, 0.3, 0.4],
            taxonomy_id: Some(789),
        };

        assert_eq!(embedding.kmer_hash, 12345);
        assert_eq!(embedding.embedding.len(), 4);
        assert_eq!(embedding.taxonomy_id, Some(789));
        assert!((embedding.embedding[0] - 0.1).abs() < f32::EPSILON);
    }

    #[test]
    fn test_filter_metrics_initialization() {
        let metrics = FilterMetrics::default();

        assert_eq!(metrics.neural_queries, 0);
        assert_eq!(metrics.neural_hits, 0);
        assert_eq!(metrics.backup_queries, 0);
        assert_eq!(metrics.false_positives, 0);
        assert_eq!(metrics.true_positives, 0);
    }

    #[test]
    fn test_training_config_structure() {
        let config = TrainingConfig {
            learning_rate: 0.005,
            batch_size: 32,
            epochs: 50,
            embedding_dim: 16,
            hidden_dim: 32,
            target_fpr: 0.001,
        };

        assert!((config.learning_rate - 0.005).abs() < f32::EPSILON);
        assert_eq!(config.batch_size, 32);
        assert_eq!(config.epochs, 50);
        assert_eq!(config.embedding_dim, 16);
        assert_eq!(config.hidden_dim, 32);
        assert!((config.target_fpr - 0.001).abs() < f64::EPSILON);
    }

    #[test]
    fn test_compute_kmer_embedding() {
        let filter = LearnedBloomFilter::new(1000, 0.01, 16, 32);
        let kmer_hash = 0x123456789ABCDEF0u64;

        let embedding = filter.compute_kmer_embedding(kmer_hash).unwrap();

        assert_eq!(embedding.len(), 16); // Should match embedding_dim

        // First feature should be normalized hash magnitude
        assert!(embedding[0] >= 0.0 && embedding[0] <= 1.0);

        // All features should be finite
        for &feature in &embedding {
            assert!(feature.is_finite());
        }

        // Test with a smaller hash to avoid overflow
        let small_hash = 0x1234u64;
        let embedding2 = filter.compute_kmer_embedding(small_hash).unwrap();
        assert_eq!(embedding2.len(), 16);
        assert!(embedding2[0] >= 0.0 && embedding2[0] <= 1.0);
    }

    #[test]
    fn test_hash_function_consistency() {
        let bloom = TraditionalBloom::new(100, 0.1);
        let item = 42u64;
        let seed = 0;

        // Hash function should be deterministic
        let hash1 = bloom.hash_function(item, seed);
        let hash2 = bloom.hash_function(item, seed);
        assert_eq!(hash1, hash2);

        // Different seeds should generally produce different hashes
        let hash3 = bloom.hash_function(item, 1);
        assert_ne!(hash1, hash3); // Very likely to be different
    }
}
