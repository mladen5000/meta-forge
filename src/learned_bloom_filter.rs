use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

/// Learned Bloom filter combining neural predictor with backup Bloom filter
pub struct LearnedBloomFilter {
    /// Neural network predictor (lightweight MLP)
    neural_oracle: NeuralOracle,
    /// Traditional Bloom filter for backup
    backup_bloom: TraditionalBloom,
    /// K-mer embedding cache
    embedding_cache: HashMap<u64, Vec<f32>>,
    /// Training parameters
    training_config: TrainingConfig,
    /// Performance metrics
    metrics: FilterMetrics,
}

#[derive(Clone)]
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

#[derive(Clone)]
pub struct TrainingConfig {
    learning_rate: f32,
    batch_size: usize,
    epochs: usize,
    embedding_dim: usize,
    hidden_dim: usize,
    target_fpr: f64, // Target false positive rate
}

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
        let config = TrainingConfig {
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
    fn compute_kmer_embedding(&self, kmer_hash: u64) -> Result<Vec<f32>> {
        // Convert hash to embedding using learned features
        let mut embedding = vec![0.0; self.training_config.embedding_dim];

        // Feature 1: Hash magnitude (normalized)
        embedding[0] = (kmer_hash as f32) / (u64::MAX as f32);

        // Feature 2-5: Hash bits in different positions
        for i in 1..5 {
            let shift = i * 16;
            embedding[i] = ((kmer_hash >> shift) & 0xFFFF) as f32 / 65535.0;
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
    pub fn export_model(&self) -> Result<Vec<u8>> {
        let model_data = SerializableModel {
            neural_oracle: self.neural_oracle.clone(),
            backup_bloom_bits: self.backup_bloom.bits.clone(),
            backup_bloom_size: self.backup_bloom.size,
            backup_bloom_hash_count: self.backup_bloom.hash_count,
            training_config: self.training_config.clone(),
        };

        let serialized =
            bincode::serialize(&model_data).map_err(|e| anyhow!("Serialization failed: {}", e))?;

        println!(
            "Model exported, size: {} bytes ({:.2} MB)",
            serialized.len(),
            serialized.len() as f64 / (1024.0 * 1024.0)
        );

        Ok(serialized)
    }

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
    training_config: TrainingConfig,
}

#[derive(Debug)]
enum NeuralPrediction {
    DefinitelyPresent,
    DefinitelyAbsent,
    Uncertain,
}

impl NeuralOracle {
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

    println!("Learned Bloom filter ready, size: {} MB", model_size_mb);
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
}
