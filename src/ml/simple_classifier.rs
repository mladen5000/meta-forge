//! Simple ML-based contig classification inspired by SemiBin2
//!
//! This module implements a lightweight, self-supervised approach for contig binning
//! using k-mer frequency features and coverage information.
//!
//! Key features:
//! - K-mer frequency-based embeddings (tetra-nucleotide frequencies)
//! - Coverage-based features (depth normalization)
//! - Simple distance-based clustering
//! - No external ML dependencies required

use ahash::AHashMap;
use anyhow::{Context, Result};
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};

use crate::core::data_structures::Contig;

/// Configuration for the simple ML classifier
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimpleClassifierConfig {
    /// K-mer size for feature extraction (typically 4 for tetranucleotides)
    pub kmer_size: usize,
    /// Minimum contig length to classify
    pub min_contig_length: usize,
    /// Number of bins/clusters to create
    pub num_bins: usize,
    /// Include coverage features
    pub use_coverage_features: bool,
    /// Normalization method
    pub normalization: NormalizationMethod,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum NormalizationMethod {
    /// Z-score normalization (mean=0, std=1)
    ZScore,
    /// Min-max normalization (0-1 range)
    MinMax,
    /// No normalization
    None,
}

impl Default for SimpleClassifierConfig {
    fn default() -> Self {
        Self {
            kmer_size: 4, // Tetranucleotide frequencies
            min_contig_length: 1000,
            num_bins: 10,
            use_coverage_features: true,
            normalization: NormalizationMethod::ZScore,
        }
    }
}

/// Simple contig classifier using k-mer frequencies and coverage
pub struct SimpleContigClassifier {
    config: SimpleClassifierConfig,
    /// Pre-computed k-mer vocabulary for indexing
    kmer_vocab: Vec<String>,
    /// Dimension of feature vectors
    feature_dim: usize,
}

/// Classification result for a contig
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContigClassification {
    pub contig_id: usize,
    pub bin_id: usize,
    pub confidence: f64,
    pub features: Vec<f64>,
}

impl SimpleContigClassifier {
    /// Create a new classifier with given configuration
    pub fn new(config: SimpleClassifierConfig) -> Result<Self> {
        // Generate all possible k-mers for the vocabulary
        let kmer_vocab = Self::generate_kmer_vocabulary(config.kmer_size);
        let mut feature_dim = kmer_vocab.len();

        // Add dimensions for coverage features if enabled
        if config.use_coverage_features {
            feature_dim += 1; // Add 1 for normalized coverage
        }

        Ok(Self {
            config,
            kmer_vocab,
            feature_dim,
        })
    }

    /// Generate all possible k-mers of a given size
    fn generate_kmer_vocabulary(k: usize) -> Vec<String> {
        let bases = ['A', 'C', 'G', 'T'];
        let mut kmers = Vec::new();

        fn generate_recursive(
            current: String,
            k: usize,
            bases: &[char],
            kmers: &mut Vec<String>,
        ) {
            if current.len() == k {
                kmers.push(current);
                return;
            }

            for &base in bases {
                let mut next = current.clone();
                next.push(base);
                generate_recursive(next, k, bases, kmers);
            }
        }

        generate_recursive(String::new(), k, &bases, &mut kmers);
        kmers.sort();
        kmers
    }

    /// Extract feature vector from a contig
    pub fn extract_features(&self, contig: &Contig) -> Result<Array1<f64>> {
        let mut features = vec![0.0; self.feature_dim];

        // Extract k-mer frequencies
        let kmer_freqs = self.compute_kmer_frequencies(&contig.sequence)?;

        // Map k-mer frequencies to feature vector
        for (kmer, freq) in kmer_freqs.iter() {
            if let Some(idx) = self.kmer_vocab.iter().position(|k| k == kmer) {
                features[idx] = *freq;
            }
        }

        // Add coverage feature if enabled
        if self.config.use_coverage_features {
            let coverage_idx = self.kmer_vocab.len();
            features[coverage_idx] = contig.coverage as f64;
        }

        Ok(Array1::from_vec(features))
    }

    /// Compute k-mer frequency distribution for a sequence
    fn compute_kmer_frequencies(&self, sequence: &str) -> Result<AHashMap<String, f64>> {
        let k = self.config.kmer_size;
        let mut kmer_counts = AHashMap::new();
        let mut total_kmers = 0;

        // Count k-mers
        for window in sequence.as_bytes().windows(k) {
            if let Ok(kmer) = std::str::from_utf8(window) {
                // Convert to uppercase and skip if contains ambiguous bases
                let kmer_upper = kmer.to_uppercase();
                if kmer_upper.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T')) {
                    *kmer_counts.entry(kmer_upper).or_insert(0) += 1;
                    total_kmers += 1;
                }
            }
        }

        // Convert counts to frequencies
        let mut kmer_freqs = AHashMap::new();
        if total_kmers > 0 {
            for (kmer, count) in kmer_counts.iter() {
                kmer_freqs.insert(kmer.clone(), *count as f64 / total_kmers as f64);
            }
        }

        Ok(kmer_freqs)
    }

    /// Classify a batch of contigs into bins
    pub fn classify_contigs(&self, contigs: &[Contig]) -> Result<Vec<ContigClassification>> {
        // Filter contigs by minimum length
        let valid_contigs: Vec<&Contig> = contigs
            .iter()
            .filter(|c| c.sequence.len() >= self.config.min_contig_length)
            .collect();

        if valid_contigs.is_empty() {
            return Ok(Vec::new());
        }

        // Extract features for all contigs
        let mut feature_matrix = Vec::new();
        let mut contig_ids = Vec::new();

        for contig in valid_contigs.iter() {
            match self.extract_features(contig) {
                Ok(features) => {
                    feature_matrix.push(features);
                    contig_ids.push(contig.id);
                }
                Err(e) => {
                    tracing::warn!("Failed to extract features for contig {}: {}", contig.id, e);
                    continue;
                }
            }
        }

        if feature_matrix.is_empty() {
            return Ok(Vec::new());
        }

        // Normalize features
        let normalized_features = self.normalize_features(&feature_matrix)?;

        // Perform simple k-means-like clustering
        let bin_assignments = self.cluster_contigs(&normalized_features)?;

        // Create classification results
        let mut classifications = Vec::new();
        for (idx, &contig_id) in contig_ids.iter().enumerate() {
            let bin_id = bin_assignments[idx];
            let confidence = self.compute_confidence(
                &normalized_features[idx],
                &normalized_features,
                &bin_assignments,
                bin_id,
            )?;

            classifications.push(ContigClassification {
                contig_id,
                bin_id,
                confidence,
                features: normalized_features[idx].to_vec(),
            });
        }

        Ok(classifications)
    }

    /// Normalize feature vectors
    fn normalize_features(&self, features: &[Array1<f64>]) -> Result<Vec<Array1<f64>>> {
        let n_samples = features.len();
        let n_features = self.feature_dim;

        if n_samples == 0 {
            return Ok(Vec::new());
        }

        match self.config.normalization {
            NormalizationMethod::ZScore => {
                // Compute mean and std for each feature
                let mut means = vec![0.0; n_features];
                let mut stds = vec![0.0; n_features];

                // Compute means
                for feature_vec in features.iter() {
                    for (i, &val) in feature_vec.iter().enumerate() {
                        means[i] += val;
                    }
                }
                for mean in means.iter_mut() {
                    *mean /= n_samples as f64;
                }

                // Compute standard deviations
                for feature_vec in features.iter() {
                    for (i, &val) in feature_vec.iter().enumerate() {
                        stds[i] += (val - means[i]).powi(2);
                    }
                }
                for std in stds.iter_mut() {
                    *std = (*std / n_samples as f64).sqrt();
                    if *std < 1e-10 {
                        *std = 1.0; // Avoid division by zero
                    }
                }

                // Normalize
                let normalized: Vec<Array1<f64>> = features
                    .iter()
                    .map(|feature_vec| {
                        let normalized_vals: Vec<f64> = feature_vec
                            .iter()
                            .enumerate()
                            .map(|(i, &val)| (val - means[i]) / stds[i])
                            .collect();
                        Array1::from_vec(normalized_vals)
                    })
                    .collect();

                Ok(normalized)
            }
            NormalizationMethod::MinMax => {
                // Compute min and max for each feature
                let mut mins = vec![f64::INFINITY; n_features];
                let mut maxs = vec![f64::NEG_INFINITY; n_features];

                for feature_vec in features.iter() {
                    for (i, &val) in feature_vec.iter().enumerate() {
                        mins[i] = mins[i].min(val);
                        maxs[i] = maxs[i].max(val);
                    }
                }

                // Normalize
                let normalized: Vec<Array1<f64>> = features
                    .iter()
                    .map(|feature_vec| {
                        let normalized_vals: Vec<f64> = feature_vec
                            .iter()
                            .enumerate()
                            .map(|(i, &val)| {
                                let range = maxs[i] - mins[i];
                                if range < 1e-10 {
                                    0.5 // If no variation, set to middle
                                } else {
                                    (val - mins[i]) / range
                                }
                            })
                            .collect();
                        Array1::from_vec(normalized_vals)
                    })
                    .collect();

                Ok(normalized)
            }
            NormalizationMethod::None => Ok(features.to_vec()),
        }
    }

    /// Simple k-means clustering algorithm
    fn cluster_contigs(&self, features: &[Array1<f64>]) -> Result<Vec<usize>> {
        let n_samples = features.len();
        let k = self.config.num_bins.min(n_samples);

        if n_samples == 0 {
            return Ok(Vec::new());
        }

        // Initialize centroids using k-means++ strategy
        let mut centroids = Vec::new();

        // First centroid: random sample
        let first_idx = fastrand::usize(..n_samples);
        centroids.push(features[first_idx].clone());

        // Remaining centroids: weighted by distance to nearest centroid
        for _ in 1..k {
            let mut distances = vec![f64::INFINITY; n_samples];

            for (i, feature) in features.iter().enumerate() {
                for centroid in centroids.iter() {
                    let dist = Self::euclidean_distance(feature, centroid);
                    distances[i] = distances[i].min(dist);
                }
            }

            // Select next centroid proportional to squared distance
            let total_dist: f64 = distances.iter().map(|d| d * d).sum();
            let mut threshold = fastrand::f64() * total_dist;

            let mut next_idx = 0;
            for (i, &dist) in distances.iter().enumerate() {
                threshold -= dist * dist;
                if threshold <= 0.0 {
                    next_idx = i;
                    break;
                }
            }

            centroids.push(features[next_idx].clone());
        }

        // Run k-means iterations
        let max_iterations = 100;
        let mut assignments = vec![0; n_samples];

        for _ in 0..max_iterations {
            let mut changed = false;

            // Assign each point to nearest centroid
            for (i, feature) in features.iter().enumerate() {
                let mut min_dist = f64::INFINITY;
                let mut best_cluster = 0;

                for (cluster_id, centroid) in centroids.iter().enumerate() {
                    let dist = Self::euclidean_distance(feature, centroid);
                    if dist < min_dist {
                        min_dist = dist;
                        best_cluster = cluster_id;
                    }
                }

                if assignments[i] != best_cluster {
                    assignments[i] = best_cluster;
                    changed = true;
                }
            }

            if !changed {
                break;
            }

            // Update centroids
            for cluster_id in 0..k {
                let cluster_points: Vec<&Array1<f64>> = features
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| assignments[*i] == cluster_id)
                    .map(|(_, f)| f)
                    .collect();

                if !cluster_points.is_empty() {
                    centroids[cluster_id] = Self::compute_centroid(&cluster_points);
                }
            }
        }

        Ok(assignments)
    }

    /// Compute euclidean distance between two feature vectors
    fn euclidean_distance(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt()
    }

    /// Compute centroid of a set of points
    fn compute_centroid(points: &[&Array1<f64>]) -> Array1<f64> {
        let n = points.len();
        let dim = points[0].len();

        let mut centroid = vec![0.0; dim];
        for point in points {
            for (i, &val) in point.iter().enumerate() {
                centroid[i] += val;
            }
        }

        for val in centroid.iter_mut() {
            *val /= n as f64;
        }

        Array1::from_vec(centroid)
    }

    /// Compute classification confidence based on distance to cluster center
    fn compute_confidence(
        &self,
        feature: &Array1<f64>,
        all_features: &[Array1<f64>],
        assignments: &[usize],
        bin_id: usize,
    ) -> Result<f64> {
        // Find all points in the same cluster
        let cluster_points: Vec<&Array1<f64>> = all_features
            .iter()
            .enumerate()
            .filter(|(i, _)| assignments[*i] == bin_id)
            .map(|(_, f)| f)
            .collect();

        if cluster_points.is_empty() {
            return Ok(0.0);
        }

        // Compute cluster centroid
        let centroid = Self::compute_centroid(&cluster_points);

        // Compute distance to centroid
        let dist_to_centroid = Self::euclidean_distance(feature, &centroid);

        // Compute average distance of cluster members to centroid
        let avg_cluster_dist: f64 = cluster_points
            .iter()
            .map(|p| Self::euclidean_distance(p, &centroid))
            .sum::<f64>()
            / cluster_points.len() as f64;

        // Confidence is inversely proportional to relative distance
        // Higher confidence if closer to centroid than average
        let confidence = if avg_cluster_dist > 0.0 {
            (1.0 - (dist_to_centroid / (avg_cluster_dist * 2.0))).max(0.0).min(1.0)
        } else {
            1.0
        };

        Ok(confidence)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_vocabulary_generation() {
        let vocab = SimpleContigClassifier::generate_kmer_vocabulary(2);
        assert_eq!(vocab.len(), 16); // 4^2 = 16 dimers
        assert!(vocab.contains(&"AA".to_string()));
        assert!(vocab.contains(&"TT".to_string()));

        let vocab4 = SimpleContigClassifier::generate_kmer_vocabulary(4);
        assert_eq!(vocab4.len(), 256); // 4^4 = 256 tetranucleotides
    }

    #[test]
    fn test_feature_extraction() {
        let config = SimpleClassifierConfig {
            kmer_size: 4,
            use_coverage_features: true,
            ..Default::default()
        };

        let classifier = SimpleContigClassifier::new(config).unwrap();

        let contig = Contig {
            id: 1,
            sequence: "ATCGATCGATCGATCG".to_string(),
            coverage: 10.0,
            length: 16,
            contig_type: crate::core::data_structures::ContigType::Linear,
            node_path: vec![],
        };

        let features = classifier.extract_features(&contig).unwrap();

        // Should have 256 k-mer features + 1 coverage feature
        assert_eq!(features.len(), 257);

        // Coverage feature should be normalized coverage value
        assert_eq!(features[256], 10.0);
    }

    #[test]
    fn test_kmer_frequency_computation() {
        let config = SimpleClassifierConfig::default();
        let classifier = SimpleContigClassifier::new(config).unwrap();

        // Sequence with known k-mer composition
        let sequence = "ATCGATCGATCG"; // Contains "ATCG" repeated
        let freqs = classifier.compute_kmer_frequencies(sequence).unwrap();

        // Should have high frequency for "ATCG"
        assert!(freqs.contains_key("ATCG"));
        assert!(freqs["ATCG"] > 0.0);
    }

    #[test]
    fn test_simple_clustering() {
        let config = SimpleClassifierConfig {
            kmer_size: 4,
            num_bins: 2,
            min_contig_length: 10,
            use_coverage_features: true,
            normalization: NormalizationMethod::ZScore,
        };

        let classifier = SimpleContigClassifier::new(config).unwrap();

        // Create test contigs with different compositions
        let contigs = vec![
            Contig {
                id: 1,
                sequence: "ATATATATAT".to_string(), // AT-rich
                coverage: 10.0,
                length: 10,
                contig_type: crate::core::data_structures::ContigType::Linear,
                node_path: vec![],
            },
            Contig {
                id: 2,
                sequence: "GCGCGCGCGC".to_string(), // GC-rich
                coverage: 15.0,
                length: 10,
                contig_type: crate::core::data_structures::ContigType::Linear,
                node_path: vec![],
            },
            Contig {
                id: 3,
                sequence: "ATATATATATAT".to_string(), // AT-rich (similar to 1)
                coverage: 12.0,
                length: 12,
                contig_type: crate::core::data_structures::ContigType::Linear,
                node_path: vec![],
            },
        ];

        let classifications = classifier.classify_contigs(&contigs).unwrap();

        assert_eq!(classifications.len(), 3);

        // Contigs 1 and 3 should be in the same bin (similar composition)
        assert_eq!(
            classifications[0].bin_id,
            classifications[2].bin_id,
            "AT-rich contigs should cluster together"
        );

        // All classifications should have confidence scores
        for classification in classifications {
            assert!(
                classification.confidence >= 0.0 && classification.confidence <= 1.0,
                "Confidence should be between 0 and 1"
            );
        }
    }

    #[test]
    fn test_normalization() {
        let config = SimpleClassifierConfig {
            normalization: NormalizationMethod::ZScore,
            ..Default::default()
        };

        let classifier = SimpleContigClassifier::new(config).unwrap();

        let features = vec![
            Array1::from_vec(vec![1.0, 2.0, 3.0]),
            Array1::from_vec(vec![4.0, 5.0, 6.0]),
            Array1::from_vec(vec![7.0, 8.0, 9.0]),
        ];

        let normalized = classifier.normalize_features(&features).unwrap();

        // Check that normalization produces valid results
        assert_eq!(normalized.len(), 3);

        // Z-score normalization should have mean ≈ 0 and std ≈ 1 for each feature
        for feature_idx in 0..3 {
            let values: Vec<f64> = normalized.iter().map(|v| v[feature_idx]).collect();
            let mean: f64 = values.iter().sum::<f64>() / values.len() as f64;
            assert!(mean.abs() < 1e-10, "Mean should be close to 0 after Z-score normalization");
        }
    }
}
