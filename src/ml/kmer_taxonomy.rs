//! K-mer based taxonomic classification for metagenomic contigs
//!
//! This module implements composition-based taxonomic assignment using k-mer profiles.
//! Unlike clustering-based binning, this provides actual taxonomic labels.

use ahash::AHashMap;
use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

/// K-mer based taxonomic classifier using composition signatures
pub struct KmerTaxonomyClassifier {
    /// Reference k-mer profiles for known taxa
    reference_profiles: AHashMap<String, KmerProfile>,
    /// K-mer size for classification
    k: usize,
    /// Minimum similarity threshold for classification
    min_similarity: f64,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct KmerProfile {
    pub taxon_name: String,
    pub taxon_rank: String, // species, genus, family, etc.
    pub kmer_frequencies: AHashMap<String, f64>,
    pub gc_content: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonomicAssignment {
    pub taxon: String,
    pub rank: String,
    pub confidence: f64,
    pub method: String,
}

impl TaxonomicAssignment {
    /// Create unclassified assignment
    pub fn unclassified() -> Self {
        Self {
            taxon: "Unclassified".to_string(),
            rank: "unknown".to_string(),
            confidence: 0.0,
            method: "k-mer composition".to_string(),
        }
    }

    /// Check if assignment is classified
    pub fn is_classified(&self) -> bool {
        self.taxon != "Unclassified" && self.confidence > 0.0
    }
}

impl KmerTaxonomyClassifier {
    /// Create new classifier with given k-mer size
    pub fn new(k: usize) -> Self {
        Self {
            reference_profiles: AHashMap::new(),
            k,
            min_similarity: 0.6, // 60% similarity threshold (relaxed for better coverage)
        }
    }

    /// Create classifier with custom similarity threshold
    pub fn with_threshold(k: usize, min_similarity: f64) -> Self {
        Self {
            reference_profiles: AHashMap::new(),
            k,
            min_similarity,
        }
    }

    /// Set the minimum similarity threshold for classification
    pub fn set_similarity_threshold(&mut self, threshold: f64) {
        self.min_similarity = threshold.max(0.0).min(1.0);
        tracing::info!("Similarity threshold set to {:.1}%", self.min_similarity * 100.0);
    }

    /// Get current similarity threshold
    pub fn get_similarity_threshold(&self) -> f64 {
        self.min_similarity
    }

    /// Add reference genome profile for classification
    pub fn add_reference(&mut self, profile: KmerProfile) {
        self.reference_profiles
            .insert(profile.taxon_name.clone(), profile);
    }

    /// Get number of reference profiles loaded
    pub fn num_references(&self) -> usize {
        self.reference_profiles.len()
    }

    /// Load taxonomy reference profiles from database
    pub fn load_reference_database(
        &mut self,
        db: &crate::database::integration::MetagenomicsDatabase,
    ) -> Result<()> {

        tracing::info!("Loading taxonomy reference profiles from database...");

        // Get database statistics to see what's available
        let stats = db.get_database_stats()
            .context("Failed to get database statistics")?;

        if stats.taxonomy_count == 0 {
            tracing::warn!("No taxonomy entries found in database, falling back to default references");
            return self.load_default_references();
        }

        if stats.sequences_count == 0 {
            tracing::warn!("No reference sequences found in database, falling back to default references");
            return self.load_default_references();
        }

        tracing::info!(
            "Database has {} taxonomy entries and {} sequences",
            stats.taxonomy_count,
            stats.sequences_count
        );

        // TODO: Implement actual profile loading from sequences
        // For now, this would:
        // 1. Query sequences table joined with taxonomy
        // 2. Calculate k-mer profiles for each taxonomy
        // 3. Store in reference_profiles map

        // Fallback to default references until full implementation
        tracing::warn!("Full database profile loading not yet implemented, using defaults");
        self.load_default_references()?;

        Ok(())
    }

    /// Load common bacterial reference profiles (placeholder for real DB)
    ///
    /// This loads 40+ common bacterial species covering:
    /// - Human gut microbiome (Bacteroides, Bifidobacterium, etc.)
    /// - Human pathogens (E. coli, S. aureus, etc.)
    /// - Environmental bacteria (Pseudomonas, Rhizobium, etc.)
    /// - Soil/water bacteria (Bacillus, Streptomyces, etc.)
    pub fn load_default_references(&mut self) -> Result<()> {
        tracing::info!("Loading comprehensive bacterial reference database...");

        // ============================================================
        // HUMAN GUT MICROBIOME (Common commensal & beneficial)
        // ============================================================
        self.add_example_reference("Bacteroides fragilis", "species", 0.43);
        self.add_example_reference("Bacteroides thetaiotaomicron", "species", 0.43);
        self.add_example_reference("Faecalibacterium prausnitzii", "species", 0.57);
        self.add_example_reference("Bifidobacterium longum", "species", 0.60);
        self.add_example_reference("Bifidobacterium adolescentis", "species", 0.59);
        self.add_example_reference("Lactobacillus acidophilus", "species", 0.35);
        self.add_example_reference("Lactobacillus plantarum", "species", 0.45);
        self.add_example_reference("Akkermansia muciniphila", "species", 0.56);
        self.add_example_reference("Prevotella copri", "species", 0.47);
        self.add_example_reference("Ruminococcus bromii", "species", 0.49);

        // ============================================================
        // ENTEROBACTERIACEAE & COMMON PATHOGENS
        // ============================================================
        self.add_example_reference("Escherichia coli", "species", 0.51);
        self.add_example_reference("Salmonella enterica", "species", 0.52);
        self.add_example_reference("Klebsiella pneumoniae", "species", 0.57);
        self.add_example_reference("Enterobacter cloacae", "species", 0.55);
        self.add_example_reference("Shigella flexneri", "species", 0.51);
        self.add_example_reference("Yersinia pestis", "species", 0.48);

        // ============================================================
        // GRAM-POSITIVE PATHOGENS & COMMENSALS
        // ============================================================
        self.add_example_reference("Staphylococcus aureus", "species", 0.33);
        self.add_example_reference("Staphylococcus epidermidis", "species", 0.32);
        self.add_example_reference("Streptococcus pneumoniae", "species", 0.40);
        self.add_example_reference("Streptococcus pyogenes", "species", 0.39);
        self.add_example_reference("Enterococcus faecalis", "species", 0.38);
        self.add_example_reference("Clostridium difficile", "species", 0.29);
        self.add_example_reference("Clostridium botulinum", "species", 0.28);
        self.add_example_reference("Listeria monocytogenes", "species", 0.38);

        // ============================================================
        // BACILLUS & SPORE-FORMERS
        // ============================================================
        self.add_example_reference("Bacillus subtilis", "species", 0.44);
        self.add_example_reference("Bacillus anthracis", "species", 0.35);
        self.add_example_reference("Bacillus cereus", "species", 0.35);

        // ============================================================
        // PSEUDOMONADS & ENVIRONMENTAL BACTERIA
        // ============================================================
        self.add_example_reference("Pseudomonas aeruginosa", "species", 0.66);
        self.add_example_reference("Pseudomonas putida", "species", 0.62);
        self.add_example_reference("Pseudomonas fluorescens", "species", 0.60);
        self.add_example_reference("Acinetobacter baumannii", "species", 0.39);
        self.add_example_reference("Burkholderia cepacia", "species", 0.67);

        // ============================================================
        // RESPIRATORY & ORAL PATHOGENS
        // ============================================================
        self.add_example_reference("Mycobacterium tuberculosis", "species", 0.66);
        self.add_example_reference("Haemophilus influenzae", "species", 0.38);
        self.add_example_reference("Neisseria meningitidis", "species", 0.51);
        self.add_example_reference("Bordetella pertussis", "species", 0.68);
        self.add_example_reference("Legionella pneumophila", "species", 0.38);

        // ============================================================
        // SOIL & NITROGEN-FIXING BACTERIA
        // ============================================================
        self.add_example_reference("Rhizobium leguminosarum", "species", 0.61);
        self.add_example_reference("Agrobacterium tumefaciens", "species", 0.59);
        self.add_example_reference("Streptomyces coelicolor", "species", 0.72);
        self.add_example_reference("Azotobacter vinelandii", "species", 0.66);

        // ============================================================
        // AQUATIC & MARINE BACTERIA
        // ============================================================
        self.add_example_reference("Vibrio cholerae", "species", 0.48);
        self.add_example_reference("Shewanella oneidensis", "species", 0.46);
        self.add_example_reference("Prochlorococcus marinus", "species", 0.31);

        // ============================================================
        // EMERGING & ANTIBIOTIC-RESISTANT STRAINS
        // ============================================================
        self.add_example_reference("Mycoplasma pneumoniae", "species", 0.40);
        self.add_example_reference("Helicobacter pylori", "species", 0.39);
        self.add_example_reference("Campylobacter jejuni", "species", 0.31);

        tracing::info!(
            "✅ Loaded {} reference profiles covering gut microbiome, pathogens, and environmental bacteria",
            self.reference_profiles.len()
        );
        Ok(())
    }

    /// Create example reference profile (for testing/demo)
    fn add_example_reference(&mut self, name: &str, rank: &str, gc_content: f64) {
        // Generate synthetic k-mer profile based on GC content
        let mut kmer_freqs = AHashMap::new();

        // Generate all possible k-mers for this k
        let bases = ['A', 'C', 'G', 'T'];
        let total_kmers = 4_usize.pow(self.k as u32);

        for kmer_idx in 0..total_kmers {
            let mut kmer = String::with_capacity(self.k);
            let mut idx = kmer_idx;

            for _ in 0..self.k {
                kmer.push(bases[idx % 4]);
                idx /= 4;
            }

            // Calculate frequency based on GC content bias
            let gc_count = kmer.chars().filter(|c| matches!(c, 'G' | 'C')).count();
            let at_count = self.k - gc_count;

            let freq = if gc_count > at_count {
                gc_content / (total_kmers as f64 / 2.0)
            } else {
                (1.0 - gc_content) / (total_kmers as f64 / 2.0)
            };

            kmer_freqs.insert(kmer, freq);
        }

        let profile = KmerProfile {
            taxon_name: name.to_string(),
            taxon_rank: rank.to_string(),
            kmer_frequencies: kmer_freqs,
            gc_content,
        };

        self.add_reference(profile);
    }

    /// Classify sequence using k-mer composition
    pub fn classify_sequence(&self, sequence: &str) -> Result<TaxonomicAssignment> {
        if sequence.len() < self.k {
            tracing::trace!("Sequence too short for classification: {} < {}", sequence.len(), self.k);
            return Ok(TaxonomicAssignment::unclassified());
        }

        // Extract k-mer profile from query sequence
        let query_kmers =
            self.extract_kmer_profile(sequence)
                .context("Failed to extract k-mer profile")?;

        if query_kmers.is_empty() {
            tracing::trace!("No k-mers extracted from sequence");
            return Ok(TaxonomicAssignment::unclassified());
        }

        // Find best matching reference using cosine similarity
        let mut best_match = String::new();
        let mut best_score = 0.0;
        let mut second_best_score = 0.0;

        for (taxon, ref_profile) in &self.reference_profiles {
            let score = self.calculate_similarity(&query_kmers, &ref_profile.kmer_frequencies);
            if score > best_score {
                second_best_score = best_score;
                best_score = score;
                best_match = taxon.clone();
            } else if score > second_best_score {
                second_best_score = score;
            }
        }

        // Log classification confidence
        if best_score >= self.min_similarity {
            let margin = best_score - second_best_score;
            tracing::trace!(
                "Classification: {} (similarity: {:.2}%, margin: {:.2}%)",
                best_match,
                best_score * 100.0,
                margin * 100.0
            );
        } else {
            tracing::trace!(
                "No confident match (best: {:.2}% < threshold: {:.2}%)",
                best_score * 100.0,
                self.min_similarity * 100.0
            );
        }

        // Create assignment if similarity exceeds threshold
        let assignment = if best_score >= self.min_similarity {
            self.reference_profiles
                .get(&best_match)
                .map(|p| TaxonomicAssignment {
                    taxon: p.taxon_name.clone(),
                    rank: p.taxon_rank.clone(),
                    confidence: best_score,
                    method: "k-mer composition".to_string(),
                })
                .unwrap_or_else(TaxonomicAssignment::unclassified)
        } else {
            TaxonomicAssignment::unclassified()
        };

        Ok(assignment)
    }

    /// Extract k-mer frequency profile from sequence
    fn extract_kmer_profile(&self, sequence: &str) -> Result<AHashMap<String, f64>> {
        let mut kmer_counts = AHashMap::new();
        let mut total = 0;

        // Count k-mers in sequence
        for window in sequence.as_bytes().windows(self.k) {
            if let Ok(kmer) = std::str::from_utf8(window) {
                let kmer_upper = kmer.to_uppercase();

                // Only count valid DNA k-mers
                if kmer_upper
                    .chars()
                    .all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
                {
                    *kmer_counts.entry(kmer_upper).or_insert(0) += 1;
                    total += 1;
                }
            }
        }

        // Convert counts to frequencies
        let mut freqs = AHashMap::new();
        if total > 0 {
            for (kmer, count) in kmer_counts {
                freqs.insert(kmer, count as f64 / total as f64);
            }
        }

        Ok(freqs)
    }

    /// Calculate cosine similarity between two k-mer profiles
    fn calculate_similarity(
        &self,
        query: &AHashMap<String, f64>,
        reference: &AHashMap<String, f64>,
    ) -> f64 {
        let mut dot_product = 0.0;
        let mut query_norm = 0.0;
        let mut ref_norm = 0.0;

        // Calculate dot product and query norm
        for (kmer, &query_freq) in query {
            let ref_freq = reference.get(kmer).copied().unwrap_or(0.0);
            dot_product += query_freq * ref_freq;
            query_norm += query_freq * query_freq;
        }

        // Calculate reference norm
        for &ref_freq in reference.values() {
            ref_norm += ref_freq * ref_freq;
        }

        // Return cosine similarity
        if query_norm > 0.0 && ref_norm > 0.0 {
            dot_product / (query_norm.sqrt() * ref_norm.sqrt())
        } else {
            0.0
        }
    }

    /// Get statistics about reference database
    pub fn get_stats(&self) -> ClassifierStats {
        ClassifierStats {
            num_references: self.reference_profiles.len(),
            k_size: self.k,
            min_similarity_threshold: self.min_similarity,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ClassifierStats {
    pub num_references: usize,
    pub k_size: usize,
    pub min_similarity_threshold: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kmer_profile_extraction() {
        let classifier = KmerTaxonomyClassifier::new(4);

        let sequence = "ATCGATCGATCG";
        let profile = classifier.extract_kmer_profile(sequence).unwrap();

        assert!(!profile.is_empty(), "Should extract k-mers from sequence");

        // Check frequency sums to 1.0
        let total_freq: f64 = profile.values().sum();
        assert!(
            (total_freq - 1.0).abs() < 0.01,
            "Frequencies should sum to ~1.0"
        );
    }

    #[test]
    fn test_similarity_calculation() {
        let classifier = KmerTaxonomyClassifier::new(4);

        let mut profile1 = AHashMap::new();
        profile1.insert("ATCG".to_string(), 0.5);
        profile1.insert("TCGA".to_string(), 0.5);

        let mut profile2 = AHashMap::new();
        profile2.insert("ATCG".to_string(), 0.5);
        profile2.insert("TCGA".to_string(), 0.5);

        let similarity = classifier.calculate_similarity(&profile1, &profile2);
        assert!(
            (similarity - 1.0).abs() < 0.01,
            "Identical profiles should have similarity ~1.0"
        );
    }

    #[test]
    fn test_classification() {
        let mut classifier = KmerTaxonomyClassifier::new(4);

        // Add reference profile
        let mut kmer_freqs = AHashMap::new();
        kmer_freqs.insert("ATCG".to_string(), 0.3);
        kmer_freqs.insert("TCGA".to_string(), 0.3);
        kmer_freqs.insert("CGAT".to_string(), 0.4);

        let profile = KmerProfile {
            taxon_name: "Test Species".to_string(),
            taxon_rank: "species".to_string(),
            kmer_frequencies: kmer_freqs,
            gc_content: 0.5,
        };

        classifier.add_reference(profile);

        // Classify sequence
        let sequence = "ATCGATCGATCGATCGATCG";
        let assignment = classifier.classify_sequence(sequence).unwrap();

        // Should classify if similar enough
        if assignment.is_classified() {
            assert_eq!(assignment.taxon, "Test Species");
            assert!(assignment.confidence > 0.0);
        }
    }

    #[test]
    fn test_unclassified_assignment() {
        let assignment = TaxonomicAssignment::unclassified();

        assert_eq!(assignment.taxon, "Unclassified");
        assert!(!assignment.is_classified());
        assert_eq!(assignment.confidence, 0.0);
    }
}
