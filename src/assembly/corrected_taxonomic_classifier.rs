//! Corrected Taxonomic Classification for Metagenomic Assembly
//!
//! This module implements biologically-sound taxonomic classification that:
//! 1. Groups contigs by sequence similarity rather than 1:1 assignment
//! 2. Uses k-mer-based LCA (Lowest Common Ancestor) classification
//! 3. Handles uncertain classifications appropriately
//! 4. Creates proper many-to-one contig:species relationships

use ahash::{AHashMap, AHashSet};
use anyhow::{anyhow, Result};
use std::collections::HashMap;

use crate::assembly::corrected_graph_construction::BiologicalContig;

/// Taxonomic classification result with confidence scoring
#[derive(Debug, Clone)]
pub struct TaxonomicClassification {
    pub contig_id: usize,
    pub taxonomy_id: u32,
    pub taxonomy_name: String,
    pub confidence: f64,
    pub lineage: String,
    pub method: String,
    pub supporting_kmers: usize,
    pub total_kmers: usize,
}

impl TaxonomicClassification {
    /// Create unclassified result
    pub fn unclassified(contig_id: usize) -> Self {
        Self {
            contig_id,
            taxonomy_id: 0,
            taxonomy_name: "Unclassified".to_string(),
            confidence: 0.0,
            lineage: "Unclassified".to_string(),
            method: "insufficient_evidence".to_string(),
            supporting_kmers: 0,
            total_kmers: 0,
        }
    }
    
    /// Check if classification meets confidence threshold
    pub fn is_confident(&self, threshold: f64) -> bool {
        self.confidence >= threshold && self.taxonomy_id != 0
    }
}

/// Taxonomic database entry
#[derive(Debug, Clone)]
pub struct TaxonomicReference {
    pub taxonomy_id: u32,
    pub taxonomy_name: String,
    pub lineage: String,
    pub rank: String,
    pub kmer_signatures: AHashSet<u64>,
}

/// K-mer based taxonomic classifier with LCA support
pub struct CorrectedTaxonomicClassifier {
    /// K-mer to taxonomy mapping
    kmer_taxonomy_db: AHashMap<u64, Vec<TaxonomicAssignment>>,
    /// Taxonomic hierarchy for LCA calculation
    taxonomy_tree: AHashMap<u32, TaxonomicReference>,
    /// Classification parameters
    confidence_threshold: f64,
    min_kmer_support: usize,
    k: usize,
}

#[derive(Debug, Clone)]
struct TaxonomicAssignment {
    taxonomy_id: u32,
    confidence: f64,
    source_organism: String,
}

impl CorrectedTaxonomicClassifier {
    /// Create new classifier with reference database
    pub fn new(k: usize, confidence_threshold: f64, min_kmer_support: usize) -> Self {
        Self {
            kmer_taxonomy_db: AHashMap::new(),
            taxonomy_tree: AHashMap::new(),
            confidence_threshold,
            min_kmer_support,
            k,
        }
    }
    
    /// Load taxonomic reference database (mock implementation for demo)
    pub fn load_reference_database(&mut self) -> Result<()> {
        println!("🔬 Loading taxonomic reference database");
        
        // Mock database with common gut microbiome species
        let references = vec![
            TaxonomicReference {
                taxonomy_id: 511145,
                taxonomy_name: "Escherichia coli str. K-12 substr. MG1655".to_string(),
                lineage: "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli".to_string(),
                rank: "strain".to_string(),
                kmer_signatures: self.generate_mock_ecoli_kmers(),
            },
            TaxonomicReference {
                taxonomy_id: 28901,
                taxonomy_name: "Salmonella enterica subsp. enterica serovar Typhimurium".to_string(),
                lineage: "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella enterica".to_string(),
                rank: "strain".to_string(),
                kmer_signatures: self.generate_mock_salmonella_kmers(),
            },
            TaxonomicReference {
                taxonomy_id: 1351,
                taxonomy_name: "Enterococcus faecalis".to_string(),
                lineage: "Bacteria;Firmicutes;Bacilli;Lactobacillales;Enterococcaceae;Enterococcus;Enterococcus faecalis".to_string(),
                rank: "species".to_string(),
                kmer_signatures: self.generate_mock_enterococcus_kmers(),
            },
            TaxonomicReference {
                taxonomy_id: 1680,
                taxonomy_name: "Bifidobacterium adolescentis".to_string(),
                lineage: "Bacteria;Actinobacteria;Actinobacteria;Bifidobacteriales;Bifidobacteriaceae;Bifidobacterium;Bifidobacterium adolescentis".to_string(),
                rank: "species".to_string(),
                kmer_signatures: self.generate_mock_bifidobacterium_kmers(),
            },
            TaxonomicReference {
                taxonomy_id: 853,
                taxonomy_name: "Faecalibacterium prausnitzii".to_string(),
                lineage: "Bacteria;Firmicutes;Clostridia;Clostridiales;Ruminococcaceae;Faecalibacterium;Faecalibacterium prausnitzii".to_string(),
                rank: "species".to_string(),
                kmer_signatures: self.generate_mock_faecalibacterium_kmers(),
            },
        ];
        
        // Build k-mer to taxonomy mapping
        for reference in references {
            for &kmer_hash in &reference.kmer_signatures {
                let assignment = TaxonomicAssignment {
                    taxonomy_id: reference.taxonomy_id,
                    confidence: 0.9, // High confidence for exact matches
                    source_organism: reference.taxonomy_name.clone(),
                };
                
                self.kmer_taxonomy_db
                    .entry(kmer_hash)
                    .or_insert_with(Vec::new)
                    .push(assignment);
            }
            
            self.taxonomy_tree.insert(reference.taxonomy_id, reference);
        }
        
        println!("   Loaded {} taxonomic references", self.taxonomy_tree.len());
        println!("   K-mer database: {} entries", self.kmer_taxonomy_db.len());
        
        Ok(())
    }
    
    /// Classify multiple contigs with proper species grouping
    pub fn classify_contigs(&self, contigs: &[BiologicalContig]) -> Result<Vec<TaxonomicClassification>> {
        println!("🔬 Classifying {} contigs using k-mer LCA approach", contigs.len());
        
        let mut classifications = Vec::new();
        let mut species_counts = AHashMap::new();
        
        for contig in contigs {
            let classification = self.classify_single_contig(contig)?;
            
            // Count species assignments
            if classification.is_confident(self.confidence_threshold) {
                *species_counts.entry(classification.taxonomy_id).or_insert(0) += 1;
            }
            
            classifications.push(classification);
        }
        
        // Log species diversity
        let confident_classifications = classifications.iter()
            .filter(|c| c.is_confident(self.confidence_threshold))
            .count();
        
        println!("   Confident classifications: {} / {} ({:.1}%)", 
                 confident_classifications, 
                 contigs.len(),
                 confident_classifications as f64 / contigs.len() as f64 * 100.0);
        
        println!("   Unique species detected: {}", species_counts.len());
        
        // Show species distribution
        let mut species_list: Vec<_> = species_counts.iter().collect();
        species_list.sort_by_key(|(_, count)| std::cmp::Reverse(**count));
        
        for (&taxonomy_id, &count) in species_list.iter().take(5) {
            if let Some(reference) = self.taxonomy_tree.get(&taxonomy_id) {
                println!("     {}: {} contigs", reference.taxonomy_name, count);
            }
        }
        
        Ok(classifications)
    }
    
    /// Classify single contig using k-mer LCA approach
    fn classify_single_contig(&self, contig: &BiologicalContig) -> Result<TaxonomicClassification> {
        let kmers = self.extract_kmers_from_contig(contig)?;
        
        if kmers.is_empty() {
            return Ok(TaxonomicClassification::unclassified(contig.id));
        }
        
        // Count taxonomic assignments for all k-mers
        let mut taxonomy_votes: AHashMap<u32, (usize, f64)> = AHashMap::new(); // (count, total_confidence)
        let mut classified_kmers = 0;
        
        for kmer_hash in &kmers {
            if let Some(assignments) = self.kmer_taxonomy_db.get(kmer_hash) {
                classified_kmers += 1;
                
                for assignment in assignments {
                    let entry = taxonomy_votes.entry(assignment.taxonomy_id).or_insert((0, 0.0));
                    entry.0 += 1;
                    entry.1 += assignment.confidence;
                }
            }
        }
        
        if classified_kmers < self.min_kmer_support {
            return Ok(TaxonomicClassification::unclassified(contig.id));
        }
        
        // Find best taxonomic assignment using LCA if needed
        let best_assignment = self.find_best_taxonomic_assignment(&taxonomy_votes)?;
        
        if let Some((taxonomy_id, (supporting_kmers, total_confidence))) = best_assignment {
            let confidence = total_confidence / supporting_kmers as f64;
            
            // Get taxonomic information
            let (taxonomy_name, lineage) = match self.taxonomy_tree.get(&taxonomy_id) {
                Some(reference) => (reference.taxonomy_name.clone(), reference.lineage.clone()),
                None => (format!("Unknown taxonomy {}", taxonomy_id), "Unknown".to_string()),
            };
            
            return Ok(TaxonomicClassification {
                contig_id: contig.id,
                taxonomy_id,
                taxonomy_name,
                confidence,
                lineage,
                method: "kmer_lca".to_string(),
                supporting_kmers,
                total_kmers: kmers.len(),
            });
        }
        
        Ok(TaxonomicClassification::unclassified(contig.id))
    }
    
    /// Extract k-mers from contig sequence
    fn extract_kmers_from_contig(&self, contig: &BiologicalContig) -> Result<Vec<u64>> {
        let mut kmers = Vec::new();
        let sequence_bytes = contig.sequence.as_bytes();
        
        for window in sequence_bytes.windows(self.k) {
            if let Ok(kmer_str) = std::str::from_utf8(window) {
                // Simple hash function for demo (would use proper canonical k-mer in production)
                let kmer_hash = self.hash_kmer(kmer_str);
                kmers.push(kmer_hash);
            }
        }
        
        Ok(kmers)
    }
    
    /// Simple k-mer hashing (would use CanonicalKmer in production)
    fn hash_kmer(&self, kmer: &str) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        let mut hasher = DefaultHasher::new();
        kmer.hash(&mut hasher);
        hasher.finish()
    }
    
    /// Find best taxonomic assignment considering confidence and support
    fn find_best_taxonomic_assignment(
        &self, 
        taxonomy_votes: &AHashMap<u32, (usize, f64)>
    ) -> Result<Option<(u32, (usize, f64))>> {
        if taxonomy_votes.is_empty() {
            return Ok(None);
        }
        
        // Find assignment with highest combined score (support count * average confidence)
        let mut best_assignment = None;
        let mut best_score = 0.0;
        
        for (&taxonomy_id, &(count, total_confidence)) in taxonomy_votes {
            let avg_confidence = total_confidence / count as f64;
            let combined_score = (count as f64) * avg_confidence;
            
            if combined_score > best_score {
                best_score = combined_score;
                best_assignment = Some((taxonomy_id, (count, total_confidence)));
            }
        }
        
        // Apply LCA if multiple high-scoring assignments exist
        if let Some(lca_assignment) = self.apply_lca_if_needed(taxonomy_votes) {
            return Ok(Some(lca_assignment));
        }
        
        Ok(best_assignment)
    }
    
    /// Apply Lowest Common Ancestor algorithm for ambiguous classifications
    fn apply_lca_if_needed(
        &self, 
        taxonomy_votes: &AHashMap<u32, (usize, f64)>
    ) -> Option<(u32, (usize, f64))> {
        // Find assignments with similar support levels
        let max_support = taxonomy_votes.values().map(|(count, _)| *count).max()?;
        let high_support_threshold = (max_support as f64 * 0.8) as usize;
        
        let high_support_taxa: Vec<u32> = taxonomy_votes
            .iter()
            .filter(|(_, (count, _))| *count >= high_support_threshold)
            .map(|(&taxonomy_id, _)| taxonomy_id)
            .collect();
        
        if high_support_taxa.len() > 1 {
            // Find LCA of high-support taxa
            if let Some(lca_taxonomy_id) = self.find_lowest_common_ancestor(&high_support_taxa) {
                let total_support: usize = high_support_taxa
                    .iter()
                    .filter_map(|&tid| taxonomy_votes.get(&tid))
                    .map(|(count, _)| count)
                    .sum();
                
                let total_confidence: f64 = high_support_taxa
                    .iter()
                    .filter_map(|&tid| taxonomy_votes.get(&tid))
                    .map(|(_, confidence)| confidence)
                    .sum();
                
                return Some((lca_taxonomy_id, (total_support, total_confidence)));
            }
        }
        
        None
    }
    
    /// Find lowest common ancestor in taxonomic tree
    fn find_lowest_common_ancestor(&self, taxonomy_ids: &[u32]) -> Option<u32> {
        if taxonomy_ids.is_empty() {
            return None;
        }
        
        if taxonomy_ids.len() == 1 {
            return Some(taxonomy_ids[0]);
        }
        
        // Get lineages for all taxonomies
        let lineages: Vec<Vec<&str>> = taxonomy_ids
            .iter()
            .filter_map(|&tid| self.taxonomy_tree.get(&tid))
            .map(|reference| reference.lineage.split(';').collect())
            .collect();
        
        if lineages.is_empty() {
            return None;
        }
        
        // Find common prefix
        let mut common_depth = 0;
        let min_depth = lineages.iter().map(|l| l.len()).min()?;
        
        for depth in 0..min_depth {
            let first_taxon = lineages[0][depth];
            if lineages.iter().all(|lineage| lineage[depth] == first_taxon) {
                common_depth = depth + 1;
            } else {
                break;
            }
        }
        
        // Find taxonomy ID corresponding to common ancestor
        if common_depth > 0 {
            let common_lineage = lineages[0][..common_depth].join(";");
            
            // Look for exact match in database
            for (taxonomy_id, reference) in &self.taxonomy_tree {
                if reference.lineage.starts_with(&common_lineage) && 
                   reference.lineage.matches(';').count() == common_depth - 1 {
                    return Some(*taxonomy_id);
                }
            }
        }
        
        None
    }
    
    // Mock k-mer generators for different species (for demonstration)
    fn generate_mock_ecoli_kmers(&self) -> AHashSet<u64> {
        let ecoli_sequences = vec![
            "ATGAAACGCATTAGCACCACCATTACCACCACCATC",
            "GGCATGCAAGCCGTTGAGTGGCGAACCTGCCACAG",
            "CGCATTAGCACCACCATTACCACCACCATCACCATTACCA",
        ];
        
        let mut kmers = AHashSet::new();
        for seq in ecoli_sequences {
            for window in seq.as_bytes().windows(self.k) {
                if let Ok(kmer_str) = std::str::from_utf8(window) {
                    kmers.insert(self.hash_kmer(kmer_str));
                }
            }
        }
        kmers
    }
    
    fn generate_mock_salmonella_kmers(&self) -> AHashSet<u64> {
        let salmonella_sequences = vec![
            "ATGAAACGCATTAGCACCACCATTACCGCCACCATC",
            "GGCATGCAAGCCGTTGAGTGGCGATCCTGCCACAG",
            "CGCATTAGCACCACCATTACCGCCACCATCGCCATTACCA",
        ];
        
        let mut kmers = AHashSet::new();
        for seq in salmonella_sequences {
            for window in seq.as_bytes().windows(self.k) {
                if let Ok(kmer_str) = std::str::from_utf8(window) {
                    kmers.insert(self.hash_kmer(kmer_str));
                }
            }
        }
        kmers
    }
    
    fn generate_mock_enterococcus_kmers(&self) -> AHashSet<u64> {
        let enterococcus_sequences = vec![
            "ATGAAATGCATTAGCACCACCATTACCACCACCATC",
            "GGCATGCAAACCGTTGAGTGGCGAACCTGCCACAG",
            "TGCATTAGCACCACCATTACCACCACCATCACCATTACCA",
        ];
        
        let mut kmers = AHashSet::new();
        for seq in enterococcus_sequences {
            for window in seq.as_bytes().windows(self.k) {
                if let Ok(kmer_str) = std::str::from_utf8(window) {
                    kmers.insert(self.hash_kmer(kmer_str));
                }
            }
        }
        kmers
    }
    
    fn generate_mock_bifidobacterium_kmers(&self) -> AHashSet<u64> {
        let bifidobacterium_sequences = vec![
            "ATGCCCGCATTAGCGCCACCATTACCGCCGCCATC",
            "GGCGTGCAAGCCGTTGCGTGGCGCACCTGCCGCAG",
            "GCATTAGCGCCACCATTACCGCCGCCATCGCCATTGCCA",
        ];
        
        let mut kmers = AHashSet::new();
        for seq in bifidobacterium_sequences {
            for window in seq.as_bytes().windows(self.k) {
                if let Ok(kmer_str) = std::str::from_utf8(window) {
                    kmers.insert(self.hash_kmer(kmer_str));
                }
            }
        }
        kmers
    }
    
    fn generate_mock_faecalibacterium_kmers(&self) -> AHashSet<u64> {
        let faecalibacterium_sequences = vec![
            "ATGAAATGCATTAGCACCACCATTACCACCACCATC",
            "GGCATGCAAACCGTTGAGTGGCGAACCTGCCACAG",
            "TGCATTAGCACCACCATTACCACCACCATCACCATTACCA",
        ];
        
        let mut kmers = AHashSet::new();
        for seq in faecalibacterium_sequences {
            for window in seq.as_bytes().windows(self.k) {
                if let Ok(kmer_str) = std::str::from_utf8(window) {
                    kmers.insert(self.hash_kmer(kmer_str));
                }
            }
        }
        kmers
    }
}

/// Generate species diversity report
pub fn generate_species_diversity_report(classifications: &[TaxonomicClassification]) -> SpeciesDiversityReport {
    let mut species_counts = HashMap::new();
    let mut genus_counts = HashMap::new();
    let mut confident_classifications = 0;
    
    for classification in classifications {
        if classification.is_confident(0.5) {
            confident_classifications += 1;
            
            // Count species
            *species_counts.entry(classification.taxonomy_name.clone()).or_insert(0) += 1;
            
            // Count genus (extract from lineage)
            if let Some(genus) = extract_genus_from_lineage(&classification.lineage) {
                *genus_counts.entry(genus).or_insert(0) += 1;
            }
        }
    }
    
    SpeciesDiversityReport {
        total_contigs: classifications.len(),
        classified_contigs: confident_classifications,
        unique_species: species_counts.len(),
        unique_genera: genus_counts.len(),
        species_distribution: species_counts,
        genus_distribution: genus_counts,
        classification_rate: confident_classifications as f64 / classifications.len() as f64,
    }
}

fn extract_genus_from_lineage(lineage: &str) -> Option<String> {
    let parts: Vec<&str> = lineage.split(';').collect();
    // Genus is typically second to last in bacterial taxonomy
    if parts.len() >= 2 {
        Some(parts[parts.len() - 2].to_string())
    } else {
        None
    }
}

#[derive(Debug)]
pub struct SpeciesDiversityReport {
    pub total_contigs: usize,
    pub classified_contigs: usize,
    pub unique_species: usize,
    pub unique_genera: usize,
    pub species_distribution: HashMap<String, usize>,
    pub genus_distribution: HashMap<String, usize>,
    pub classification_rate: f64,
}

impl SpeciesDiversityReport {
    pub fn print_summary(&self) {
        println!("\n🧬 Species Diversity Report");
        println!("═══════════════════════════");
        println!("Total contigs: {}", self.total_contigs);
        println!("Classified contigs: {} ({:.1}%)", 
                 self.classified_contigs, 
                 self.classification_rate * 100.0);
        println!("Unique species: {}", self.unique_species);
        println!("Unique genera: {}", self.unique_genera);
        
        println!("\n🏆 Top Species:");
        let mut species_list: Vec<_> = self.species_distribution.iter().collect();
        species_list.sort_by_key(|(_, count)| std::cmp::Reverse(**count));
        
        for (species, count) in species_list.iter().take(10) {
            let percentage = (*count as f64 / self.classified_contigs as f64) * 100.0;
            println!("  {}: {} contigs ({:.1}%)", species, count, percentage);
        }
        
        // Calculate diversity index (Shannon)
        let shannon_diversity = self.calculate_shannon_diversity();
        println!("\n📊 Shannon Diversity Index: {:.3}", shannon_diversity);
    }
    
    fn calculate_shannon_diversity(&self) -> f64 {
        if self.classified_contigs == 0 {
            return 0.0;
        }
        
        let mut shannon = 0.0;
        for &count in self.species_distribution.values() {
            let proportion = count as f64 / self.classified_contigs as f64;
            if proportion > 0.0 {
                shannon -= proportion * proportion.ln();
            }
        }
        shannon
    }
}