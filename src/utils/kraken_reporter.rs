use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;
use std::io::Write;
use chrono::{DateTime, Utc};

/// Kraken-style taxonomic classification result
/// Follows the standard Kraken output format for compatibility
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KrakenClassification {
    /// Classification status: 'C' for classified, 'U' for unclassified
    pub classification_status: ClassificationStatus,
    /// Sequence identifier (from FASTA/FASTQ header)
    pub sequence_id: String,
    /// NCBI Taxonomy ID (0 if unclassified)
    pub taxonomy_id: u32,
    /// Sequence length in base pairs
    pub sequence_length: usize,
    /// Space-delimited list of LCA mappings for each k-mer
    pub kmer_lca_mappings: Vec<u32>,
    /// Additional confidence metrics
    pub confidence_metrics: ConfidenceMetrics,
    /// Full taxonomic lineage
    pub taxonomic_lineage: TaxonomicLineage,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum ClassificationStatus {
    Classified,
    Unclassified,
}

impl ClassificationStatus {
    pub fn as_char(&self) -> char {
        match self {
            ClassificationStatus::Classified => 'C',
            ClassificationStatus::Unclassified => 'U',
        }
    }
}

/// Enhanced confidence metrics beyond basic Kraken format
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfidenceMetrics {
    /// Fraction of k-mers that support the classification
    pub kmer_support_fraction: f64,
    /// Number of k-mers classified
    pub classified_kmers: usize,
    /// Total number of k-mers in sequence
    pub total_kmers: usize,
    /// Confidence score (0.0 - 1.0)
    pub confidence_score: f64,
    /// Classification method used
    pub classification_method: String,
}

/// Hierarchical taxonomic lineage following NCBI taxonomy
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonomicLineage {
    /// Domain/Superkingdom (e.g., Bacteria, Archaea, Eukaryota)
    pub domain: Option<TaxonomicRank>,
    /// Kingdom
    pub kingdom: Option<TaxonomicRank>,
    /// Phylum
    pub phylum: Option<TaxonomicRank>,
    /// Class
    pub class: Option<TaxonomicRank>,
    /// Order
    pub order: Option<TaxonomicRank>,
    /// Family
    pub family: Option<TaxonomicRank>,
    /// Genus
    pub genus: Option<TaxonomicRank>,
    /// Species
    pub species: Option<TaxonomicRank>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonomicRank {
    pub taxonomy_id: u32,
    pub name: String,
    pub rank: String,
}

impl TaxonomicLineage {
    /// Convert to traditional pipe-delimited lineage string
    pub fn to_lineage_string(&self) -> String {
        let mut components = Vec::new();
        
        if let Some(domain) = &self.domain {
            components.push(domain.name.clone());
        }
        if let Some(kingdom) = &self.kingdom {
            components.push(kingdom.name.clone());
        }
        if let Some(phylum) = &self.phylum {
            components.push(phylum.name.clone());
        }
        if let Some(class) = &self.class {
            components.push(class.name.clone());
        }
        if let Some(order) = &self.order {
            components.push(order.name.clone());
        }
        if let Some(family) = &self.family {
            components.push(family.name.clone());
        }
        if let Some(genus) = &self.genus {
            components.push(genus.name.clone());
        }
        if let Some(species) = &self.species {
            components.push(species.name.clone());
        }
        
        components.join("|")
    }

    /// Get the most specific classification available
    pub fn get_most_specific(&self) -> Option<&TaxonomicRank> {
        [&self.species, &self.genus, &self.family, &self.order, 
         &self.class, &self.phylum, &self.kingdom, &self.domain]
            .iter()
            .find_map(|rank| rank.as_ref())
    }
}

/// Kraken-style report generator
pub struct KrakenReporter {
    /// Classifications to report
    classifications: Vec<KrakenClassification>,
    /// Sample name
    sample_name: String,
    /// Analysis timestamp
    timestamp: DateTime<Utc>,
    /// Summary statistics
    summary_stats: KrakenSummaryStats,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KrakenSummaryStats {
    pub total_sequences: usize,
    pub classified_sequences: usize,
    pub unclassified_sequences: usize,
    pub classification_rate: f64,
    pub unique_taxa: usize,
    pub total_kmers: usize,
    pub classified_kmers: usize,
}

impl KrakenReporter {
    /// Create a new Kraken reporter
    pub fn new(sample_name: String) -> Self {
        Self {
            classifications: Vec::new(),
            sample_name,
            timestamp: Utc::now(),
            summary_stats: KrakenSummaryStats {
                total_sequences: 0,
                classified_sequences: 0,
                unclassified_sequences: 0,
                classification_rate: 0.0,
                unique_taxa: 0,
                total_kmers: 0,
                classified_kmers: 0,
            },
        }
    }

    /// Add a classification result
    pub fn add_classification(&mut self, classification: KrakenClassification) {
        self.classifications.push(classification);
        self.update_summary_stats();
    }

    /// Update summary statistics
    fn update_summary_stats(&mut self) {
        self.summary_stats.total_sequences = self.classifications.len();
        self.summary_stats.classified_sequences = self.classifications
            .iter()
            .filter(|c| c.classification_status == ClassificationStatus::Classified)
            .count();
        self.summary_stats.unclassified_sequences = self.summary_stats.total_sequences - self.summary_stats.classified_sequences;
        
        self.summary_stats.classification_rate = if self.summary_stats.total_sequences > 0 {
            self.summary_stats.classified_sequences as f64 / self.summary_stats.total_sequences as f64
        } else {
            0.0
        };

        let unique_taxa: std::collections::HashSet<u32> = self.classifications
            .iter()
            .filter(|c| c.classification_status == ClassificationStatus::Classified)
            .map(|c| c.taxonomy_id)
            .collect();
        self.summary_stats.unique_taxa = unique_taxa.len();

        self.summary_stats.total_kmers = self.classifications
            .iter()
            .map(|c| c.confidence_metrics.total_kmers)
            .sum();
        self.summary_stats.classified_kmers = self.classifications
            .iter()
            .map(|c| c.confidence_metrics.classified_kmers)
            .sum();
    }

    /// Generate standard Kraken output file
    pub fn write_kraken_output(&self, output_path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create Kraken output file: {}", output_path.display()))?;

        for classification in &self.classifications {
            let kmer_mappings = classification.kmer_lca_mappings
                .iter()
                .map(|id| id.to_string())
                .collect::<Vec<_>>()
                .join(" ");

            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                classification.classification_status.as_char(),
                classification.sequence_id,
                classification.taxonomy_id,
                classification.sequence_length,
                kmer_mappings
            )?;
        }

        Ok(())
    }

    /// Generate Kraken-style report with abundance information
    pub fn write_kraken_report(&self, output_path: &Path) -> Result<()> {
        let mut file = std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create Kraken report file: {}", output_path.display()))?;

        // Write header
        writeln!(file, "# Kraken-style Metagenomic Classification Report")?;
        writeln!(file, "# Sample: {}", self.sample_name)?;
        writeln!(file, "# Timestamp: {}", self.timestamp.to_rfc3339())?;
        writeln!(file, "# Total sequences: {}", self.summary_stats.total_sequences)?;
        writeln!(file, "# Classified sequences: {} ({:.2}%)", 
                self.summary_stats.classified_sequences, 
                self.summary_stats.classification_rate * 100.0)?;
        writeln!(file, "# Unique taxa: {}", self.summary_stats.unique_taxa)?;
        writeln!(file)?;

        // Calculate abundance for each taxon
        let mut taxon_counts: HashMap<u32, TaxonAbundance> = HashMap::new();
        
        for classification in &self.classifications {
            if classification.classification_status == ClassificationStatus::Classified {
                let entry = taxon_counts
                    .entry(classification.taxonomy_id)
                    .or_insert_with(|| TaxonAbundance {
                        taxonomy_id: classification.taxonomy_id,
                        sequence_count: 0,
                        kmer_count: 0,
                        total_sequence_length: 0,
                        lineage: classification.taxonomic_lineage.clone(),
                        confidence_sum: 0.0,
                    });
                
                entry.sequence_count += 1;
                entry.kmer_count += classification.confidence_metrics.classified_kmers;
                entry.total_sequence_length += classification.sequence_length;
                entry.confidence_sum += classification.confidence_metrics.confidence_score;
            }
        }

        // Sort by sequence count (most abundant first)
        let mut sorted_taxa: Vec<_> = taxon_counts.into_values().collect();
        sorted_taxa.sort_by(|a, b| b.sequence_count.cmp(&a.sequence_count));

        // Write report header
        writeln!(file, "Percentage\tSequences\tK-mers\tTaxID\tRank\tName\tLineage")?;

        for taxon in &sorted_taxa {
            let percentage = (taxon.sequence_count as f64 / self.summary_stats.classified_sequences as f64) * 100.0;
            let avg_confidence = taxon.confidence_sum / taxon.sequence_count as f64;
            let most_specific = taxon.lineage.get_most_specific();
            
            let (rank, name) = if let Some(specific) = most_specific {
                (specific.rank.clone(), specific.name.clone())
            } else {
                ("unranked".to_string(), "Unknown".to_string())
            };

            writeln!(
                file,
                "{:.2}\t{}\t{}\t{}\t{}\t{}\t{}",
                percentage,
                taxon.sequence_count,
                taxon.kmer_count,
                taxon.taxonomy_id,
                rank,
                name,
                taxon.lineage.to_lineage_string()
            )?;
        }

        Ok(())
    }

    /// Generate enhanced JSON report with detailed classification information
    pub fn write_enhanced_json_report(&self, output_path: &Path) -> Result<()> {
        let report = EnhancedKrakenReport {
            sample_name: self.sample_name.clone(),
            timestamp: self.timestamp,
            summary_stats: self.summary_stats.clone(),
            classifications: self.classifications.clone(),
            taxonomic_summary: self.generate_taxonomic_summary(),
        };

        let json_content = serde_json::to_string_pretty(&report)
            .context("Failed to serialize enhanced Kraken report")?;
        
        std::fs::write(output_path, json_content)
            .with_context(|| format!("Failed to write enhanced JSON report: {}", output_path.display()))?;

        Ok(())
    }

    /// Generate taxonomic summary with hierarchical breakdown
    fn generate_taxonomic_summary(&self) -> TaxonomicSummary {
        let mut domain_counts = HashMap::new();
        let mut phylum_counts = HashMap::new();
        let mut class_counts = HashMap::new();
        let mut order_counts = HashMap::new();
        let mut family_counts = HashMap::new();
        let mut genus_counts = HashMap::new();
        let mut species_counts = HashMap::new();

        for classification in &self.classifications {
            if classification.classification_status == ClassificationStatus::Classified {
                let lineage = &classification.taxonomic_lineage;
                
                if let Some(domain) = &lineage.domain {
                    *domain_counts.entry(domain.name.clone()).or_insert(0) += 1;
                }
                if let Some(phylum) = &lineage.phylum {
                    *phylum_counts.entry(phylum.name.clone()).or_insert(0) += 1;
                }
                if let Some(class) = &lineage.class {
                    *class_counts.entry(class.name.clone()).or_insert(0) += 1;
                }
                if let Some(order) = &lineage.order {
                    *order_counts.entry(order.name.clone()).or_insert(0) += 1;
                }
                if let Some(family) = &lineage.family {
                    *family_counts.entry(family.name.clone()).or_insert(0) += 1;
                }
                if let Some(genus) = &lineage.genus {
                    *genus_counts.entry(genus.name.clone()).or_insert(0) += 1;
                }
                if let Some(species) = &lineage.species {
                    *species_counts.entry(species.name.clone()).or_insert(0) += 1;
                }
            }
        }

        TaxonomicSummary {
            domain_counts,
            phylum_counts,
            class_counts,
            order_counts,
            family_counts,
            genus_counts,
            species_counts,
        }
    }
}

#[derive(Debug, Clone)]
struct TaxonAbundance {
    taxonomy_id: u32,
    sequence_count: usize,
    kmer_count: usize,
    total_sequence_length: usize,
    lineage: TaxonomicLineage,
    confidence_sum: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnhancedKrakenReport {
    pub sample_name: String,
    pub timestamp: DateTime<Utc>,
    pub summary_stats: KrakenSummaryStats,
    pub classifications: Vec<KrakenClassification>,
    pub taxonomic_summary: TaxonomicSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonomicSummary {
    pub domain_counts: HashMap<String, usize>,
    pub phylum_counts: HashMap<String, usize>,
    pub class_counts: HashMap<String, usize>,
    pub order_counts: HashMap<String, usize>,
    pub family_counts: HashMap<String, usize>,
    pub genus_counts: HashMap<String, usize>,
    pub species_counts: HashMap<String, usize>,
}

/// Convert existing taxonomic classification to Kraken format
impl From<crate::pipeline::complete_integration::TaxonomicClassification> for KrakenClassification {
    fn from(classification: crate::pipeline::complete_integration::TaxonomicClassification) -> Self {
        // Parse the lineage string (assumes format like "Kingdom|Phylum|Class")
        let lineage_parts: Vec<&str> = classification.lineage.split('|').collect();
        
        let taxonomic_lineage = TaxonomicLineage {
            domain: None, // Would need to be parsed from lineage or database
            kingdom: lineage_parts.get(0).map(|name| TaxonomicRank {
                taxonomy_id: 0, // Would need lookup
                name: name.to_string(),
                rank: "kingdom".to_string(),
            }),
            phylum: lineage_parts.get(1).map(|name| TaxonomicRank {
                taxonomy_id: 0, // Would need lookup
                name: name.to_string(),
                rank: "phylum".to_string(),
            }),
            class: lineage_parts.get(2).map(|name| TaxonomicRank {
                taxonomy_id: 0, // Would need lookup
                name: name.to_string(),
                rank: "class".to_string(),
            }),
            order: None,
            family: None,
            genus: None,
            species: Some(TaxonomicRank {
                taxonomy_id: classification.taxonomy_id,
                name: classification.taxonomy_name.clone(),
                rank: "species".to_string(),
            }),
        };

        Self {
            classification_status: if classification.confidence > 0.0 {
                ClassificationStatus::Classified
            } else {
                ClassificationStatus::Unclassified
            },
            sequence_id: format!("contig_{}", classification.contig_id),
            taxonomy_id: classification.taxonomy_id,
            sequence_length: 1000, // Would need actual sequence length
            kmer_lca_mappings: vec![classification.taxonomy_id; 10], // Mock k-mer mappings
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: classification.confidence,
                classified_kmers: 100, // Mock value
                total_kmers: 100,      // Mock value
                confidence_score: classification.confidence,
                classification_method: classification.method,
            },
            taxonomic_lineage,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_kraken_classification_creation() {
        let classification = KrakenClassification {
            classification_status: ClassificationStatus::Classified,
            sequence_id: "test_seq_1".to_string(),
            taxonomy_id: 511145,
            sequence_length: 1000,
            kmer_lca_mappings: vec![511145, 511145, 562, 511145],
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: 0.85,
                classified_kmers: 85,
                total_kmers: 100,
                confidence_score: 0.85,
                classification_method: "exact_kmer_match".to_string(),
            },
            taxonomic_lineage: TaxonomicLineage {
                domain: Some(TaxonomicRank {
                    taxonomy_id: 2,
                    name: "Bacteria".to_string(),
                    rank: "superkingdom".to_string(),
                }),
                kingdom: None,
                phylum: Some(TaxonomicRank {
                    taxonomy_id: 1224,
                    name: "Proteobacteria".to_string(),
                    rank: "phylum".to_string(),
                }),
                class: Some(TaxonomicRank {
                    taxonomy_id: 1236,
                    name: "Gammaproteobacteria".to_string(),
                    rank: "class".to_string(),
                }),
                order: Some(TaxonomicRank {
                    taxonomy_id: 91347,
                    name: "Enterobacterales".to_string(),
                    rank: "order".to_string(),
                }),
                family: Some(TaxonomicRank {
                    taxonomy_id: 543,
                    name: "Enterobacteriaceae".to_string(),
                    rank: "family".to_string(),
                }),
                genus: Some(TaxonomicRank {
                    taxonomy_id: 561,
                    name: "Escherichia".to_string(),
                    rank: "genus".to_string(),
                }),
                species: Some(TaxonomicRank {
                    taxonomy_id: 511145,
                    name: "Escherichia coli str. K-12 substr. MG1655".to_string(),
                    rank: "strain".to_string(),
                }),
            },
        };

        assert_eq!(classification.classification_status, ClassificationStatus::Classified);
        assert_eq!(classification.taxonomy_id, 511145);
        
        let lineage_string = classification.taxonomic_lineage.to_lineage_string();
        assert!(lineage_string.contains("Bacteria"));
        assert!(lineage_string.contains("Proteobacteria"));
        assert!(lineage_string.contains("Escherichia"));
    }

    #[test]
    fn test_kraken_reporter() -> Result<()> {
        let temp_dir = tempdir()?;
        let mut reporter = KrakenReporter::new("test_sample".to_string());

        // Add test classifications
        let classification1 = KrakenClassification {
            classification_status: ClassificationStatus::Classified,
            sequence_id: "seq_1".to_string(),
            taxonomy_id: 562,
            sequence_length: 800,
            kmer_lca_mappings: vec![562, 562, 562],
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: 0.9,
                classified_kmers: 90,
                total_kmers: 100,
                confidence_score: 0.9,
                classification_method: "exact_match".to_string(),
            },
            taxonomic_lineage: TaxonomicLineage {
                domain: Some(TaxonomicRank {
                    taxonomy_id: 2,
                    name: "Bacteria".to_string(),
                    rank: "superkingdom".to_string(),
                }),
                kingdom: None,
                phylum: None,
                class: None,
                order: None,
                family: None,
                genus: Some(TaxonomicRank {
                    taxonomy_id: 561,
                    name: "Escherichia".to_string(),
                    rank: "genus".to_string(),
                }),
                species: Some(TaxonomicRank {
                    taxonomy_id: 562,
                    name: "Escherichia coli".to_string(),
                    rank: "species".to_string(),
                }),
            },
        };

        let classification2 = KrakenClassification {
            classification_status: ClassificationStatus::Unclassified,
            sequence_id: "seq_2".to_string(),
            taxonomy_id: 0,
            sequence_length: 600,
            kmer_lca_mappings: vec![],
            confidence_metrics: ConfidenceMetrics {
                kmer_support_fraction: 0.0,
                classified_kmers: 0,
                total_kmers: 80,
                confidence_score: 0.0,
                classification_method: "none".to_string(),
            },
            taxonomic_lineage: TaxonomicLineage {
                domain: None,
                kingdom: None,
                phylum: None,
                class: None,
                order: None,
                family: None,
                genus: None,
                species: None,
            },
        };

        reporter.add_classification(classification1);
        reporter.add_classification(classification2);

        // Test statistics
        assert_eq!(reporter.summary_stats.total_sequences, 2);
        assert_eq!(reporter.summary_stats.classified_sequences, 1);
        assert_eq!(reporter.summary_stats.unclassified_sequences, 1);
        assert_eq!(reporter.summary_stats.classification_rate, 0.5);

        // Test output generation
        let kraken_output_path = temp_dir.path().join("output.kraken");
        reporter.write_kraken_output(&kraken_output_path)?;
        assert!(kraken_output_path.exists());

        let report_path = temp_dir.path().join("report.txt");
        reporter.write_kraken_report(&report_path)?;
        assert!(report_path.exists());

        let json_path = temp_dir.path().join("enhanced.json");
        reporter.write_enhanced_json_report(&json_path)?;
        assert!(json_path.exists());

        Ok(())
    }

    #[test]
    fn test_lineage_parsing() {
        let lineage = TaxonomicLineage {
            domain: Some(TaxonomicRank {
                taxonomy_id: 2,
                name: "Bacteria".to_string(),
                rank: "superkingdom".to_string(),
            }),
            kingdom: None,
            phylum: Some(TaxonomicRank {
                taxonomy_id: 1224,
                name: "Proteobacteria".to_string(),
                rank: "phylum".to_string(),
            }),
            class: None,
            order: None,
            family: None,
            genus: Some(TaxonomicRank {
                taxonomy_id: 561,
                name: "Escherichia".to_string(),
                rank: "genus".to_string(),
            }),
            species: Some(TaxonomicRank {
                taxonomy_id: 562,
                name: "Escherichia coli".to_string(),
                rank: "species".to_string(),
            }),
        };

        let lineage_string = lineage.to_lineage_string();
        assert_eq!(lineage_string, "Bacteria|Proteobacteria|Escherichia|Escherichia coli");

        let most_specific = lineage.get_most_specific();
        assert!(most_specific.is_some());
        assert_eq!(most_specific.unwrap().name, "Escherichia coli");
    }
}