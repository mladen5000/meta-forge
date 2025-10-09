//! Assembly orchestrator for pipeline-level assembly coordination
//!
//! Provides high-level assembly operations including:
//! - Configuration management
//! - Database integration
//! - Statistics calculation
//! - Coverage-based binning

use anyhow::Result;
use ahash::AHashMap;
use tracing::info;

use crate::core::data_structures::{
    AssemblyStats, Contig, CorrectedRead, CoverageStats, GraphFragment,
};
use crate::pipeline::complete_integration::AssemblyResults;
use crate::database::integration::MetagenomicsDatabase;
use super::{LaptopAssembler, LaptopConfig};

/// Assembly orchestrator for managing the assembly process
pub struct AssemblyOrchestrator {
    config: LaptopConfig,
    database: Option<MetagenomicsDatabase>,
}

impl AssemblyOrchestrator {
    /// Create a new assembly orchestrator with auto-detected configuration
    pub fn new() -> Self {
        Self {
            config: LaptopConfig::auto_detect(),
            database: None,
        }
    }

    /// Create with custom configuration
    pub fn with_config(config: LaptopConfig) -> Self {
        Self {
            config,
            database: None,
        }
    }

    /// Attach a database for storing results
    pub fn with_database(mut self, database: MetagenomicsDatabase) -> Self {
        self.database = Some(database);
        self
    }

    /// Run assembly with verbose progress reporting
    pub fn assemble(&self, reads: &[CorrectedRead]) -> Result<AssemblyResults> {
        println!("\n🧬 === Starting Metagenomic Assembly ===");
        println!("📈 Dataset: {} reads", reads.len());
        println!("💻 Using configuration: {} MB RAM, {} CPU cores",
            self.config.memory_budget_mb, self.config.cpu_cores);
        println!("═══════════════════════════════════════════════\n");

        // Use laptop-optimized assembler
        let assembler = LaptopAssembler::new(self.config.clone());

        info!("🚀 Using laptop-optimized assembly pipeline");
        let contigs = assembler.assemble(reads)?;

        // Calculate assembly statistics
        let assembly_stats = Self::calculate_assembly_stats(&contigs);

        // Create basic graph fragment (simplified for laptop assembler)
        let graph_fragment = GraphFragment {
            nodes: AHashMap::new(),
            edges: vec![],
            fragment_id: 0,
            coverage_stats: CoverageStats::default(),
        };

        let results = AssemblyResults {
            contigs,
            assembly_stats,
            graph_fragment,
        };

        // Store in database if available
        if let Some(ref db) = self.database {
            // Create a simple config string for database storage
            let config_string = format!(
                "memory_budget_mb: {}, cpu_cores: {}, max_k: {}",
                self.config.memory_budget_mb,
                self.config.cpu_cores,
                self.config.max_k
            );

            let assembly_id = db.store_assembly_results(
                &results.assembly_stats,
                "current_sample",
                &config_string,
            )?;

            db.store_contigs(assembly_id, &results.contigs)?;
        }

        Ok(results)
    }

    /// Calculate comprehensive assembly statistics
    pub fn calculate_assembly_stats(contigs: &[Contig]) -> AssemblyStats {
        AssemblyStats {
            total_length: contigs.iter().map(|c| c.length).sum(),
            num_contigs: contigs.len(),
            n50: Self::calculate_n50(contigs),
            n90: Self::calculate_n90(contigs),
            largest_contig: contigs.iter().map(|c| c.length).max().unwrap_or(0),
            gc_content: Self::calculate_average_gc(contigs),
            coverage_mean: if contigs.is_empty() {
                0.0
            } else {
                contigs.iter().map(|c| c.coverage).sum::<f64>() / contigs.len() as f64
            },
            coverage_std: Self::calculate_coverage_std(contigs),
        }
    }

    /// Calculate N50 statistic
    pub fn calculate_n50(contigs: &[Contig]) -> usize {
        if contigs.is_empty() {
            return 0;
        }

        let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
        lengths.sort_unstable_by(|a, b| b.cmp(a)); // Sort descending

        let total_length: usize = lengths.iter().sum();
        let half_length = total_length / 2;

        let mut cumulative = 0;
        for &length in &lengths {
            cumulative += length;
            if cumulative >= half_length {
                return length;
            }
        }

        0
    }

    /// Calculate N90 statistic
    pub fn calculate_n90(contigs: &[Contig]) -> usize {
        if contigs.is_empty() {
            return 0;
        }

        let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
        lengths.sort_unstable_by(|a, b| b.cmp(a));

        let total_length: usize = lengths.iter().sum();
        let ninety_percent = (total_length * 9) / 10;

        let mut cumulative = 0;
        for &length in &lengths {
            cumulative += length;
            if cumulative >= ninety_percent {
                return length;
            }
        }

        0
    }

    /// Calculate average GC content across all contigs
    pub fn calculate_average_gc(contigs: &[Contig]) -> f64 {
        if contigs.is_empty() {
            return 0.0;
        }

        let total_gc: usize = contigs
            .iter()
            .map(|c| {
                c.sequence
                    .chars()
                    .filter(|&ch| ch == 'G' || ch == 'C' || ch == 'g' || ch == 'c')
                    .count()
            })
            .sum();

        let total_bases: usize = contigs.iter().map(|c| c.length).sum();

        if total_bases == 0 {
            0.0
        } else {
            (total_gc as f64 / total_bases as f64) * 100.0
        }
    }

    /// Calculate standard deviation of coverage
    pub fn calculate_coverage_std(contigs: &[Contig]) -> f64 {
        if contigs.len() <= 1 {
            return 0.0;
        }

        let mean: f64 = contigs.iter().map(|c| c.coverage).sum::<f64>() / contigs.len() as f64;

        let variance: f64 = contigs
            .iter()
            .map(|c| {
                let diff = c.coverage - mean;
                diff * diff
            })
            .sum::<f64>()
            / (contigs.len() - 1) as f64;

        variance.sqrt()
    }

    /// Bin contigs by coverage uniformity for species identification
    ///
    /// Based on STRONG (Quince et al., 2021) methodology:
    /// Contigs from the same genome have similar sequencing depth
    pub fn bin_contigs_by_coverage(contigs: &[Contig]) -> Vec<Vec<usize>> {
        if contigs.is_empty() {
            return Vec::new();
        }

        // Calculate median coverage
        let mut coverages: Vec<f64> = contigs.iter().map(|c| c.coverage).collect();
        coverages.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_coverage = coverages[coverages.len() / 2];

        // Group contigs by coverage similarity (±50% tolerance)
        let coverage_tolerance = 0.5; // ±50%
        let min_coverage = median_coverage * (1.0 - coverage_tolerance);
        let max_coverage = median_coverage * (1.0 + coverage_tolerance);

        let mut main_bin: Vec<usize> = Vec::new();
        let mut outlier_bins: Vec<Vec<usize>> = Vec::new();

        for (idx, contig) in contigs.iter().enumerate() {
            if contig.coverage >= min_coverage && contig.coverage <= max_coverage {
                // Main genome bin (uniform coverage)
                main_bin.push(idx);
            } else {
                // Potential contamination, plasmid, or assembly error
                // Group outliers by coverage (e.g., plasmids at different copy number)
                let mut placed = false;
                for outlier_bin in outlier_bins.iter_mut() {
                    if let Some(&first_idx) = outlier_bin.first() {
                        let bin_coverage = contigs[first_idx].coverage;
                        let ratio = contig.coverage / bin_coverage;
                        if (0.7..=1.3).contains(&ratio) {
                            outlier_bin.push(idx);
                            placed = true;
                            break;
                        }
                    }
                }
                if !placed {
                    outlier_bins.push(vec![idx]);
                }
            }
        }

        let mut bins = Vec::new();
        if !main_bin.is_empty() {
            bins.push(main_bin);
        }
        bins.extend(outlier_bins);

        bins
    }
}

impl Default for AssemblyOrchestrator {
    fn default() -> Self {
        Self::new()
    }
}
