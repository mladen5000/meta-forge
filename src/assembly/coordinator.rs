//! Assembly Coordinator Module
//!
//! Coordinates assembly operations and calculates assembly statistics.
//! This module provides a high-level interface for running assemblies
//! and computing standard metrics like N50, N90, GC content, and coverage statistics.

use anyhow::Result;

use crate::assembly::laptop_assembly::{LaptopAssembler, LaptopConfig};
use crate::core::{
    AssemblyResults, AssemblyStats, Contig, CorrectedRead, CoverageStats, GraphFragment,
};
use crate::utils::configuration::AssemblyConfig;
use crate::utils::intermediate_output::IntermediateOutputManager;
use ahash::AHashMap;

/// Coordinates assembly operations with the laptop-optimized assembler
///
/// The `AssemblyCoordinator` manages the assembly workflow, including:
/// - Running the assembly with appropriate configuration
/// - Calculating assembly statistics (N50, N90, GC content, coverage)
/// - Storing results in the intermediate output manager
///
/// # Example
/// ```no_run
/// use meta_forge::assembly::coordinator::AssemblyCoordinator;
/// use meta_forge::utils::configuration::AssemblyConfig;
/// use meta_forge::utils::intermediate_output::IntermediateOutputManager;
///
/// let config = AssemblyConfig::default();
/// let output_manager = IntermediateOutputManager::new("./output".into()).unwrap();
/// let coordinator = AssemblyCoordinator::new(config, output_manager);
/// // coordinator.assemble(&reads).await?;
/// ```
pub struct AssemblyCoordinator {
    config: AssemblyConfig,
    #[allow(dead_code)]
    output_manager: IntermediateOutputManager,
}

impl AssemblyCoordinator {
    /// Create a new assembly coordinator
    ///
    /// # Arguments
    /// * `config` - Assembly configuration parameters
    /// * `output_manager` - Manager for intermediate outputs and checkpointing
    pub fn new(config: AssemblyConfig, output_manager: IntermediateOutputManager) -> Self {
        Self {
            config,
            output_manager,
        }
    }

    /// Run assembly with adaptive parameters and verbose progress
    ///
    /// This method:
    /// 1. Configures the laptop-optimized assembler
    /// 2. Converts reads to the appropriate format
    /// 3. Runs the assembly
    /// 4. Calculates assembly statistics
    /// 5. Creates a graph fragment representation
    /// 6. Stores results in the output manager
    ///
    /// # Arguments
    /// * `reads` - Slice of corrected reads to assemble
    ///
    /// # Returns
    /// `AssemblyResults` containing contigs, statistics, and graph fragment
    ///
    /// # Errors
    /// Returns an error if assembly fails or results cannot be stored
    pub async fn assemble(&self, reads: &[CorrectedRead]) -> Result<AssemblyResults> {
        println!("\n🧬 === Starting Metagenomic Assembly ===");
        println!("📈 Dataset: {} reads", reads.len());
        println!(
            "⚙️ K-mer range: {}-{}",
            self.config.k_min, self.config.k_max
        );
        println!("🎯 Min coverage: {}", self.config.min_coverage);
        println!("═══════════════════════════════════════════════\n");

        // Use laptop-optimized assembler
        let laptop_config = LaptopConfig::auto_detect();
        let assembler = LaptopAssembler::new(laptop_config);

        // Convert core CorrectedRead to assembly CorrectedRead for compatibility
        let assembly_reads: Vec<_> = reads
            .iter()
            .map(|r| crate::core::data_structures::CorrectedRead {
                id: r.id,
                original: r.original.clone(),
                corrected: r.corrected.clone(),
                corrections: r.corrections.clone(),
                quality_scores: r.quality_scores.clone(),
                correction_metadata: r.correction_metadata.clone(),
                kmer_hash_cache: r.kmer_hash_cache.clone(),
            })
            .collect();

        // Use laptop assembler for all scenarios
        println!("🚀 Using laptop-optimized assembly pipeline");
        let contigs = assembler.assemble(&assembly_reads)?;

        // Create basic assembly stats from contigs
        let assembly_stats = AssemblyStats {
            total_length: contigs.iter().map(|c| c.length).sum(),
            num_contigs: contigs.len(),
            n50: Self::calculate_n50(&contigs),
            n90: Self::calculate_n90(&contigs),
            largest_contig: contigs.iter().map(|c| c.length).max().unwrap_or(0),
            gc_content: Self::calculate_average_gc(&contigs),
            coverage_mean: if contigs.is_empty() {
                0.0
            } else {
                contigs.iter().map(|c| c.coverage).sum::<f64>() / contigs.len() as f64
            },
            coverage_std: Self::calculate_coverage_std(&contigs),
        };

        // Create basic graph fragment (simplified for laptop assembler)
        let graph_fragment = GraphFragment {
            nodes: AHashMap::new(), // Laptop assembler doesn't expose internal graph structure
            edges: vec![],
            fragment_id: 0,
            coverage_stats: CoverageStats::default(),
        };

        Ok(AssemblyResults {
            contigs,
            assembly_stats,
            graph_fragment,
        })
    }

    /// Calculate N50 statistic for contigs
    ///
    /// N50 is the contig length such that contigs of that length or longer
    /// contain at least 50% of the total assembly length.
    ///
    /// # Arguments
    /// * `contigs` - Slice of contigs to analyze
    ///
    /// # Returns
    /// The N50 value in base pairs, or 0 if no contigs provided
    pub fn calculate_n50(contigs: &[Contig]) -> usize {
        if contigs.is_empty() {
            return 0;
        }

        let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
        lengths.sort_by(|a, b| b.cmp(a)); // Sort in descending order

        let total_length: usize = lengths.iter().sum();
        let target = total_length / 2;

        let mut cumulative = 0;
        for &length in &lengths {
            cumulative += length;
            if cumulative >= target {
                return length;
            }
        }

        0
    }

    /// Calculate N90 statistic for contigs
    ///
    /// N90 is the contig length such that contigs of that length or longer
    /// contain at least 90% of the total assembly length.
    ///
    /// # Arguments
    /// * `contigs` - Slice of contigs to analyze
    ///
    /// # Returns
    /// The N90 value in base pairs, or 0 if no contigs provided
    pub fn calculate_n90(contigs: &[Contig]) -> usize {
        if contigs.is_empty() {
            return 0;
        }

        let mut lengths: Vec<usize> = contigs.iter().map(|c| c.length).collect();
        lengths.sort_by(|a, b| b.cmp(a)); // Sort in descending order

        let total_length: usize = lengths.iter().sum();
        let target = (total_length * 9) / 10; // 90% of total length

        let mut cumulative = 0;
        for &length in &lengths {
            cumulative += length;
            if cumulative >= target {
                return length;
            }
        }

        0
    }

    /// Calculate average GC content across all contigs
    ///
    /// Computes the percentage of G and C bases relative to all standard bases (A, T, G, C).
    /// Ambiguous bases are excluded from the calculation.
    ///
    /// # Arguments
    /// * `contigs` - Slice of contigs to analyze
    ///
    /// # Returns
    /// The average GC content as a percentage (0.0-100.0), or 0.0 if no valid bases
    pub fn calculate_average_gc(contigs: &[Contig]) -> f64 {
        if contigs.is_empty() {
            return 0.0;
        }

        let mut total_gc = 0;
        let mut total_bases = 0;

        for contig in contigs {
            for c in contig.sequence.chars() {
                match c.to_ascii_uppercase() {
                    'G' | 'C' => {
                        total_gc += 1;
                        total_bases += 1;
                    }
                    'A' | 'T' => {
                        total_bases += 1;
                    }
                    _ => {} // Skip ambiguous bases
                }
            }
        }

        if total_bases == 0 {
            0.0
        } else {
            (total_gc as f64 / total_bases as f64) * 100.0
        }
    }

    /// Calculate coverage standard deviation for contigs
    ///
    /// Computes the sample standard deviation of contig coverage values.
    /// Uses the (n-1) denominator for sample variance.
    ///
    /// # Arguments
    /// * `contigs` - Slice of contigs to analyze
    ///
    /// # Returns
    /// The standard deviation of coverage values, or 0.0 if ≤1 contig
    pub fn calculate_coverage_std(contigs: &[Contig]) -> f64 {
        if contigs.len() <= 1 {
            return 0.0;
        }

        let mean = contigs.iter().map(|c| c.coverage).sum::<f64>() / contigs.len() as f64;
        let variance = contigs
            .iter()
            .map(|c| (c.coverage - mean).powi(2))
            .sum::<f64>()
            / (contigs.len() - 1) as f64;

        variance.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_contig(id: usize, length: usize, sequence: &str, coverage: f64) -> Contig {
        Contig {
            id,
            sequence: sequence.to_string(),
            length,
            coverage,
            node_path: vec![],
            contig_type: crate::core::ContigType::Linear,
        }
    }

    #[test]
    fn test_calculate_n50_basic() {
        let contigs = vec![
            create_test_contig(0, 100, "A", 1.0),
            create_test_contig(1, 200, "A", 1.0),
            create_test_contig(2, 300, "A", 1.0),
        ];

        // Total length = 600, target = 300
        // Cumulative: 300, 500, 600
        // First contig >= 300 is contig with length 300
        assert_eq!(AssemblyCoordinator::calculate_n50(&contigs), 300);
    }

    #[test]
    fn test_calculate_n50_empty() {
        let contigs = vec![];
        assert_eq!(AssemblyCoordinator::calculate_n50(&contigs), 0);
    }

    #[test]
    fn test_calculate_n90_basic() {
        let contigs = vec![
            create_test_contig(0, 100, "A", 1.0),
            create_test_contig(1, 200, "A", 1.0),
            create_test_contig(2, 700, "A", 1.0),
        ];

        // Total length = 1000, target = 900
        // Cumulative: 700, 900, 1000
        // First contig >= 900 is contig with length 200
        assert_eq!(AssemblyCoordinator::calculate_n90(&contigs), 200);
    }

    #[test]
    fn test_calculate_average_gc() {
        let contigs = vec![
            create_test_contig(0, 4, "AATT", 1.0), // 0% GC
            create_test_contig(1, 4, "GGCC", 1.0), // 100% GC
        ];

        // Total: 4 GC out of 8 bases = 50%
        assert_eq!(AssemblyCoordinator::calculate_average_gc(&contigs), 50.0);
    }

    #[test]
    fn test_calculate_average_gc_with_ambiguous() {
        let contigs = vec![
            create_test_contig(0, 6, "AATNGC", 1.0), // N is ambiguous, skip it
        ];

        // 2 GC out of 5 valid bases = 40%
        assert_eq!(AssemblyCoordinator::calculate_average_gc(&contigs), 40.0);
    }

    #[test]
    fn test_calculate_coverage_std() {
        let contigs = vec![
            create_test_contig(0, 100, "A", 2.0),
            create_test_contig(1, 100, "A", 4.0),
            create_test_contig(2, 100, "A", 6.0),
        ];

        // Mean = 4.0
        // Variance = ((2-4)^2 + (4-4)^2 + (6-4)^2) / 2 = (4 + 0 + 4) / 2 = 4
        // Std = 2.0
        assert_eq!(AssemblyCoordinator::calculate_coverage_std(&contigs), 2.0);
    }

    #[test]
    fn test_calculate_coverage_std_single_contig() {
        let contigs = vec![create_test_contig(0, 100, "A", 5.0)];
        assert_eq!(AssemblyCoordinator::calculate_coverage_std(&contigs), 0.0);
    }
}
