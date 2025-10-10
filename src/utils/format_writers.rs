/// Standard bioinformatics format writers for intermediate outputs
/// Provides FASTQ, GFA, and other standard format exports
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use tracing::info;

use crate::core::data_structures::{Contig, CorrectedRead};

/// Write corrected reads to FASTQ format (standard preprocessing output)
pub fn write_fastq<P: AsRef<Path>>(reads: &[CorrectedRead], output_path: P) -> Result<()> {
    let path = output_path.as_ref();
    info!(
        "🔍 write_fastq called with {} reads, writing to: {}",
        reads.len(),
        path.display()
    );

    if reads.is_empty() {
        info!("⚠️  WARNING: No reads to write!");
        return Ok(());
    }

    let file = File::create(path)
        .with_context(|| format!("Failed to create FASTQ file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    for (idx, read) in reads.iter().enumerate() {
        // @read_id
        writeln!(writer, "@read_{}", idx)?;

        // Sequence (use corrected sequence)
        writeln!(writer, "{}", read.corrected)?;

        // + separator
        writeln!(writer, "+")?;

        // Quality scores (generate dummy high-quality scores for corrected reads)
        // 'I' = Phred quality 40 (99.99% accuracy)
        let quality: String = std::iter::repeat_n('I', read.corrected.len()).collect();
        writeln!(writer, "{}", quality)?;
    }

    writer.flush()?;
    info!(
        "✅ Successfully wrote {} reads to FASTQ: {}",
        reads.len(),
        path.display()
    );
    Ok(())
}

/// Write assembly graph to GFA format (Graphical Fragment Assembly)
/// GFA is the standard format for assembly graphs used by MetaSPAdes, Bandage, etc.
pub fn write_gfa<P: AsRef<Path>>(contigs: &[Contig], output_path: P) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create GFA file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    // GFA Header
    writeln!(writer, "H\tVN:Z:1.0")?;

    // Segment lines (S) - one per contig
    for contig in contigs {
        writeln!(
            writer,
            "S\t{}\t{}\tLN:i:{}\tRC:i:{}\tDP:f:{:.2}",
            contig.id,
            contig.sequence,
            contig.length,
            (contig.coverage * contig.length as f64) as u32, // Read count estimate
            contig.coverage
        )?;
    }

    // Link lines (L) - connections between contigs
    // For now, we output a simple linear structure
    // In a full implementation, this would come from the actual graph edges
    for i in 0..(contigs.len().saturating_sub(1)) {
        writeln!(
            writer,
            "L\t{}\t+\t{}\t+\t0M",
            contigs[i].id,
            contigs[i + 1].id
        )?;
    }

    writer.flush()?;
    info!(
        "📊 Wrote {} segments to GFA: {}",
        contigs.len(),
        path.display()
    );
    Ok(())
}

/// Write contigs with detailed metadata to FASTA format
pub fn write_contigs_fasta_detailed<P: AsRef<Path>>(
    contigs: &[Contig],
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create FASTA file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    for contig in contigs {
        // Header with detailed metadata
        writeln!(
            writer,
            ">contig_{} length={} coverage={:.2}x gc={:.2}%",
            contig.id,
            contig.length,
            contig.coverage,
            calculate_gc_content(&contig.sequence) * 100.0
        )?;

        // Sequence (wrapped at 80 characters per line - FASTA standard)
        for chunk in contig.sequence.as_bytes().chunks(80) {
            writeln!(writer, "{}", std::str::from_utf8(chunk)?)?;
        }
    }

    writer.flush()?;
    info!(
        "📝 Wrote {} contigs to FASTA: {}",
        contigs.len(),
        path.display()
    );
    Ok(())
}

/// Write assembly statistics to human-readable text file
pub fn write_assembly_stats<P: AsRef<Path>>(
    contigs: &[Contig],
    n50: usize,
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create stats file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    let num_contigs = contigs.len();
    let total_length: usize = contigs.iter().map(|c| c.length).sum();
    let avg_length = if num_contigs > 0 {
        total_length / num_contigs
    } else {
        0
    };
    let avg_coverage: f64 =
        contigs.iter().map(|c| c.coverage).sum::<f64>() / num_contigs.max(1) as f64;
    let max_length = contigs.iter().map(|c| c.length).max().unwrap_or(0);
    let min_length = contigs.iter().map(|c| c.length).min().unwrap_or(0);

    let avg_gc: f64 = contigs
        .iter()
        .map(|c| calculate_gc_content(&c.sequence))
        .sum::<f64>()
        / num_contigs.max(1) as f64;

    writeln!(writer, "Assembly Statistics")?;
    writeln!(writer, "===================")?;
    writeln!(writer)?;
    writeln!(writer, "Contig Metrics:")?;
    writeln!(writer, "  Number of contigs: {}", num_contigs)?;
    writeln!(
        writer,
        "  Total assembly length: {} bp ({:.2} Mb)",
        total_length,
        total_length as f64 / 1_000_000.0
    )?;
    writeln!(writer, "  Average contig length: {} bp", avg_length)?;
    writeln!(writer, "  Longest contig: {} bp", max_length)?;
    writeln!(writer, "  Shortest contig: {} bp", min_length)?;
    writeln!(writer, "  N50: {} bp", n50)?;
    writeln!(writer)?;
    writeln!(writer, "Coverage Metrics:")?;
    writeln!(writer, "  Average coverage: {:.2}x", avg_coverage)?;
    writeln!(
        writer,
        "  Max coverage: {:.2}x",
        contigs.iter().map(|c| c.coverage).fold(0.0, f64::max)
    )?;
    writeln!(
        writer,
        "  Min coverage: {:.2}x",
        contigs
            .iter()
            .map(|c| c.coverage)
            .fold(f64::INFINITY, f64::min)
    )?;
    writeln!(writer)?;
    writeln!(writer, "Composition:")?;
    writeln!(writer, "  Average GC content: {:.2}%", avg_gc * 100.0)?;
    writeln!(writer)?;
    writeln!(writer, "Contig Length Distribution:")?;

    let mut length_bins = [0; 10];
    for contig in contigs {
        let bin = (contig.length / (max_length / 10 + 1)).min(9);
        length_bins[bin] += 1;
    }

    for (i, count) in length_bins.iter().enumerate() {
        if *count > 0 {
            let min_bin = i * max_length / 10;
            let max_bin = (i + 1) * max_length / 10;
            writeln!(writer, "  {}-{} bp: {} contigs", min_bin, max_bin, count)?;
        }
    }

    writer.flush()?;
    info!("📈 Wrote assembly statistics: {}", path.display());
    Ok(())
}

/// Write comprehensive quality control report with detailed statistics
pub fn write_qc_report<P: AsRef<Path>>(
    stats: &crate::qc::qc_stats::QCStats,
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create QC report: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    let retained_pct = if stats.reads_input > 0 {
        (stats.reads_passed as f64 / stats.reads_input as f64) * 100.0
    } else {
        0.0
    };

    let filtered_pct = 100.0 - retained_pct;
    let bases_retained_pct = if stats.total_bases_before > 0 {
        (stats.total_bases_after as f64 / stats.total_bases_before as f64) * 100.0
    } else {
        0.0
    };

    writeln!(
        writer,
        "╔══════════════════════════════════════════════════════════════════╗"
    )?;
    writeln!(
        writer,
        "║            QUALITY CONTROL REPORT - MetaForge                    ║"
    )?;
    writeln!(
        writer,
        "╚══════════════════════════════════════════════════════════════════╝"
    )?;
    writeln!(writer)?;

    // ═══ Input Statistics ═══
    writeln!(writer, "═══ INPUT STATISTICS ═══")?;
    writeln!(writer, "  Total reads:              {}", stats.reads_input)?;
    writeln!(
        writer,
        "  Total bases:              {} bp ({:.2} Mb)",
        stats.total_bases_before,
        stats.total_bases_before as f64 / 1_000_000.0
    )?;
    writeln!(
        writer,
        "  Average read length:      {:.1} bp",
        stats.mean_length_before
    )?;
    writeln!(
        writer,
        "  Average quality score:    Q{:.1}",
        stats.mean_quality_before
    )?;
    writeln!(
        writer,
        "  Q20 bases:                {:.2}%",
        stats.q20_percentage_before
    )?;
    writeln!(
        writer,
        "  Q30 bases:                {:.2}%",
        stats.q30_percentage_before
    )?;
    writeln!(writer)?;

    // ═══ Filtering Results ═══
    writeln!(writer, "═══ FILTERING RESULTS ═══")?;
    writeln!(
        writer,
        "  ✅ Reads passed:          {} ({:.2}%)",
        stats.reads_passed, retained_pct
    )?;
    writeln!(
        writer,
        "  ❌ Reads filtered:        {} ({:.2}%)",
        stats.reads_failed, filtered_pct
    )?;
    writeln!(writer)?;

    writeln!(writer, "  Failure breakdown:")?;
    writeln!(
        writer,
        "    • Low quality:          {} ({:.1}%)",
        stats.reads_failed_quality,
        (stats.reads_failed_quality as f64 / stats.reads_input.max(1) as f64) * 100.0
    )?;
    writeln!(
        writer,
        "    • Too short:            {} ({:.1}%)",
        stats.reads_failed_length,
        (stats.reads_failed_length as f64 / stats.reads_input.max(1) as f64) * 100.0
    )?;
    writeln!(
        writer,
        "    • Adapter issues:       {} ({:.1}%)",
        stats.reads_failed_adapter,
        (stats.reads_failed_adapter as f64 / stats.reads_input.max(1) as f64) * 100.0
    )?;
    writeln!(writer)?;

    // ═══ Adapter Detection ═══
    if stats.adapters_detected > 0 {
        writeln!(writer, "═══ ADAPTER DETECTION & TRIMMING ═══")?;
        writeln!(
            writer,
            "  Adapters detected:        {} ({:.1}% of input reads)",
            stats.adapters_detected,
            (stats.adapters_detected as f64 / stats.reads_input.max(1) as f64) * 100.0
        )?;
        writeln!(
            writer,
            "  Bases trimmed (adapter):  {} bp",
            stats.bases_trimmed_adapter
        )?;
        writeln!(writer)?;

        if !stats.adapter_types.is_empty() {
            writeln!(writer, "  Adapter types found:")?;
            let mut adapter_vec: Vec<_> = stats.adapter_types.iter().collect();
            adapter_vec.sort_by_key(|(_, count)| std::cmp::Reverse(*count));
            for (adapter, count) in adapter_vec {
                let adapter_name = if adapter.len() > 30 {
                    format!("{}...", &adapter[..27])
                } else {
                    adapter.clone()
                };
                writeln!(
                    writer,
                    "    • {}: {} occurrences ({:.1}%)",
                    adapter_name,
                    count,
                    (*count as f64 / stats.adapters_detected.max(1) as f64) * 100.0
                )?;
            }
            writeln!(writer)?;
        }
    }

    // ═══ Quality Trimming ═══
    writeln!(writer, "═══ QUALITY TRIMMING ═══")?;
    writeln!(
        writer,
        "  Bases trimmed (quality):  {} bp",
        stats.bases_trimmed_quality
    )?;
    writeln!(
        writer,
        "  Total bases trimmed:      {} bp",
        stats.bases_trimmed_quality + stats.bases_trimmed_adapter
    )?;
    writeln!(
        writer,
        "  Bases retained:           {} bp ({:.2}%)",
        stats.total_bases_after, bases_retained_pct
    )?;
    writeln!(writer)?;

    // ═══ Output Statistics ═══
    writeln!(writer, "═══ OUTPUT STATISTICS ═══")?;
    writeln!(writer, "  Total reads:              {}", stats.reads_passed)?;
    writeln!(
        writer,
        "  Total bases:              {} bp ({:.2} Mb)",
        stats.total_bases_after,
        stats.total_bases_after as f64 / 1_000_000.0
    )?;
    writeln!(
        writer,
        "  Average read length:      {:.1} bp",
        stats.mean_length_after
    )?;
    writeln!(
        writer,
        "  Average quality score:    Q{:.1}",
        stats.mean_quality_after
    )?;
    writeln!(
        writer,
        "  Q20 bases:                {:.2}%",
        stats.q20_percentage_after
    )?;
    writeln!(
        writer,
        "  Q30 bases:                {:.2}%",
        stats.q30_percentage_after
    )?;
    writeln!(writer)?;

    // ═══ Quality Improvement ═══
    let quality_improvement = stats.mean_quality_after - stats.mean_quality_before;
    let length_change = stats.mean_length_after - stats.mean_length_before;
    let q20_improvement = stats.q20_percentage_after - stats.q20_percentage_before;
    let q30_improvement = stats.q30_percentage_after - stats.q30_percentage_before;

    writeln!(writer, "═══ QUALITY IMPROVEMENT ═══")?;
    writeln!(
        writer,
        "  Quality score change:     {:+.1} (Q{:.1} → Q{:.1})",
        quality_improvement, stats.mean_quality_before, stats.mean_quality_after
    )?;
    writeln!(
        writer,
        "  Read length change:       {:+.1} bp ({:.1} → {:.1})",
        length_change, stats.mean_length_before, stats.mean_length_after
    )?;
    writeln!(
        writer,
        "  Q20 improvement:          {:+.2}% ({:.2}% → {:.2}%)",
        q20_improvement, stats.q20_percentage_before, stats.q20_percentage_after
    )?;
    writeln!(
        writer,
        "  Q30 improvement:          {:+.2}% ({:.2}% → {:.2}%)",
        q30_improvement, stats.q30_percentage_before, stats.q30_percentage_after
    )?;
    writeln!(writer)?;

    // ═══ Summary ═══
    writeln!(writer, "═══ SUMMARY ═══")?;
    let quality_status = if quality_improvement > 0.0 {
        "✅ IMPROVED"
    } else {
        "⚠️  UNCHANGED"
    };
    let retention_status = if retained_pct >= 80.0 {
        "✅ GOOD"
    } else if retained_pct >= 60.0 {
        "⚠️  MODERATE"
    } else {
        "❌ LOW"
    };

    writeln!(writer, "  Overall quality:          {}", quality_status)?;
    writeln!(
        writer,
        "  Read retention:           {} ({:.1}%)",
        retention_status, retained_pct
    )?;
    writeln!(
        writer,
        "  Data retained:            {:.1}% of input bases",
        bases_retained_pct
    )?;

    if stats.mean_quality_after >= 30.0 {
        writeln!(writer, "  📊 Output quality:        Excellent (Q30+)")?;
    } else if stats.mean_quality_after >= 20.0 {
        writeln!(writer, "  📊 Output quality:        Good (Q20+)")?;
    } else {
        writeln!(writer, "  📊 Output quality:        Needs review (<Q20)")?;
    }

    writeln!(writer)?;
    writeln!(
        writer,
        "╔══════════════════════════════════════════════════════════════════╗"
    )?;
    writeln!(
        writer,
        "║  Report generated by MetaForge Quality Control Pipeline          ║"
    )?;
    writeln!(
        writer,
        "╚══════════════════════════════════════════════════════════════════╝"
    )?;

    writer.flush()?;
    info!("✅ Wrote comprehensive QC report: {}", path.display());
    Ok(())
}

/// Calculate GC content of a sequence
fn calculate_gc_content(sequence: &str) -> f64 {
    let gc_count = sequence
        .chars()
        .filter(|&c| c == 'G' || c == 'C' || c == 'g' || c == 'c')
        .count();
    let total = sequence.len();

    if total > 0 {
        gc_count as f64 / total as f64
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::data_structures::{BaseCorrection, CorrectionMetadata};
    use crate::pipeline::complete_integration::TaxonomicClassification;

    #[test]
    fn test_gc_content_calculation() {
        assert_eq!(calculate_gc_content("ATCG"), 0.5);
        assert_eq!(calculate_gc_content("GGCC"), 1.0);
        assert_eq!(calculate_gc_content("AATT"), 0.0);
        assert_eq!(calculate_gc_content(""), 0.0);
    }

    #[test]
    fn test_fastq_write() {
        let reads = vec![CorrectedRead {
            id: 0,
            original: "ATCG".to_string(),
            corrected: "ATCG".to_string(),
            corrections: Vec::new(),
            quality_scores: vec![40, 40, 40, 40],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.9,
                context_window: 5,
                correction_time_ms: 0,
            },
            kmer_hash_cache: Vec::new(),
        }];

        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("test.fastq");

        write_fastq(&reads, &output_path).unwrap();

        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains("@read_0"));
        assert!(content.contains("ATCG"));
        assert!(content.contains("IIII")); // Quality scores
    }

    #[test]
    fn test_classification_tsv_write() {
        let classifications = vec![
            TaxonomicClassification {
                contig_id: 1,
                taxonomy_id: 100,
                taxonomy_name: "Escherichia coli".to_string(),
                confidence: 0.95,
                lineage: "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli".to_string(),
                method: "hybrid_kmer_taxonomy".to_string(),
            },
            TaxonomicClassification {
                contig_id: 2,
                taxonomy_id: 101,
                taxonomy_name: "Salmonella enterica".to_string(),
                confidence: 0.87,
                lineage: "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Salmonella;Salmonella enterica".to_string(),
                method: "ml_kmer_clustering".to_string(),
            },
        ];

        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("test_classifications.tsv");

        write_classification_tsv(&classifications, &output_path).unwrap();

        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains("contig_id\tbin_id\ttaxonomy_name")); // Header
        assert!(content.contains("Escherichia coli"));
        assert!(content.contains("Salmonella enterica"));
        assert!(content.contains("0.9500")); // Confidence format
        assert!(content.contains("hybrid_kmer_taxonomy"));
    }

    #[test]
    fn test_classification_summary_write() {
        let classifications = vec![
            TaxonomicClassification {
                contig_id: 1,
                taxonomy_id: 100,
                taxonomy_name: "Escherichia coli".to_string(),
                confidence: 0.95,
                lineage: "Bacteria;Proteobacteria".to_string(),
                method: "hybrid_kmer_taxonomy".to_string(),
            },
            TaxonomicClassification {
                contig_id: 2,
                taxonomy_id: 100,
                taxonomy_name: "Escherichia coli".to_string(),
                confidence: 0.92,
                lineage: "Bacteria;Proteobacteria".to_string(),
                method: "hybrid_kmer_taxonomy".to_string(),
            },
            TaxonomicClassification {
                contig_id: 3,
                taxonomy_id: 101,
                taxonomy_name: "Salmonella enterica".to_string(),
                confidence: 0.87,
                lineage: "Bacteria;Proteobacteria".to_string(),
                method: "ml_kmer_clustering".to_string(),
            },
        ];

        let temp_dir = tempfile::tempdir().unwrap();
        let output_path = temp_dir.path().join("test_classification_summary.txt");

        write_classification_summary(&classifications, &output_path).unwrap();

        let content = std::fs::read_to_string(&output_path).unwrap();
        assert!(content.contains("Total Contigs: 3"));
        assert!(content.contains("Taxa Distribution"));
        assert!(content.contains("Escherichia coli: 2 contigs (66.7%)")); // Most abundant
        assert!(content.contains("Salmonella enterica: 1 contigs (33.3%)"));
        assert!(content.contains("Classification Methods"));
    }
}

/// Write feature extraction results to TSV format
pub fn write_features_tsv<P: AsRef<Path>>(
    features: &[(String, Vec<f64>)], // (contig_id, feature_vector)
    feature_names: &[String],
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create features TSV: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    // Header row
    write!(writer, "contig_id")?;
    for name in feature_names {
        write!(writer, "\t{}", name)?;
    }
    writeln!(writer)?;

    // Data rows
    for (contig_id, feature_vec) in features {
        write!(writer, "{}", contig_id)?;
        for value in feature_vec {
            write!(writer, "\t{:.6}", value)?;
        }
        writeln!(writer)?;
    }

    writer.flush()?;
    info!(
        "📊 Wrote {} feature vectors to TSV: {}",
        features.len(),
        path.display()
    );
    Ok(())
}

/// Write feature extraction summary report
pub fn write_features_summary<P: AsRef<Path>>(
    num_sequences: usize,
    num_features: usize,
    feature_names: &[String],
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create features summary: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "Feature Extraction Summary")?;
    writeln!(writer, "=========================")?;
    writeln!(writer)?;
    writeln!(writer, "Sequences analyzed: {}", num_sequences)?;
    writeln!(writer, "Features per sequence: {}", num_features)?;
    writeln!(writer)?;
    writeln!(writer, "Feature Types:")?;

    for (i, name) in feature_names.iter().enumerate() {
        writeln!(writer, "  {}. {}", i + 1, name)?;
    }

    writer.flush()?;
    info!("📝 Wrote feature summary: {}", path.display());
    Ok(())
}

/// Write taxonomic classifications to TSV format (Kraken-compatible)
pub fn write_class_tsv<P: AsRef<Path>>(
    classifications: &[(String, String, String, f64, String)], // (contig_id, tax_id, tax_name, confidence, lineage)
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create classifications TSV: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    // Header
    writeln!(
        writer,
        "contig_id\ttaxonomy_id\ttaxonomy_name\tconfidence\tlineage"
    )?;

    // Data rows
    for (contig_id, tax_id, tax_name, confidence, lineage) in classifications {
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.4}\t{}",
            contig_id, tax_id, tax_name, confidence, lineage
        )?;
    }

    writer.flush()?;
    info!(
        "🏷️  Wrote {} classifications to TSV: {}",
        classifications.len(),
        path.display()
    );
    Ok(())
}

/// Write Kraken-style report (compatible with Krona)
pub fn write_kraken_report<P: AsRef<Path>>(
    classifications: &[(String, String, f64)], // (taxonomy_name, lineage, count)
    total_sequences: usize,
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create Kraken report: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    for (tax_name, _lineage, count) in classifications {
        let percentage = if total_sequences > 0 {
            (count / total_sequences as f64) * 100.0
        } else {
            0.0
        };

        writeln!(
            writer,
            "{:.2}\t{}\t{}\tS\t{}",
            percentage,
            *count as usize,
            *count as usize, // Rank (Species)
            tax_name
        )?;
    }

    writer.flush()?;
    info!("📊 Wrote Kraken report: {}", path.display());
    Ok(())
}

/// Write abundance profile to TSV format
pub fn write_abund_tsv<P: AsRef<Path>>(
    abundances: &[(String, f64, f64, usize)], // (taxon_name, relative_abundance, avg_coverage, num_contigs)
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create abundance TSV: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "taxon_name\trelative_abundance(%)\taverage_coverage\tnum_contigs"
    )?;

    for (taxon, rel_abund, avg_cov, num_contigs) in abundances {
        writeln!(
            writer,
            "{}\t{:.4}\t{:.2}\t{}",
            taxon, rel_abund, avg_cov, num_contigs
        )?;
    }

    writer.flush()?;
    info!(
        "📈 Wrote {} abundance entries to TSV: {}",
        abundances.len(),
        path.display()
    );
    Ok(())
}

/// Write Krona-compatible taxonomy file
pub fn write_krona_tax<P: AsRef<Path>>(
    abundances: &[(f64, String)], // (count, full_lineage)
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create Krona file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    for (count, lineage) in abundances {
        // Krona format: count<TAB>taxonomy1<TAB>taxonomy2<TAB>...
        write!(writer, "{:.0}", count)?;
        for taxon in lineage.split(';') {
            write!(writer, "\t{}", taxon.trim())?;
        }
        writeln!(writer)?;
    }

    writer.flush()?;
    info!("🌸 Wrote Krona taxonomy file: {}", path.display());
    Ok(())
}

/// Write abundance summary report
pub fn write_abund_summary<P: AsRef<Path>>(
    total_sequences: usize,
    num_taxa: usize,
    dominant_taxon: &str,
    dominant_abundance: f64,
    diversity_index: f64,
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create abundance summary: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "Abundance Profile Summary")?;
    writeln!(writer, "=========================")?;
    writeln!(writer)?;
    writeln!(writer, "Total Sequences: {}", total_sequences)?;
    writeln!(writer, "Unique Taxa: {}", num_taxa)?;
    writeln!(writer)?;
    writeln!(writer, "Dominant Taxon:")?;
    writeln!(writer, "  Name: {}", dominant_taxon)?;
    writeln!(writer, "  Relative Abundance: {:.2}%", dominant_abundance)?;
    writeln!(writer)?;
    writeln!(writer, "Diversity Metrics:")?;
    writeln!(writer, "  Shannon Diversity Index: {:.4}", diversity_index)?;
    writeln!(writer, "  (H' = 0: no diversity, H' > 3: high diversity)")?;

    writer.flush()?;
    info!("📊 Wrote abundance summary: {}", path.display());
    Ok(())
}

/// Write taxonomic classifications to TSV format (compatible with MetaBAT2/MaxBin2)
pub fn write_classification_tsv<P: AsRef<Path>>(
    classifications: &[crate::pipeline::complete_integration::TaxonomicClassification],
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    info!(
        "📊 Writing {} classifications to TSV: {}",
        classifications.len(),
        path.display()
    );

    if classifications.is_empty() {
        info!("⚠️  WARNING: No classifications to write!");
        return Ok(());
    }

    let file = File::create(path)
        .with_context(|| format!("Failed to create classification TSV: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    // TSV Header
    writeln!(
        writer,
        "contig_id\tbin_id\ttaxonomy_name\tconfidence\tlineage\tmethod"
    )?;

    // Write each classification
    for classification in classifications {
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.4}\t{}\t{}",
            classification.contig_id,
            classification.taxonomy_id,
            classification.taxonomy_name,
            classification.confidence,
            classification.lineage,
            classification.method
        )?;
    }

    writer.flush()?;
    info!(
        "✅ Successfully wrote {} classifications to: {}",
        classifications.len(),
        path.display()
    );
    Ok(())
}

/// Write human-readable classification summary
pub fn write_classification_summary<P: AsRef<Path>>(
    classifications: &[crate::pipeline::complete_integration::TaxonomicClassification],
    output_path: P,
) -> Result<()> {
    use std::collections::HashMap;

    let path = output_path.as_ref();
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Count taxa
    let mut taxa_counts: HashMap<String, usize> = HashMap::new();
    let mut method_counts: HashMap<String, usize> = HashMap::new();

    for c in classifications {
        *taxa_counts.entry(c.taxonomy_name.clone()).or_insert(0) += 1;
        *method_counts.entry(c.method.clone()).or_insert(0) += 1;
    }

    writeln!(writer, "# Classification Summary")?;
    writeln!(writer, "Total Contigs: {}\n", classifications.len())?;

    writeln!(writer, "## Taxa Distribution")?;
    let mut taxa_vec: Vec<_> = taxa_counts.into_iter().collect();
    taxa_vec.sort_by_key(|(_, count)| std::cmp::Reverse(*count));

    for (taxon, count) in taxa_vec {
        let pct = (count as f64 / classifications.len() as f64) * 100.0;
        writeln!(writer, "{}: {} contigs ({:.1}%)", taxon, count, pct)?;
    }

    writeln!(writer, "\n## Classification Methods")?;
    for (method, count) in method_counts {
        let pct = (count as f64 / classifications.len() as f64) * 100.0;
        writeln!(writer, "{}: {} contigs ({:.1}%)", method, count, pct)?;
    }

    writer.flush()?;
    info!("✅ Wrote classification summary to: {}", path.display());
    Ok(())
}

/// Write comprehensive analysis report in Markdown format
pub fn write_md_report<P: AsRef<Path>>(
    sample_name: &str,
    num_reads: usize,
    num_contigs: usize,
    assembly_n50: usize,
    num_taxa: usize,
    dominant_taxon: &str,
    processing_time_sec: f64,
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create Markdown report: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "# Metagenomics Analysis Report")?;
    writeln!(writer)?;
    writeln!(writer, "## Sample Information")?;
    writeln!(writer)?;
    writeln!(writer, "- **Sample Name**: {}", sample_name)?;
    writeln!(
        writer,
        "- **Analysis Date**: {}",
        chrono::Utc::now().format("%Y-%m-%d %H:%M:%S UTC")
    )?;
    writeln!(
        writer,
        "- **Processing Time**: {:.1} seconds",
        processing_time_sec
    )?;
    writeln!(writer)?;
    writeln!(writer, "## Input Data")?;
    writeln!(writer)?;
    writeln!(writer, "- **Number of Reads**: {}", num_reads)?;
    writeln!(writer)?;
    writeln!(writer, "## Assembly Results")?;
    writeln!(writer)?;
    writeln!(writer, "- **Number of Contigs**: {}", num_contigs)?;
    writeln!(writer, "- **N50**: {} bp", assembly_n50)?;
    writeln!(writer)?;
    writeln!(writer, "## Taxonomic Classification")?;
    writeln!(writer)?;
    writeln!(writer, "- **Unique Taxa Identified**: {}", num_taxa)?;
    writeln!(writer, "- **Dominant Organism**: {}", dominant_taxon)?;
    writeln!(writer)?;
    writeln!(writer, "## Output Files")?;
    writeln!(writer)?;
    writeln!(
        writer,
        "Use `ktImportText abundance/krona_input.txt -o krona_chart.html` for visualization."
    )?;

    writer.flush()?;
    info!("📄 Wrote Markdown report: {}", path.display());
    Ok(())
}

/// Write HTML analysis report
pub fn write_html_report<P: AsRef<Path>>(
    sample_name: &str,
    num_reads: usize,
    num_contigs: usize,
    assembly_n50: usize,
    num_taxa: usize,
    dominant_taxon: &str,
    dominant_abundance: f64,
    processing_time_sec: f64,
    output_path: P,
) -> Result<()> {
    let path = output_path.as_ref();
    let file = File::create(path)
        .with_context(|| format!("Failed to create HTML report: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "<!DOCTYPE html><html><head><title>{}</title>",
        sample_name
    )?;
    writeln!(writer, "<style>body{{font-family:Arial;margin:40px;background:#f5f5f5}}.container{{max-width:1200px;margin:0 auto;background:white;padding:30px;border-radius:8px}}</style></head><body><div class=\"container\">")?;
    writeln!(writer, "<h1>Metagenomics Analysis Report</h1>")?;
    writeln!(
        writer,
        "<p>Sample: {} | Reads: {} | Contigs: {} | N50: {} bp</p>",
        sample_name, num_reads, num_contigs, assembly_n50
    )?;
    writeln!(
        writer,
        "<p>Taxa: {} | Dominant: {} ({:.1}%)</p>",
        num_taxa, dominant_taxon, dominant_abundance
    )?;
    writeln!(
        writer,
        "<p>Processing Time: {:.1}s</p>",
        processing_time_sec
    )?;
    writeln!(writer, "</div></body></html>")?;

    writer.flush()?;
    info!("🌐 Wrote HTML report: {}", path.display());
    Ok(())
}
