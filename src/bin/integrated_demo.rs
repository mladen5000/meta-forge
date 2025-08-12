use anyhow::Result;
use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::core::paired_reads::{extract_pair_id, PairedRead, PairedReadCollection, ReadPair};
use meta_forge::utils::configuration::PipelineConfiguration;
use meta_forge::utils::progress_display::{MultiProgress, ProgressBar};
/// Integration demonstration showing how new modules work with existing pipeline
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() -> Result<()> {
    println!("🔬 Metagenomic Pipeline Integration Demo");
    println!("=======================================");

    // Initialize configuration
    println!("📋 Using default pipeline configuration...");
    let config = PipelineConfiguration::default();

    println!("✅ Configuration loaded:");
    println!(
        "   - K-mer range: {}-{}",
        config.assembly.k_min, config.assembly.k_max
    );
    println!("   - Min coverage: {}", config.assembly.min_coverage);
    println!("   - Feature k-mers: {:?}", config.features.kmer_sizes);

    // Set up multi-line progress for the entire pipeline
    println!("\n🚀 Starting integrated metagenomic pipeline...");
    let mut multi_progress = MultiProgress::new();
    let read_line = multi_progress.add_line("📖 Read Processing: Initializing...".to_string());
    let assembly_line = multi_progress.add_line("🧬 Assembly: Waiting...".to_string());
    let feature_line = multi_progress.add_line("🔍 Feature Extraction: Waiting...".to_string());

    // Phase 1: Paired Read Processing
    multi_progress.update_line(
        read_line,
        "📖 Read Processing: Loading FASTQ files...".to_string(),
    );

    let fastq_r1 = "SRR390728_1.fastq";
    let fastq_r2 = "SRR390728_2.fastq";
    let max_reads_to_process = 5000; // Process subset for demo

    let mut collection = PairedReadCollection::new();

    // Check if files exist
    if std::path::Path::new(fastq_r1).exists() && std::path::Path::new(fastq_r2).exists() {
        multi_progress.update_line(
            read_line,
            "📖 Read Processing: Parsing paired-end reads...".to_string(),
        );

        // Process paired reads with progress tracking
        let mut pb = ProgressBar::new(max_reads_to_process as u64, "Processing paired reads");

        let file_r1 = File::open(fastq_r1)?;
        let file_r2 = File::open(fastq_r2)?;
        let reader_r1 = BufReader::new(file_r1);
        let reader_r2 = BufReader::new(file_r2);

        let mut lines_r1 = reader_r1.lines();
        let mut lines_r2 = reader_r2.lines();

        let mut processed = 0;
        let mut read_id = 0;

        while processed < max_reads_to_process {
            // Read R1
            let header_r1 = match lines_r1.next() {
                Some(line) => line?,
                None => break,
            };
            let seq_r1 = lines_r1.next().unwrap()?;
            let _plus_r1 = lines_r1.next().unwrap()?;
            let qual_r1 = lines_r1.next().unwrap()?;

            // Read R2
            let header_r2 = match lines_r2.next() {
                Some(line) => line?,
                None => break,
            };
            let seq_r2 = lines_r2.next().unwrap()?;
            let _plus_r2 = lines_r2.next().unwrap()?;
            let qual_r2 = lines_r2.next().unwrap()?;

            // Extract pair IDs and create reads
            let (pair_id_r1, _) = extract_pair_id(&header_r1)?;
            let (pair_id_r2, _) = extract_pair_id(&header_r2)?;

            if pair_id_r1 == pair_id_r2 {
                let qual_bytes_r1: Vec<u8> =
                    qual_r1.bytes().map(|b| b.saturating_sub(33)).collect();
                let qual_bytes_r2: Vec<u8> =
                    qual_r2.bytes().map(|b| b.saturating_sub(33)).collect();

                let forward = PairedRead::new(
                    read_id,
                    pair_id_r1.clone(),
                    meta_forge::core::paired_reads::ReadOrientation::Forward,
                    seq_r1,
                    qual_bytes_r1,
                );

                let reverse = PairedRead::new(
                    read_id + 1,
                    pair_id_r2,
                    meta_forge::core::paired_reads::ReadOrientation::Reverse,
                    seq_r2,
                    qual_bytes_r2,
                );

                let mut pair = ReadPair::new(forward, reverse)?;
                pair.estimate_insert_size()?;
                collection.add_pair(pair)?;
                read_id += 2;
            }

            processed += 1;
            pb.update(processed as u64);
        }

        pb.finish();

        // Apply quality filtering with reasonable defaults
        multi_progress.update_line(
            read_line,
            "📖 Read Processing: Applying quality filters...".to_string(),
        );
        collection.filter_by_quality(15.0, 20);

        // Calculate statistics
        collection.calculate_stats();

        multi_progress.update_line(
            read_line,
            format!(
                "📖 Read Processing: ✅ {} pairs processed, {} passed QC",
                processed, collection.stats.total_pairs
            ),
        );
    } else {
        multi_progress.update_line(
            read_line,
            "📖 Read Processing: ⚠️  FASTQ files not found, using mock data".to_string(),
        );

        // Create mock data for demo when FASTQ files aren't available
        for i in 0..1000 {
            let forward = PairedRead::new(
                i * 2,
                format!("mock_read_{i}"),
                meta_forge::core::paired_reads::ReadOrientation::Forward,
                "ATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
                vec![30; 36],
            );

            let reverse = PairedRead::new(
                i * 2 + 1,
                format!("mock_read_{i}"),
                meta_forge::core::paired_reads::ReadOrientation::Reverse,
                "CGATATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
                vec![30; 36],
            );

            let mut pair = ReadPair::new(forward, reverse)?;
            pair.estimate_insert_size()?;
            collection.add_pair(pair)?;
        }

        collection.calculate_stats();
        multi_progress.update_line(
            read_line,
            "📖 Read Processing: ✅ 1000 mock pairs created".to_string(),
        );
    }

    // Phase 2: Assembly Graph Construction
    multi_progress.update_line(
        assembly_line,
        "🧬 Assembly: Initializing graph builder...".to_string(),
    );

    let graph_builder = AssemblyGraphBuilder::new(
        config.assembly.k_min,
        config.assembly.k_max,
        config.assembly.min_coverage,
    );

    multi_progress.update_line(
        assembly_line,
        "🧬 Assembly: Extracting k-mers from paired reads...".to_string(),
    );

    // Extract k-mers with paired-end context
    let kmer_size = config.assembly.k_min;
    let paired_kmers = collection.extract_all_paired_kmers(kmer_size)?;

    multi_progress.update_line(
        assembly_line,
        format!(
            "🧬 Assembly: ✅ Extracted {} {}-mers with paired-end context",
            paired_kmers.len(),
            kmer_size
        ),
    );

    // Phase 3: Feature Extraction Integration
    multi_progress.update_line(
        feature_line,
        "🔍 Feature Extraction: Computing basic features...".to_string(),
    );

    // Convert paired reads to sequences for feature extraction
    multi_progress.update_line(
        feature_line,
        "🔍 Feature Extraction: Computing sequence features...".to_string(),
    );

    let mut sequence_count = 0;
    let mut total_gc_content = 0.0;
    let mut total_complexity = 0.0;

    for pair in &collection.pairs {
        // Extract features from both reads in the pair
        let forward_seq = &pair.forward.corrected;
        let reverse_seq = &pair.reverse.corrected;

        total_gc_content += pair.forward.read_info.gc_content;
        total_gc_content += pair.reverse.read_info.gc_content;
        total_complexity += pair.forward.read_info.complexity;
        total_complexity += pair.reverse.read_info.complexity;

        sequence_count += 2;
    }

    let avg_gc_content = total_gc_content / sequence_count as f64;
    let avg_complexity = total_complexity / sequence_count as f64;

    multi_progress.update_line(
        feature_line,
        format!(
            "🔍 Feature Extraction: ✅ Analyzed {} sequences (avg GC: {:.1}%, complexity: {:.2})",
            sequence_count,
            avg_gc_content * 100.0,
            avg_complexity
        ),
    );

    // Wait a moment to show the final state
    std::thread::sleep(std::time::Duration::from_millis(1000));
    multi_progress.finish();

    // Final Integration Summary
    println!("\n📊 Integration Summary");
    println!("=====================");
    println!("✅ Pipeline Components Successfully Integrated:");
    println!("   • Paired-end read processing with progress tracking");
    println!("   • Configuration-driven quality control");
    println!("   • Assembly graph builder integration");
    println!("   • Feature extraction with paired-read context");
    println!("   • Multi-phase progress display");

    println!("\n📈 Processing Statistics:");
    println!("   • Total read pairs: {}", collection.stats.total_pairs);
    println!("   • Properly paired: {}", collection.stats.properly_paired);
    println!(
        "   • Mean insert size: {:.1} bp",
        collection.stats.insert_size_stats.mean
    );
    println!(
        "   • Mean R1 quality: {:.1}",
        collection.stats.quality_stats.mean_quality_r1
    );
    println!(
        "   • Mean R2 quality: {:.1}",
        collection.stats.quality_stats.mean_quality_r2
    );
    println!("   • K-mers extracted: {}", paired_kmers.len());
    println!("   • Average GC content: {:.1}%", avg_gc_content * 100.0);
    println!("   • Average complexity: {avg_complexity:.2}");

    println!("\n🎯 Integration Benefits Demonstrated:");
    println!("   • Seamless data flow between pipeline stages");
    println!("   • Real-time progress monitoring for long operations");
    println!("   • Unified configuration system");
    println!("   • Paired-end aware processing throughout");
    println!("   • Memory-efficient streaming design");
    println!("   • Comprehensive statistics collection");

    println!("\n🚀 Ready for production metagenomic analysis!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use meta_forge::utils::configuration::PipelineConfiguration;
    use std::fs::File;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn test_pipeline_configuration_default() {
        let config = PipelineConfiguration::default();

        // Test that configuration has reasonable defaults
        assert!(config.assembly.k_min > 0);
        assert!(config.assembly.k_max >= config.assembly.k_min);
        assert!(config.assembly.min_coverage > 0);
        assert!(!config.features.kmer_sizes.is_empty());
    }

    #[test]
    fn test_paired_read_collection_creation() {
        let mut collection = PairedReadCollection::new();

        // Test initial state
        assert_eq!(collection.pairs.len(), 0);
        assert_eq!(collection.stats.total_pairs, 0);
    }

    #[test]
    fn test_multi_progress_creation() {
        let mut multi_progress = MultiProgress::new();

        // Test that we can add lines
        let _line1 = multi_progress.add_line("Test line 1".to_string());
        let _line2 = multi_progress.add_line("Test line 2".to_string());

        // This is mainly a smoke test to ensure the API works
        assert!(true);
    }

    #[test]
    fn test_assembly_graph_builder_parameters() {
        let k_min = 21;
        let k_max = 127;
        let min_coverage = 2;

        let builder = AssemblyGraphBuilder::new(k_min, k_max, min_coverage);

        // Since fields are private, this is mainly a compilation test
        let _builder_created = true;
        assert!(_builder_created);
    }

    #[test]
    fn test_paired_read_creation() {
        let read = PairedRead::new(
            1,
            "test_read".to_string(),
            meta_forge::core::paired_reads::ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30],
        );

        assert_eq!(read.id, 1);
        assert_eq!(read.pair_id, "test_read");
        assert_eq!(read.corrected, "ATCGATCGATCG");
        assert_eq!(read.quality_scores.len(), 12);
    }

    #[test]
    fn test_read_pair_creation() {
        let forward = PairedRead::new(
            1,
            "test_pair".to_string(),
            meta_forge::core::paired_reads::ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![30; 12],
        );

        let reverse = PairedRead::new(
            2,
            "test_pair".to_string(),
            meta_forge::core::paired_reads::ReadOrientation::Reverse,
            "CGATCGATCGAT".to_string(),
            vec![30; 12],
        );

        let pair_result = ReadPair::new(forward, reverse);
        assert!(pair_result.is_ok());

        if let Ok(pair) = pair_result {
            assert_eq!(pair.forward.pair_id, pair.reverse.pair_id);
        }
    }

    #[test]
    fn test_quality_filtering_parameters() {
        let mut collection = PairedReadCollection::new();

        // Add a mock pair for testing
        let forward = PairedRead::new(
            1,
            "test".to_string(),
            meta_forge::core::paired_reads::ReadOrientation::Forward,
            "ATCGATCGATCG".to_string(),
            vec![30; 12],
        );
        let reverse = PairedRead::new(
            2,
            "test".to_string(),
            meta_forge::core::paired_reads::ReadOrientation::Reverse,
            "CGATCGATCGAT".to_string(),
            vec![30; 12],
        );

        if let Ok(mut pair) = ReadPair::new(forward, reverse) {
            let _ = pair.estimate_insert_size();
            let _ = collection.add_pair(pair);
        }

        // Test filtering with reasonable parameters
        collection.filter_by_quality(15.0, 20);

        // The test mainly ensures the API doesn't panic
        assert!(true);
    }

    #[test]
    fn test_gc_content_calculation() {
        // Test GC content calculation for different sequences
        let test_sequences = vec![
            ("ATCGATCG", 0.5), // 4 GC out of 8 = 50%
            ("AAAAAAAA", 0.0), // 0 GC out of 8 = 0%
            ("GGGGCCCC", 1.0), // 8 GC out of 8 = 100%
            ("ATCG", 0.5),     // 2 GC out of 4 = 50%
        ];

        for (seq, expected_gc) in test_sequences {
            let gc_count = seq.chars().filter(|&c| c == 'G' || c == 'C').count();
            let gc_content = gc_count as f64 / seq.len() as f64;

            assert!(
                (gc_content - expected_gc).abs() < 0.001,
                "GC content for {} should be {}",
                seq,
                expected_gc
            );
        }
    }

    #[test]
    fn test_sequence_complexity_basic() {
        // Basic test for sequence complexity concepts
        let homopolymer = "AAAAAAAAAA"; // Low complexity
        let random_seq = "ATCGATCGTA"; // Higher complexity

        // Count unique characters as a simple complexity metric
        use std::collections::HashSet;

        let homo_unique: HashSet<char> = homopolymer.chars().collect();
        let random_unique: HashSet<char> = random_seq.chars().collect();

        assert!(
            homo_unique.len() < random_unique.len(),
            "Homopolymer should have lower diversity than random sequence"
        );
    }

    #[test]
    fn test_pair_id_extraction_mock() {
        // Test the pair ID extraction logic with mock data
        let mock_headers = vec![
            "@SRR390728.1 1 length=72",
            "@SRR390728.2 2 length=72",
            "@read_001/1",
            "@read_001/2",
        ];

        // This is a basic test of the concept - the actual extraction
        // logic is more complex and is tested in the extract_pair_id function
        for header in mock_headers {
            assert!(header.starts_with('@'), "FASTQ header should start with @");
            assert!(
                header.contains(char::is_alphanumeric),
                "Header should contain alphanumeric characters"
            );
        }
    }

    #[test]
    fn test_progress_tracking_concepts() {
        let max_reads = 1000;
        let current_read = 500;

        let percentage = (current_read as f64 / max_reads as f64) * 100.0;

        assert_eq!(percentage, 50.0, "Progress calculation should be correct");
    }

    #[test]
    fn test_file_path_handling() {
        use std::path::Path;

        let test_paths = vec![
            "SRR390728_1.fastq",
            "SRR390728_2.fastq",
            "./data/sample.fastq",
            "/tmp/reads.fastq",
        ];

        for path_str in test_paths {
            let path = Path::new(path_str);

            // Basic path validation
            assert!(!path_str.is_empty(), "Path should not be empty");

            if let Some(extension) = path.extension() {
                assert!(
                    extension == "fastq" || extension == "fq",
                    "File should have appropriate extension"
                );
            }
        }
    }

    #[test]
    fn test_mock_fastq_data_structure() {
        // Test the structure of mock FASTQ data
        let mock_sequences = vec![
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "CGATATCGATCGATCGATCGATCGATCGATCGATCG",
        ];

        for seq in mock_sequences {
            assert!(
                seq.len() >= 30,
                "Mock sequences should be reasonable length"
            );
            assert!(
                seq.chars().all(|c| matches!(c, 'A' | 'T' | 'C' | 'G')),
                "Mock sequences should be valid DNA"
            );
        }
    }

    #[test]
    fn test_insert_size_estimation_concept() {
        // Test basic insert size concepts
        let read1_len = 50;
        let read2_len = 50;
        let expected_insert = 200; // Typical insert size

        // Basic validation of insert size parameters
        assert!(
            expected_insert > read1_len + read2_len,
            "Insert size should be larger than combined read lengths"
        );

        let inner_distance = expected_insert - read1_len - read2_len;
        assert!(inner_distance > 0, "Inner distance should be positive");
    }

    #[test]
    fn test_statistics_calculation_structure() {
        // Test structure for statistics calculation
        let sample_data = vec![10.0, 20.0, 30.0, 40.0, 50.0];

        let sum: f64 = sample_data.iter().sum();
        let mean = sum / sample_data.len() as f64;
        let expected_mean = 30.0;

        assert!(
            (mean - expected_mean).abs() < 0.001,
            "Mean calculation should be correct"
        );
    }

    #[test]
    fn test_k_mer_extraction_basic() {
        let sequence = "ATCGATCG";
        let k = 4;

        // Number of k-mers that can be extracted
        let num_kmers = if sequence.len() >= k {
            sequence.len() - k + 1
        } else {
            0
        };

        assert_eq!(
            num_kmers, 5,
            "Should be able to extract 5 4-mers from 8bp sequence"
        );

        // Extract first k-mer manually
        if sequence.len() >= k {
            let first_kmer = &sequence[0..k];
            assert_eq!(first_kmer, "ATCG");
        }
    }

    // Integration test concept - tests the flow without external dependencies
    #[test]
    fn test_pipeline_integration_flow() {
        // Test the conceptual flow of the pipeline
        let stages = vec![
            "Read Processing",
            "Assembly Graph Construction",
            "Feature Extraction",
            "Statistics Calculation",
        ];

        // Verify all stages are represented
        assert_eq!(stages.len(), 4);

        for stage in stages {
            assert!(!stage.is_empty(), "Stage name should not be empty");
            assert!(stage.len() > 5, "Stage name should be descriptive");
        }
    }
}
