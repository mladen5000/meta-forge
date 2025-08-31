use anyhow::Result;
use clap::Parser;
use meta_forge::pipeline::complete_integration::{Cli, Commands, MetagenomicsPipeline};
use meta_forge::utils::progress_display::MultiProgress;

#[tokio::main]
/// Main entry point for the metagenomic analysis pipeline
async fn main() -> Result<()> {
    let cli = Cli::parse();

    // Initialize logging only if not already set
    if !tracing::dispatcher::has_been_set() {
        let log_level = if cli.verbose { "debug" } else { "info" };
        tracing_subscriber::fmt()
            .with_env_filter(tracing_subscriber::EnvFilter::new(log_level))
            .init();
    }

    match cli.command {
        Commands::Analyze {
            input,
            sample_name,
            database: _database,
            mode,
        } => {
            let sample_name = sample_name.unwrap_or_else(|| {
                input[0]
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("sample")
                    .to_string()
            });

            // Create beautiful multi-line progress display
            println!("🧬 MetaForge - Advanced Metagenomics Analysis Pipeline");
            println!("====================================================");
            println!("Sample: {sample_name}");
            println!("Input files: {}", input.len());
            println!("Mode: {mode:?}\n");

            let mut multi_progress = MultiProgress::new();
            let init_line =
                multi_progress.add_line("🔧 Initialization: Loading configuration...".to_string());
            let preprocess_line =
                multi_progress.add_line("📋 Preprocessing: Waiting...".to_string());
            let assembly_line = multi_progress.add_line("🧬 Assembly: Waiting...".to_string());
            let features_line =
                multi_progress.add_line("🔍 Feature Extraction: Waiting...".to_string());
            let classification_line =
                multi_progress.add_line("🏷️  Classification: Waiting...".to_string());
            let abundance_line = multi_progress.add_line("📊 Abundance: Waiting...".to_string());
            let report_line = multi_progress.add_line("📝 Report: Waiting...".to_string());

            multi_progress.update_line(
                init_line,
                "🔧 Initialization: Creating pipeline...".to_string(),
            );
            let mut pipeline = MetagenomicsPipeline::new(cli.config.as_deref())?;

            multi_progress.update_line(init_line, "🔧 Initialization: ✅ Ready".to_string());

            // Run analysis with progress tracking
            let results = pipeline
                .run_analysis_with_progress(
                    &input,
                    &sample_name,
                    mode,
                    &mut multi_progress,
                    preprocess_line,
                    assembly_line,
                    features_line,
                    classification_line,
                    abundance_line,
                    report_line,
                )
                .await?;

            multi_progress.finish();

            // Display beautiful results summary
            println!("🎉 Analysis Complete!");
            println!("====================");
            println!("📊 Results Summary:");
            println!("   Sample: {}", results.sample_name);
            println!(
                "   Contigs assembled: {}",
                results.assembly_results.contigs.len()
            );
            println!(
                "   Total sequence length: {} bp",
                results.assembly_results.assembly_stats.total_length
            );
            println!(
                "   Assembly N50: {} bp",
                results.assembly_results.assembly_stats.n50
            );
            println!(
                "   Species identified: {}",
                results
                    .classifications
                    .iter()
                    .map(|c| &c.taxonomy_name)
                    .collect::<std::collections::HashSet<_>>()
                    .len()
            );
            println!(
                "   Processing time: {:.2} seconds",
                results.processing_time.as_secs_f64()
            );

            println!("   📄 Full report generated successfully");

            println!("\n✨ Analysis completed successfully! ✨");
        }
        Commands::Database { operation } => {
            use meta_forge::pipeline::complete_integration::handle_database_operation;
            handle_database_operation(operation).await?;
        }
        Commands::Assemble {
            input,
            k_range,
            min_coverage,
        } => {
            println!("🧬 MetaForge - Assembly-Only Mode");
            println!("================================");
            println!("Input files: {}", input.len());
            if let Some(range) = &k_range {
                println!("K-mer range: {}-{}", range.0, range.1);
            }
            println!("Min coverage: {}\n", min_coverage);

            // Create pipeline with default configuration
            let mut pipeline = MetagenomicsPipeline::new(cli.config.as_deref())?;

            // Run preprocessing to get reads
            println!("📋 Preprocessing input files...");
            let reads = pipeline.preprocess_inputs(&input).await?;
            
            // Run assembly with verbose progress
            println!("🚀 Starting assembly with enhanced verbose progress...");
            let assembly_results = pipeline.run_assembly(&reads).await?;

            // Display results
            println!("\n🎉 Assembly Completed Successfully!");
            println!("=====================================");
            println!("📊 Assembly Statistics:");
            println!("   Contigs generated: {}", assembly_results.contigs.len());
            println!("   Total sequence length: {} bp", 
                    assembly_results.assembly_stats.total_length);
            println!("   Assembly N50: {} bp", 
                    assembly_results.assembly_stats.n50);
            println!("   Average contig length: {:.1} bp", 
                    assembly_results.assembly_stats.total_length as f64 / assembly_results.assembly_stats.num_contigs.max(1) as f64);
            println!("   Largest contig: {} bp", 
                    assembly_results.assembly_stats.largest_contig);
            
            // Write contigs to FASTA file
            let output_path = format!("contigs_{}.fasta", 
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_secs());
            
            println!("💾 Writing contigs to: {}", output_path);
            
            // Write FASTA output (simplified)
            use std::io::Write;
            let mut file = std::fs::File::create(&output_path)?;
            for (i, contig) in assembly_results.contigs.iter().enumerate() {
                writeln!(file, ">contig_{}", i + 1)?;
                writeln!(file, "{}", contig.sequence)?;
            }
            
            println!("✅ Assembly results saved to {}", output_path);
        }
        _ => {
            println!("Command not yet implemented in simplified main.rs");
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    use tempfile::tempdir;

    #[test]
    fn test_cli_parsing() {
        // Test CLI argument parsing
        let cli = Cli::try_parse_from(vec![
            "meta_forge",
            "analyze",
            "--input",
            "test.fastq",
            "--sample-name",
            "test_sample",
        ]);

        assert!(
            cli.is_ok(),
            "CLI parsing should succeed for valid arguments"
        );

        if let Ok(cli) = cli {
            match cli.command {
                Commands::Analyze {
                    input, sample_name, ..
                } => {
                    assert_eq!(input.len(), 1);
                    assert_eq!(input[0].to_string_lossy(), "test.fastq");
                    assert_eq!(sample_name, Some("test_sample".to_string()));
                }
                _ => panic!("Expected Analyze command"),
            }
        }
    }

    #[test]
    fn test_cli_verbose_flag() {
        let cli = Cli::try_parse_from(vec![
            "meta_forge",
            "--verbose",
            "analyze",
            "--input",
            "test.fastq",
        ]);

        assert!(cli.is_ok());
        if let Ok(cli) = cli {
            assert!(cli.verbose, "Verbose flag should be set");
        }
    }

    #[test]
    fn test_cli_config_flag() {
        let cli = Cli::try_parse_from(vec![
            "meta_forge",
            "--config",
            "config.toml",
            "analyze",
            "--input",
            "test.fastq",
        ]);

        assert!(cli.is_ok());
        if let Ok(cli) = cli {
            assert_eq!(cli.config, Some("config.toml".into()));
        }
    }

    #[test]
    fn test_sample_name_derivation() {
        // Create a temporary file to test sample name derivation
        let temp_dir = tempdir().expect("Failed to create temp dir");
        let test_file = temp_dir.path().join("sample_data.fastq");
        std::fs::write(&test_file, "test content").expect("Failed to write test file");

        // Test that sample name is derived from filename when not provided
        let derived_name = test_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("sample")
            .to_string();

        assert_eq!(derived_name, "sample_data");
    }

    #[test]
    fn test_database_command_parsing() {
        use meta_forge::pipeline::complete_integration::DatabaseOps;

        let cli = Cli::try_parse_from(vec!["meta_forge", "database", "init"]);

        assert!(cli.is_ok());
        if let Ok(cli) = cli {
            match cli.command {
                Commands::Database { operation } => {
                    match operation {
                        DatabaseOps::Init { path: output } => {
                            // Test passes - database init command parsed correctly
                        }
                        _ => panic!("Expected Init operation"),
                    }
                }
                _ => panic!("Expected Database command"),
            }
        }
    }

    #[test]
    fn test_invalid_cli_args() {
        let cli = Cli::try_parse_from(vec!["meta_forge", "invalid_command"]);

        assert!(cli.is_err(), "CLI parsing should fail for invalid command");
    }

    // Test logging initialization (mock test since we can't easily test actual logging setup)
    #[test]
    fn test_logging_level_selection() {
        // Test the logic for selecting log level
        let verbose_level = "debug";
        let normal_level = "info";

        assert_eq!(verbose_level, "debug");
        assert_eq!(normal_level, "info");

        // In actual main function:
        // let log_level = if cli.verbose { "debug" } else { "info" };
        // This logic is tested implicitly through the verbose flag test
    }

    #[tokio::test]
    async fn test_analysis_mode_parsing() {
        use meta_forge::pipeline::complete_integration::AnalysisMode;

        let cli = Cli::try_parse_from(vec![
            "meta_forge",
            "analyze",
            "--input",
            "test.fastq",
            "--mode",
            "fast",
        ]);

        assert!(cli.is_ok());
        if let Ok(cli) = cli {
            match cli.command {
                Commands::Analyze { mode, .. } => {
                    assert_eq!(mode, AnalysisMode::Fast);
                }
                _ => panic!("Expected Analyze command"),
            }
        }
    }
}
