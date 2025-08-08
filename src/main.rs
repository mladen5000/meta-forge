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
            let init_line = multi_progress.add_line("🔧 Initialization: Loading configuration...".to_string());
            let preprocess_line = multi_progress.add_line("📋 Preprocessing: Waiting...".to_string());
            let assembly_line = multi_progress.add_line("🧬 Assembly: Waiting...".to_string());
            let features_line = multi_progress.add_line("🔍 Feature Extraction: Waiting...".to_string());
            let classification_line = multi_progress.add_line("🏷️  Classification: Waiting...".to_string());
            let abundance_line = multi_progress.add_line("📊 Abundance: Waiting...".to_string());
            let report_line = multi_progress.add_line("📝 Report: Waiting...".to_string());

            multi_progress.update_line(init_line, "🔧 Initialization: Creating pipeline...".to_string());
            let mut pipeline = MetagenomicsPipeline::new(cli.config.as_deref())?;

            multi_progress.update_line(init_line, "🔧 Initialization: ✅ Ready".to_string());
            
            // Run analysis with progress tracking
            let results = pipeline.run_analysis_with_progress(&input, &sample_name, mode, &mut multi_progress, 
                                                              preprocess_line, assembly_line, features_line, 
                                                              classification_line, abundance_line, report_line).await?;

            multi_progress.finish();

            // Display beautiful results summary
            println!("🎉 Analysis Complete!");
            println!("====================");
            println!("📊 Results Summary:");
            println!("   Sample: {}", results.sample_name);
            println!("   Contigs assembled: {}", results.assembly_results.contigs.len());
            println!("   Total sequence length: {} bp", results.assembly_results.assembly_stats.total_length);
            println!("   Assembly N50: {} bp", results.assembly_results.assembly_stats.n50);
            println!("   Species identified: {}", results.classifications.len());
            println!("   Processing time: {:.2} seconds", results.processing_time.as_secs_f64());
            
            println!("   📄 Full report generated successfully");
            
            println!("\n✨ Analysis completed successfully! ✨");
        }
        Commands::Database { operation } => {
            use meta_forge::pipeline::complete_integration::handle_database_operation;
            handle_database_operation(operation).await?;
        }
        _ => {
            println!("Command not yet implemented in simplified main.rs");
        }
    }

    Ok(())
}
