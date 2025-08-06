use anyhow::Result;
use clap::Parser;
use meta_forge::pipeline::complete_integration::{Cli, Commands, MetagenomicsPipeline};

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();

    // Initialize logging
    let log_level = if cli.verbose { "debug" } else { "info" };
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::new(log_level))
        .init();

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

            let mut pipeline = MetagenomicsPipeline::new(cli.config.as_deref())?;

            // Note: CLI parameter overrides would need to be implemented
            // in the pipeline configuration system

            let results = pipeline.run_analysis(&input, &sample_name, mode).await?;

            println!("✅ Analysis completed successfully!");
            println!("📊 Results:");
            println!("   Sample: {}", results.sample_name);
            println!("   Contigs: {}", results.assembly_results.contigs.len());
            println!(
                "   Total length: {} bp",
                results.assembly_results.assembly_stats.total_length
            );
            println!("   N50: {} bp", results.assembly_results.assembly_stats.n50);
            println!(
                "   Processing time: {:.2} seconds",
                results.processing_time.as_secs_f64()
            );
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
