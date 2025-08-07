use anyhow::Result;
use clap::Parser;

/// Meta-Forge TUI - Interactive Terminal Interface
#[derive(Parser)]
#[command(name = "meta-tui")]
#[command(about = "Interactive terminal interface for metagenomic analysis")]
#[command(version = "0.1.0")]
struct Cli {
    /// Enable debug logging
    #[arg(short, long)]
    debug: bool,
    
    /// Configuration file path
    #[arg(short, long)]
    config: Option<std::path::PathBuf>,
}

#[tokio::main]
async fn main() -> Result<()> {
    let cli = Cli::parse();
    
    // Initialize logging
    let log_level = if cli.debug { "debug" } else { "info" };
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::new(log_level))
        .init();
    
    // Initialize and run TUI
    meta_forge::tui::app::run_tui().await?;
    
    Ok(())
}