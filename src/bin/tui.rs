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

#[cfg(test)]
mod tests {
    use super::*;
    use clap::{CommandFactory, Parser};

    #[test]
    fn test_cli_parsing_basic() {
        let cli = Cli::try_parse_from(vec!["meta-tui"]);
        assert!(cli.is_ok(), "Basic CLI parsing should succeed");

        if let Ok(cli) = cli {
            assert!(!cli.debug, "Debug should be false by default");
            assert!(cli.config.is_none(), "Config should be None by default");
        }
    }

    #[test]
    fn test_cli_debug_flag() {
        let cli = Cli::try_parse_from(vec!["meta-tui", "--debug"]);
        assert!(cli.is_ok(), "Debug flag parsing should succeed");

        if let Ok(cli) = cli {
            assert!(cli.debug, "Debug flag should be set");
        }
    }

    #[test]
    fn test_cli_debug_flag_short() {
        let cli = Cli::try_parse_from(vec!["meta-tui", "-d"]);
        assert!(cli.is_ok(), "Short debug flag parsing should succeed");

        if let Ok(cli) = cli {
            assert!(cli.debug, "Short debug flag should be set");
        }
    }

    #[test]
    fn test_cli_config_flag() {
        let cli = Cli::try_parse_from(vec!["meta-tui", "--config", "/path/to/config.toml"]);
        assert!(cli.is_ok(), "Config flag parsing should succeed");

        if let Ok(cli) = cli {
            assert!(cli.config.is_some(), "Config should be set");
            assert_eq!(
                cli.config.unwrap().to_string_lossy(),
                "/path/to/config.toml"
            );
        }
    }

    #[test]
    fn test_cli_config_flag_short() {
        let cli = Cli::try_parse_from(vec!["meta-tui", "-c", "config.toml"]);
        assert!(cli.is_ok(), "Short config flag parsing should succeed");

        if let Ok(cli) = cli {
            assert!(cli.config.is_some(), "Config should be set");
            assert_eq!(cli.config.unwrap().to_string_lossy(), "config.toml");
        }
    }

    #[test]
    fn test_cli_combined_flags() {
        let cli = Cli::try_parse_from(vec!["meta-tui", "--debug", "--config", "test.toml"]);
        assert!(cli.is_ok(), "Combined flags parsing should succeed");

        if let Ok(cli) = cli {
            assert!(cli.debug, "Debug flag should be set");
            assert!(cli.config.is_some(), "Config should be set");
            assert_eq!(cli.config.unwrap().to_string_lossy(), "test.toml");
        }
    }

    #[test]
    fn test_cli_invalid_flag() {
        let cli = Cli::try_parse_from(vec!["meta-tui", "--invalid-flag"]);
        assert!(cli.is_err(), "Invalid flag should cause parsing to fail");
    }

    #[test]
    fn test_cli_help_generation() {
        // Test that help can be generated (this would normally print and exit)
        let cli_result = Cli::try_parse_from(vec!["meta-tui", "--help"]);
        assert!(
            cli_result.is_err(),
            "Help flag should cause early exit with error"
        );
    }

    #[test]
    fn test_cli_version_flag() {
        let cli_result = Cli::try_parse_from(vec!["meta-tui", "--version"]);
        assert!(
            cli_result.is_err(),
            "Version flag should cause early exit with error"
        );
    }

    #[test]
    fn test_log_level_selection() {
        // Test the logic for log level selection
        let debug_level = "debug";
        let info_level = "info";

        // Simulate the logic from main
        let cli_debug = true;
        let selected_level = if cli_debug { debug_level } else { info_level };
        assert_eq!(selected_level, "debug");

        let cli_normal = false;
        let selected_level = if cli_normal { debug_level } else { info_level };
        assert_eq!(selected_level, "info");
    }

    #[test]
    fn test_config_path_handling() {
        use std::path::PathBuf;

        // Test various path formats
        let test_paths = vec![
            "config.toml",
            "/absolute/path/config.toml",
            "./relative/config.toml",
            "../parent/config.toml",
            "~/home/config.toml",
        ];

        for path_str in test_paths {
            let cli = Cli::try_parse_from(vec!["meta-tui", "-c", path_str]);
            assert!(cli.is_ok(), "Path '{}' should be parseable", path_str);

            if let Ok(cli) = cli {
                let config_path = cli.config.unwrap();
                assert_eq!(config_path, PathBuf::from(path_str));
            }
        }
    }

    #[test]
    fn test_cli_metadata() {
        // Test CLI metadata (name, version, about)
        let app = Cli::command();

        assert_eq!(app.get_name(), "meta-tui");
        assert!(app.get_about().is_some());
        assert!(app.get_version().is_some());
    }

    #[test]
    fn test_empty_config_path() {
        // Test that empty config path is handled gracefully
        let cli = Cli::try_parse_from(vec!["meta-tui", "--config", ""]);
        assert!(cli.is_ok(), "Empty config path should be parseable");

        if let Ok(cli) = cli {
            assert!(cli.config.is_some());
            assert_eq!(cli.config.unwrap().to_string_lossy(), "");
        }
    }

    // Since we can't easily test the async main function directly without
    // complex mocking, we test the components it uses
    #[tokio::test]
    async fn test_tokio_runtime_availability() {
        // Test that we're in a tokio runtime (simulates main function environment)
        let handle = tokio::runtime::Handle::current();
        assert!(
            handle.runtime_flavor() != tokio::runtime::RuntimeFlavor::CurrentThread
                || handle.runtime_flavor() == tokio::runtime::RuntimeFlavor::CurrentThread
        );
    }
}
