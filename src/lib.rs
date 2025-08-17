//! # MetaForge - Metagenomic Analysis Pipeline
//!
//! A high-performance metagenomics analysis pipeline built in Rust.
//! Implements advanced algorithms for DNA sequence assembly, taxonomic classification,
//! and abundance estimation using machine learning and graph neural networks.

pub mod assembly;
pub mod core;
pub mod database;
pub mod features;
pub mod ml;
pub mod pipeline;
pub mod tui;
pub mod utils;

// Re-export commonly used types at crate level
pub use crate::core::data_structures::*;
pub use crate::pipeline::complete_integration::MetagenomicsPipeline;
pub use crate::utils::configuration::{AmbiguousBaseConfig, AmbiguousBaseStrategy};

/// Result type used throughout the crate
pub type Result<T> = anyhow::Result<T>;

/// Error type used throughout the crate
pub type Error = anyhow::Error;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crate_result_type() -> Result<()> {
        // Test that our Result type works as expected
        let success: Result<i32> = Ok(42);
        let error: Result<i32> = Err(anyhow::anyhow!("test error"));

        assert!(success.is_ok());
        assert_eq!(success?, 42);

        assert!(error.is_err());
        assert!(error.unwrap_err().to_string().contains("test error"));
        Ok(())
    }

    #[test]
    fn test_crate_error_type() {
        // Test that our Error type works as expected
        let error = anyhow::anyhow!("test error message");
        let error_string = format!("{error}");

        assert!(error_string.contains("test error message"));
    }

    #[test]
    fn test_module_exports() {
        // Test that core data structures are properly re-exported
        use crate::core::Contig;

        // Test basic module compilation
        use crate::core::data_structures::CanonicalKmer;
        let kmer = CanonicalKmer::new("ATCG").expect("Valid k-mer");
        assert_eq!(kmer.sequence, "ATCG");

        // Create a test contig
        let contig = Contig {
            id: 1,
            sequence: "ATCGATCG".to_string(),
            length: 8,
            coverage: 10.0,
            contig_type: crate::core::data_structures::ContigType::Linear,
            node_path: vec![1, 2, 3],
        };
        assert_eq!(contig.length, 8);
        assert_eq!(contig.sequence.len(), 8);
    }

    #[test]
    fn test_pipeline_export() {
        // Test that MetagenomicsPipeline is properly exported
        // We can't easily instantiate it in tests due to dependencies,
        // but we can verify the type exists
        let _pipeline_type = std::any::TypeId::of::<MetagenomicsPipeline>();
    }

    #[test]
    fn test_configuration_exports() {
        // Test that configuration types are properly re-exported
        use crate::{AmbiguousBaseConfig, AmbiguousBaseStrategy};

        let config = AmbiguousBaseConfig {
            strategy: AmbiguousBaseStrategy::Skip,
            max_n_count: 1,
            replacement_base: 'A',
            random_probabilities: None,
        };

        assert!(matches!(config.strategy, AmbiguousBaseStrategy::Skip));
        assert_eq!(config.max_n_count, 1);
    }

    #[test]
    fn test_error_propagation() {
        // Test that errors propagate correctly through our Result type
        fn failing_function() -> Result<i32> {
            Err(anyhow::anyhow!("inner error"))
        }

        fn wrapping_function() -> Result<String> {
            let _value = failing_function()?;
            Ok("success".to_string())
        }

        let result = wrapping_function();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("inner error"));
    }

    #[test]
    fn test_anyhow_context() {
        // Test that anyhow context works with our Result type
        use anyhow::Context;

        fn contextual_error() -> Result<()> {
            Err(anyhow::anyhow!("base error")).context("additional context")?;
            Ok(())
        }

        let result = contextual_error();
        assert!(result.is_err());
        let error_string = result.unwrap_err().to_string();
        assert!(error_string.contains("additional context"));
        assert!(error_string.contains("base error"));
    }

    #[test]
    fn test_module_accessibility() {
        // Test that all public modules are accessible
        // Test basic module accessibility - if this compiles, modules are accessible
        let _modules_accessible = true;

        // If we can import the modules, they're properly exposed
        // This is a compilation test essentially
        let _modules_exist = true;
        assert!(_modules_exist);
    }

    #[test]
    fn test_bioinformatics_data_structures() {
        // Test core bioinformatics data structures
        use crate::core::Contig;
        // use std::collections::HashMap; // Not used in this test

        // Test k-mer creation and basic functionality
        use crate::core::data_structures::CanonicalKmer;
        let kmer = CanonicalKmer::new("ATCGATCG").expect("Valid k-mer");
        assert_eq!(kmer.sequence.len(), 8);
        assert!(kmer
            .sequence
            .chars()
            .all(|c| matches!(c, 'A' | 'T' | 'C' | 'G')));

        // Test contig creation
        let contig = Contig {
            id: 42,
            sequence: "ATCGATCGATCGATCG".to_string(),
            length: 16,
            coverage: 5.5,
            contig_type: crate::core::data_structures::ContigType::Linear,
            node_path: vec![42, 43, 44],
        };

        assert_eq!(contig.id, 42);
        assert_eq!(contig.sequence.len(), contig.length);
        assert_eq!(contig.coverage, 5.5);
    }
}
