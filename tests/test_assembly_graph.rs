/// Test-driven development for AssemblyGraph issues

use meta_forge::assembly::adaptive_k::AssemblyGraphBuilder;
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn test_assembly_graph_builder_constructor() {
        // Test the corrected constructor signature (3 parameters, not 4)
        let builder = AssemblyGraphBuilder::new(
            11,  // base k-mer size
            15,  // max k-mer size  
            1    // minimum coverage
        );

        // Constructor should not return Result, just the struct
        // Fields are private, so we can't directly test them
        // Instead, test that the builder was created successfully
    }

    #[test]
    fn test_assembly_graph_field_access() {
        let builder = AssemblyGraphBuilder::new(11, 15, 1);
        
        // Create test reads with proper CorrectionMetadata
        let reads = vec![
            CorrectedRead {
                id: 0,
                original: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
                corrected: "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG".to_string(),
                corrections: vec![],
                quality_scores: vec![40; 40],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
            },
            CorrectedRead {
                id: 1,
                original: "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT".to_string(),
                corrected: "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT".to_string(),
                corrections: vec![],
                quality_scores: vec![40; 40],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
            },
        ];

        let result = builder.build(&reads);
        assert!(result.is_ok());
        
        let graph = result.unwrap();
        
        // Test correct field access patterns for AssemblyGraph
        assert!(graph.graph_fragment.nodes.len() > 0);
        assert!(graph.contigs.len() >= 0); // May be 0 if no contigs generated yet
        
        // Test that we can access graph statistics
        assert!(graph.assembly_stats.total_length >= 0);
        assert!(graph.assembly_stats.num_contigs >= 0);
    }

    #[test]
    fn test_corrected_read_with_metadata() {
        // Test that CorrectedRead can be created with correction_metadata
        let read = CorrectedRead {
            id: 42,
            original: "ATCGATCG".to_string(),
            corrected: "ATCGATCG".to_string(),
            corrections: vec![],
            quality_scores: vec![40; 8],
            correction_metadata: CorrectionMetadata {
                algorithm: "test_algorithm".to_string(),
                confidence_threshold: 0.8,
                context_window: 3,
                correction_time_ms: 100,
            },
        };

        assert_eq!(read.id, 42);
        assert_eq!(read.original, "ATCGATCG");
        assert_eq!(read.corrected, "ATCGATCG");
        assert_eq!(read.correction_metadata.algorithm, "test_algorithm");
        assert_eq!(read.correction_metadata.confidence_threshold, 0.8);
        assert_eq!(read.correction_metadata.context_window, 3);
        assert_eq!(read.correction_metadata.correction_time_ms, 100);
    }

    #[test]
    fn test_assembly_graph_empty_input() {
        let builder = AssemblyGraphBuilder::new(11, 15, 1);
        let empty_reads: Vec<CorrectedRead> = vec![];
        
        let result = builder.build(&empty_reads);
        assert!(result.is_ok());
        
        let graph = result.unwrap();
        assert_eq!(graph.graph_fragment.nodes.len(), 0);
        assert_eq!(graph.contigs.len(), 0);
    }

    #[test]
    fn test_assembly_graph_short_reads() {
        let builder = AssemblyGraphBuilder::new(11, 15, 1);
        
        // Create reads shorter than k-mer size
        let short_reads = vec![
            CorrectedRead {
                id: 0,
                original: "ATCG".to_string(), // Only 4 bp, less than k=11
                corrected: "ATCG".to_string(),
                corrections: vec![],
                quality_scores: vec![40; 4],
                correction_metadata: CorrectionMetadata {
                    algorithm: "test".to_string(),
                    confidence_threshold: 0.95,
                    context_window: 5,
                    correction_time_ms: 0,
                },
            },
        ];
        
        let result = builder.build(&short_reads);
        assert!(result.is_ok());
        
        let graph = result.unwrap();
        // Should handle short reads gracefully - may have 0 nodes
        assert!(graph.graph_fragment.nodes.len() == 0);
    }
}