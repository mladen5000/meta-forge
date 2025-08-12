use anyhow::Result;
use meta_forge::utils::intermediate_output::{
    IntermediateOutputManager, OutputConfig, PipelineSection,
};
use serde_json::json;
use std::fs;
use tempfile::tempdir;

#[tokio::test]
async fn test_complete_intermediate_workflow() -> Result<()> {
    let temp_dir = tempdir()?;
    let config = OutputConfig {
        enable_json: true,
        enable_binary: true,
        enable_fasta: true,
        enable_tsv: true,
        compress_files: false, // Disable compression for easier testing
        max_file_size_mb: 10,
    };

    // Create output manager
    let manager = IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config)?;

    // Test 1: Save and load JSON intermediate data
    let test_data = json!({
        "reads_processed": 1000,
        "corrections_made": 25,
        "quality_score": 0.95
    });

    let metadata = manager.save_intermediate(
        PipelineSection::Preprocessing,
        "test_data",
        &test_data,
        json!({"test_type": "preprocessing"}),
    )?;

    assert_eq!(metadata.section, PipelineSection::Preprocessing);
    assert_eq!(metadata.format, "json");
    assert!(!metadata.checksums.md5.is_empty());
    assert!(!metadata.checksums.sha256.is_empty());

    // Load the data back
    let loaded_data: serde_json::Value =
        manager.load_intermediate(PipelineSection::Preprocessing, "test_data")?;

    assert_eq!(loaded_data["reads_processed"], 1000);
    assert_eq!(loaded_data["corrections_made"], 25);

    // Test 2: Save FASTA sequences
    let sequences = vec![
        ("contig_1".to_string(), "ATCGATCGATCGATCG".to_string()),
        ("contig_2".to_string(), "GCTAGCTAGCTAGCTA".to_string()),
        ("contig_3".to_string(), "TTAAGGCCTTAAGGCC".to_string()),
    ];

    let fasta_metadata = manager.save_sequences(
        PipelineSection::Assembly,
        "contigs",
        &sequences,
        json!({"num_contigs": sequences.len()}),
    )?;

    assert_eq!(fasta_metadata.section, PipelineSection::Assembly);
    assert_eq!(fasta_metadata.format, "fasta");

    // Test 3: Save TSV data
    let headers = vec![
        "contig_id".to_string(),
        "taxonomy".to_string(),
        "confidence".to_string(),
    ];
    let rows = vec![
        vec![
            "contig_1".to_string(),
            "E. coli".to_string(),
            "0.95".to_string(),
        ],
        vec![
            "contig_2".to_string(),
            "B. subtilis".to_string(),
            "0.87".to_string(),
        ],
        vec![
            "contig_3".to_string(),
            "S. aureus".to_string(),
            "0.92".to_string(),
        ],
    ];

    let tsv_metadata = manager.save_tsv(
        PipelineSection::Classification,
        "classifications",
        &headers,
        &rows,
        json!({"num_classifications": rows.len()}),
    )?;

    assert_eq!(tsv_metadata.section, PipelineSection::Classification);
    assert_eq!(tsv_metadata.format, "tsv");

    // Test 4: Check file existence
    assert!(manager.has_intermediate(PipelineSection::Preprocessing, "test_data"));
    assert!(manager.has_intermediate(PipelineSection::Assembly, "contigs"));
    assert!(manager.has_intermediate(PipelineSection::Classification, "classifications"));
    assert!(!manager.has_intermediate(PipelineSection::Features, "nonexistent"));

    // Test 5: Generate run summary
    let summary = manager.generate_run_summary()?;
    assert_eq!(summary.run_id, manager.run_id);
    assert!(!summary.sections.is_empty());

    // Should have files in preprocessing, assembly, and classification sections
    assert!(summary
        .sections
        .contains_key(&PipelineSection::Preprocessing));
    assert!(summary.sections.contains_key(&PipelineSection::Assembly));
    assert!(summary
        .sections
        .contains_key(&PipelineSection::Classification));

    // Test 6: Verify directory structure
    let run_dir = &manager.run_dir;
    assert!(run_dir.exists());
    assert!(run_dir.join("preprocessing").exists());
    assert!(run_dir.join("assembly").exists());
    assert!(run_dir.join("features").exists());
    assert!(run_dir.join("classification").exists());
    assert!(run_dir.join("abundance").exists());
    assert!(run_dir.join("report").exists());

    // Test 7: Check actual file contents
    let json_file = run_dir.join("preprocessing/test_data.json");
    assert!(json_file.exists());
    let json_content = fs::read_to_string(&json_file)?;
    let parsed: serde_json::Value = serde_json::from_str(&json_content)?;
    assert_eq!(parsed["reads_processed"], 1000);

    let fasta_file = run_dir.join("assembly/contigs.fasta");
    assert!(fasta_file.exists());
    let fasta_content = fs::read_to_string(&fasta_file)?;
    assert!(fasta_content.contains(">contig_1"));
    assert!(fasta_content.contains("ATCGATCGATCGATCG"));
    assert!(fasta_content.contains(">contig_2"));
    assert!(fasta_content.contains("GCTAGCTAGCTAGCTA"));

    let tsv_file = run_dir.join("classification/classifications.tsv");
    assert!(tsv_file.exists());
    let tsv_content = fs::read_to_string(&tsv_file)?;
    assert!(tsv_content.contains("contig_id\ttaxonomy\tconfidence"));
    assert!(tsv_content.contains("contig_1\tE. coli\t0.95"));
    assert!(tsv_content.contains("contig_2\tB. subtilis\t0.87"));

    println!("✅ All intermediate output tests passed!");
    println!("📁 Run directory: {}", run_dir.display());
    println!(
        "📊 Run summary: {} sections with files",
        summary.sections.len()
    );

    Ok(())
}

#[tokio::test]
async fn test_compressed_output() -> Result<()> {
    let temp_dir = tempdir()?;
    let config = OutputConfig {
        enable_json: true,
        enable_binary: false,
        enable_fasta: true,
        enable_tsv: false,
        compress_files: true, // Enable compression
        max_file_size_mb: 10,
    };

    let manager = IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config)?;

    // Save compressed data
    let test_data = json!({
        "large_dataset": vec![1, 2, 3, 4, 5],
        "metadata": "This is test data for compression"
    });

    manager.save_intermediate(
        PipelineSection::Features,
        "compressed_data",
        &test_data,
        json!({"compressed": true}),
    )?;

    // Check that compressed file was created
    let compressed_file = manager.run_dir.join("features/compressed_data.json.gz");
    assert!(compressed_file.exists());

    // Should NOT have uncompressed file
    let uncompressed_file = manager.run_dir.join("features/compressed_data.json");
    assert!(!uncompressed_file.exists());

    // Load compressed data
    let loaded_data: serde_json::Value =
        manager.load_intermediate(PipelineSection::Features, "compressed_data")?;

    assert_eq!(loaded_data["metadata"], "This is test data for compression");
    assert_eq!(loaded_data["large_dataset"].as_array().unwrap().len(), 5);

    println!("✅ Compressed output test passed!");
    Ok(())
}

#[tokio::test]
async fn test_multiple_runs() -> Result<()> {
    let temp_dir = tempdir()?;
    let config = OutputConfig::default();

    // Create first run
    let manager1 = IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config.clone())?;
    manager1.save_intermediate(
        PipelineSection::Assembly,
        "run1_data",
        &json!({"run": 1}),
        json!({}),
    )?;

    // Wait a moment to ensure different timestamps
    tokio::time::sleep(tokio::time::Duration::from_millis(1001)).await;

    // Create second run
    let manager2 = IntermediateOutputManager::new(temp_dir.path().to_path_buf(), config)?;
    manager2.save_intermediate(
        PipelineSection::Assembly,
        "run2_data",
        &json!({"run": 2}),
        json!({}),
    )?;

    // List run directories
    let run_dirs = manager1.list_run_directories()?;
    assert_eq!(run_dirs.len(), 2);

    // Check that both runs exist
    assert!(manager1.has_intermediate(PipelineSection::Assembly, "run1_data"));
    assert!(manager2.has_intermediate(PipelineSection::Assembly, "run2_data"));

    // Cross-check that runs are isolated
    assert!(!manager1.has_intermediate(PipelineSection::Assembly, "run2_data"));
    assert!(!manager2.has_intermediate(PipelineSection::Assembly, "run1_data"));

    println!("✅ Multiple runs test passed!");
    println!("📁 Created {} run directories", run_dirs.len());

    Ok(())
}

#[test]
fn test_pipeline_section_conversions() {
    assert_eq!(PipelineSection::Preprocessing.as_str(), "preprocessing");
    assert_eq!(PipelineSection::QualityControl.as_str(), "quality_control");
    assert_eq!(PipelineSection::Assembly.as_str(), "assembly");
    assert_eq!(PipelineSection::Features.as_str(), "features");
    assert_eq!(PipelineSection::Classification.as_str(), "classification");
    assert_eq!(PipelineSection::Abundance.as_str(), "abundance");
    assert_eq!(PipelineSection::Report.as_str(), "report");
}
