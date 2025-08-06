# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a high-performance metagenomics analysis pipeline built in Rust. The project implements advanced algorithms for DNA sequence assembly, taxonomic classification, and abundance estimation using machine learning and graph neural networks.

## Core Architecture

The project is structured around several key modules:

- **core_data_structures.rs** - Fundamental data structures including CanonicalKmer, GraphNode, GraphEdge, GraphFragment, and AssemblyChunk
- **complete_pipeline_integration.rs** - Main CLI application and pipeline orchestrator (MetagenomicsPipeline struct)
- **assembly_graph_construction.rs** - Assembly graph building and contig generation
- **adaptive_k_assembly.rs** - Adaptive k-mer selection algorithms (AdaptiveGraph)
- **gnn_repeat_resolution.rs** - Graph neural network for repeat resolution (RepeatResolverGNN)
- **feature_extraction.rs** - Advanced feature extraction from sequences and graphs
- **database_integration.rs** - SQLite database operations and k-mer indexing
- **configuration_management.rs** - Configuration system with TOML support
- **comprehensive_test_suite.rs** - Testing framework and benchmarks
- **streaming_abundance_estimator.rs** - HyperLogLog-based abundance estimation
- **learned_bloom_filter.rs** - ML-enhanced bloom filters

## Common Commands

### Build and Run
```bash
# Build the project
cargo build --release

# Run the main CLI application
cargo run --release -- --help

# Run complete analysis
cargo run --release -- analyze input.fastq --sample-name "sample1" --mode standard

# Assembly only
cargo run --release -- assemble reads.fastq --k-range 21-31 --min-coverage 3

# Feature extraction
cargo run --release -- features sequences.fasta --types composition,patterns --format json
```

### Testing
```bash
# Run all tests
cargo test

# Run specific test module
cargo test core_data_structures

# Run benchmarks
cargo test --release --features criterion

# Run integration tests
cargo test integration_tests
```

### Database Operations
```bash
# Initialize database
cargo run --release -- database init ./data/metagenomics.db

# Import taxonomy data
cargo run --release -- database import-taxonomy ./data/metagenomics.db taxonomy.txt

# Build k-mer index
cargo run --release -- database build-index ./data/metagenomics.db --k 21
```

## Configuration

The pipeline uses TOML configuration files managed by the ConfigurationManager:

```bash
# Generate configuration templates
cargo run --release -- config standard --output config.toml
cargo run --release -- config high-performance --output high_perf.toml
cargo run --release -- config low-memory --output low_mem.toml
```

Configuration sections include:
- `general` - working directories, output paths, debug mode
- `assembly` - k-mer ranges, coverage thresholds, graph simplification
- `features` - feature extraction parameters and dimensions
- `database` - database connection and indexing settings
- `ml` - machine learning model configurations
- `performance` - threading, memory limits, monitoring
- `logging` - log levels and output formats
- `io` - input/output formats and compression

## Key Data Flow

1. **Input Processing**: FASTQ/FASTA files → CorrectedRead structures
2. **Assembly**: Reads → GraphFragment → AssemblyGraph → Contigs
3. **Feature Extraction**: Sequences → FeatureVector (composition, complexity, topology)
4. **Classification**: Features → TaxonomicClassification via ML models
5. **Abundance Estimation**: K-mers → AbundanceProfile via streaming algorithms
6. **Reporting**: Combined results → AnalysisReport (JSON/HTML/TSV)

## Development Notes

- The project uses Rust 2024 edition
- Async/await for I/O operations with tokio runtime
- Heavy use of traits for modularity (e.g., AssemblyGraphBuilder)
- Comprehensive error handling with anyhow and thiserror
- Serialization with serde for data persistence
- Performance monitoring with tracing
- CLI built with clap 4.x
- Graph operations using petgraph
- Linear algebra with ndarray
- Database operations with rusqlite

## Testing Framework

The comprehensive_test_suite.rs provides:
- TestDataGenerator for creating test sequences
- Unit tests for each module
- Integration tests for full pipeline
- Benchmarks using criterion
- TestRunner for orchestrating tests

Run tests with different data sizes:
```bash
cargo test -- --ignored  # Long-running tests
cargo test test_large_dataset  # Specific large tests
```

## Performance Considerations

- Uses memory-mapped files for large inputs
- Streaming algorithms for abundance estimation
- Multi-threaded assembly and feature extraction
- Configurable memory limits and thread counts
- Resource monitoring during execution

## Machine Learning Components

- ONNX runtime for neural network inference
- Graph neural networks for repeat resolution
- Learned data structures (bloom filters)
- Feature extraction for ML classification
- Model serialization and loading

## File Formats Supported

**Input**: FASTQ, FASTA (plain or gzipped)
**Output**: JSON, HTML reports, TSV summaries, FASTA contigs
**Configuration**: TOML files
**Database**: SQLite with custom schema