# MetaForge Development Guide

## Project Overview

MetaForge is a high-performance Rust-based metagenomics pipeline featuring:
- Adaptive k-mer assembly with graph optimization
- K-mer-based taxonomic classification with ML integration
- Real-time error correction during assembly
- Comprehensive sequence and graph feature extraction
- Modern terminal UI with real-time progress tracking
- SQLite-backed k-mer and taxonomy databases
- Streaming output with multiple format support (TSV, JSON, GFA, Kraken2)

**Note**: This is a clean, maintained codebase with recent housekeeping (Oct 2025) removing dead code and unused dependencies.

## Quick Start

### Build
```bash
cargo build --release
```

### Run Analysis
```bash
# Basic analysis
./target/release/meta-forge analyze reads.fastq

# With memory and thread control
./target/release/meta-forge -m 4096 -j 8 analyze reads.fastq

# With all options
./target/release/meta-forge -m 4096 -j 8 -o ./results analyze reads.fastq --sample-name my_sample
```

## Important CLI Rules

### ⚠️ Flag Order Matters!
**Global flags MUST come BEFORE the subcommand:**

```bash
# ✅ CORRECT
meta-forge -m 4096 -j 8 analyze input.fastq

# ❌ WRONG (will cause "Invalid input" error)
meta-forge analyze input.fastq -m 4096
```

### Global Flags (before subcommand)
- `-m <MB>` - Memory budget in MB
- `-j <N>` - Number of threads
- `-o <DIR>` - Output directory
- `-c <FILE>` - Config file
- `-v` - Verbose logging

### Subcommands
- `analyze` - Full analysis pipeline
- `assemble` - Assembly only
- `features` - Feature extraction only
- `database` - Database operations
- `config` - Generate config template

## Development Practices

### Code Organization
- **Never save files to root** - use appropriate subdirectories:
  - `/src` - Source code
  - `/tests` - Test files
  - `/docs` - Documentation
  - `/benches` - Benchmarks
  - `/examples` - Example code

### Testing
```bash
# Run all tests
cargo test

# Run specific test
cargo test test_name

# Run with output
cargo test -- --nocapture
```

### Performance
```bash
# Release build (optimized)
cargo build --release

# Run benchmarks
cargo bench

# Profile with flamegraph
cargo flamegraph --bin meta-forge -- analyze input.fastq
```

## Architecture

### Core Modules
- **`assembly/`** - Adaptive k-mer assembly with optimized graph construction
  - `laptop_assembly.rs` - **Production assembler** (use this!)
  - `adaptive_k.rs` - K-mer size selection
  - `orchestrator.rs` - Assembly workflow coordination
  - `optimized/` - Experimental optimizations (not in main pipeline)

- **`ml/`** - K-mer taxonomy and classification
  - `kmer_taxonomy.rs` - K-mer-based taxonomic classification
  - `simple_classifier.rs` - Contig classification
  - `classification_reporter.rs` - Result formatting
  - `semibin_logger.rs` - Logging utilities

- **`core/`** - Core data structures
  - `data_structures.rs` - Graph, reads, contigs
  - `paired_reads.rs` - Paired-end read handling

- **`database/`** - SQLite backend for k-mers and taxonomy
  - `integration.rs` - Database operations

- **`features/`** - Sequence and graph feature extraction
  - `extraction.rs` - Feature extraction pipeline

- **`qc/`** - Quality control and filtering
  - `quality_filter.rs` - Read quality filtering
  - `adapter_trimmer.rs` - Adapter removal
  - `qc_stats.rs` - Quality metrics
  - `presets.rs` - QC presets

- **`pipeline/`** - Main pipeline orchestration
  - `complete_integration.rs` - Full analysis pipeline
  - `fast_pipeline.rs` - Fast mode implementation

- **`reporting/`** - Output generation
  - `report_generator.rs` - Analysis reports

- **`utils/`** - Shared utilities
  - `configuration.rs` - Config management
  - `kraken_reporter.rs` - Kraken2-format output
  - `output_writers.rs` - Format writers (TSV, JSON, FASTA, GFA)
  - `streaming_abundance.rs` - Streaming abundance estimation
  - `async_output_manager.rs` - Async output handling
  - `progress_display.rs` - Terminal UI
  - `genomic_validator.rs` - Sequence validation
  - `intermediate_output.rs` - Checkpoint management

### Key Components

#### Assembly Pipeline
```rust
// Laptop-optimized assembly
use meta_forge::assembly::laptop_assembly::LaptopAssembler;

let assembler = LaptopAssembler::new(config);
let contigs = assembler.assemble(reads)?;
```

#### ML Classification
```rust
// K-mer based classification
use meta_forge::ml::simple_classifier::SimpleContigClassifier;
use meta_forge::ml::kmer_taxonomy::KmerTaxonomyClassifier;

let classifier = SimpleContigClassifier::new(config)?;
let results = classifier.classify_contigs(&contigs)?;
```

#### Feature Extraction
```rust
use meta_forge::features::extraction::AdvancedFeatureExtractor;

let extractor = AdvancedFeatureExtractor::new(config)?;
let features = extractor.extract_sequence_features(sequence)?;
```

## Output Files

The pipeline generates:
- `classifications.tsv` - Classification results in TSV format (MetaBAT2/MaxBin2 compatible)
- `classification_summary.txt` - Human-readable summary with taxa distribution
- `contigs.fasta` - Assembled contigs with coverage information
- `corrected_reads.fastq` - Quality-filtered and error-corrected reads
- `features.tsv` - Extracted sequence and graph features
- `final_report.json` - Complete analysis report with metrics
- `assembly_graph.gfa` - Assembly graph in GFA format
- `kraken_report.txt` - Kraken2-compatible classification report

## Performance Optimizations

### Recent Improvements
1. **Zero-copy k-mer processing** (6-8x speedup)
   - Rolling hash for O(1) updates
   - Direct buffer writing without allocations

2. **Lock-free parallel graph construction** (3-4x speedup)
   - DashMap for concurrent node/edge operations
   - Atomic coverage updates

3. **Cache-optimized CSR graph** (4-6x speedup)
   - Compressed sparse row format
   - Sequential memory access patterns

### Memory Management
```rust
// Auto-detect system resources
let config = LaptopConfig::auto_detect();

// Or manual configuration
let config = LaptopConfig::custom(
    memory_budget_mb: 4096,
    cpu_cores: 8,
    max_k: 31
)?;
```

## Configuration

### Config File (config.toml)
```toml
[general]
output = "./results"
threads = 8
memory_mb = 4096

[assembly]
k_min = 21
k_max = 31
min_coverage = 3

[classification]
kmer_size = 4
min_contig_length = 500
num_bins = 10

[features]
use_ml = true
extract_graph_metrics = true
```

### Load Config
```bash
meta-forge -c config.toml analyze input.fastq
```

## Database Operations

### Initialize
```bash
meta-forge database init data/metadb.sqlite
```

### Import Taxonomy
```bash
meta-forge database import-taxonomy data/metadb.sqlite \
  --names nodes.dmp \
  --nodes names.dmp
```

### Build K-mer Index
```bash
meta-forge database build-index data/metadb.sqlite \
  --k 21 \
  --input-fasta refs.fasta
```

## Common Workflows

### Standard Analysis
```bash
meta-forge -m 4096 -j 8 analyze reads.fastq --sample-name project1
```

### High-Performance Mode
```bash
meta-forge -m 16384 -j 32 analyze reads.fastq --mode sensitive
```

### Assembly Only
```bash
meta-forge assemble reads.fastq --k-range 21-51 --min-coverage 3
```

### Resume from Checkpoint
```bash
meta-forge resume --run-id 20240102_123456 --section classification
```

## Troubleshooting

### Error: "Invalid input: '-m' looks like a command-line flag"
**Cause**: Flags placed after subcommand
**Fix**: Move flags before subcommand
```bash
# Wrong: meta-forge analyze input.fastq -m 4096
# Right: meta-forge -m 4096 analyze input.fastq
```

### Error: "No input files provided"
**Fix**: Ensure input file path is correct and file exists
```bash
ls -lh input.fastq  # Check file exists
meta-forge analyze input.fastq
```

### Performance Issues
1. Check system resources: `htop` or `Activity Monitor`
2. Adjust memory: `-m` flag (1024-8192 MB typical)
3. Adjust threads: `-j` flag (2-32 cores typical)
4. Use fast mode: `--mode fast`

## Contributing

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Make changes and test: `cargo test`
4. Commit: `git commit -m "Description"`
5. Push: `git push origin feature-name`
6. Create pull request

## Testing Guidelines

- Write tests for new features
- Ensure biological correctness
- Benchmark performance changes
- Document edge cases

## Code Cleanliness

This codebase follows strict organizational principles:

### Recently Removed (Oct 2025 Housekeeping)
- ❌ Dead CLI module (`src/cli/fast_mode.rs`) - unused
- ❌ Orphaned test files referencing non-existent modules (9 files removed)
- ❌ Unused ML dependencies (`candle-core`, `candle-nn`, `ort`, `smartcore`) - 100+ MB savings
- ❌ Duplicate dependencies consolidated
- ✅ Test dependencies moved to `[dev-dependencies]`

### File Organization Rules
- **Never save files to project root** - use subdirectories
- Source code → `/src`
- Tests → `/tests`
- Examples → `/examples`
- Benchmarks → `/benches`
- Documentation → `/docs` or module-level `.md` files

### Dependency Philosophy
- **Minimal**: Only include dependencies actively used in production
- **Performant**: Prefer zero-cost abstractions and low-overhead crates
- **Maintained**: All dependencies actively maintained
- **Consolidated**: No duplicate crate versions unless unavoidable

## Resources

- GitHub: https://github.com/mladen5000/meta-forge
- Issues: https://github.com/mladen5000/meta-forge/issues
- Documentation: `/docs` directory

---

**Remember**:
1. Always put global flags BEFORE the subcommand! ✅
2. Use `LaptopAssembler` for production assembly (not the `optimized/` module)
3. Keep dependencies lean - remove unused crates immediately
