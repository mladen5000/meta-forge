# MetaForge Development Guide

## Project Overview

MetaForge is a high-performance Rust-based metagenomics pipeline featuring:
- Adaptive k-mer assembly with graph optimization
- ML-based taxonomic classification
- Real-time error correction
- Comprehensive feature extraction
- Modern terminal UI with progress tracking

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
- `assembly/` - k-mer assembly and graph construction
- `ml/` - Machine learning classification
- `database/` - SQLite backend for k-mers and taxonomy
- `features/` - Sequence and graph feature extraction
- `pipeline/` - Main pipeline orchestration
- `utils/` - Shared utilities

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
- `classifications.tsv` - Classification results (TSV format)
- `classification_summary.txt` - Human-readable summary
- `contigs.fasta` - Assembled contigs
- `corrected_reads.fastq` - Error-corrected reads
- `final_report.json` - Complete analysis report
- `assembly_graph.gfa` - Assembly graph (GFA format)

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

## Resources

- GitHub: https://github.com/mladen5000/meta-forge
- Issues: https://github.com/mladen5000/meta-forge/issues
- Documentation: `/docs` directory

---

**Remember**: Always put global flags BEFORE the subcommand! ✅
