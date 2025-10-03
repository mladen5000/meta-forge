# MetaForge – Rust Metagenomics Pipeline

A high-performance, laptop-friendly metagenomics analysis toolkit built in Rust. MetaForge handles everything from adaptive k-mer assembly to ML-based taxonomic classification with a modern terminal UI.

---

## Features

* **Adaptive k-mer assembly** - Intelligently selects k-mer sizes based on sequence complexity
* **ML-enhanced classification** - K-mer taxonomy + neural clustering for accurate species identification
* **Interactive TUI** - Real-time monitoring with progress tracking and live metrics
* **SQLite backend** - Embedded database for k-mer indices and taxonomy (no external database)
* **Optimized performance** - Lock-free parallelization, zero-copy k-mer processing, cache-friendly data structures
* **Production-ready** - Comprehensive error handling, checkpointing, and resume capabilities

---

## Quick Start

### 1. Build from source
```bash
git clone https://github.com/mladen5000/meta-forge.git
cd meta-forge
cargo build --release
```

Binary will be at `target/release/meta-forge`

### 2. Run analysis
```bash
# Basic analysis (auto-detects system resources)
./target/release/meta-forge analyze reads.fastq

# With memory and thread control
./target/release/meta-forge -m 4096 -j 8 analyze reads.fastq --sample-name my_sample

# Full command with all options
./target/release/meta-forge -m 4096 -j 8 -o ./results analyze reads.fastq \
  --sample-name project1 \
  --mode standard
```

### 3. Check outputs
Results saved to `./output/` (or specified directory):
- `classifications.tsv` - Classification results (TSV format)
- `classification_summary.txt` - Human-readable summary
- `contigs.fasta` - Assembled contigs
- `corrected_reads.fastq` - Error-corrected reads
- `final_report.json` - Complete analysis report
- `assembly_graph.gfa` - Assembly graph (GFA format)

---

## ⚠️ Important: Command-Line Usage

### Global Flags MUST Come BEFORE Subcommand

```bash
# ✅ CORRECT - Flags before subcommand
meta-forge -m 4096 -j 8 analyze input.fastq

# ❌ WRONG - Will cause error
meta-forge analyze input.fastq -m 4096
```

### Global Flags
- `-m <MB>` or `--memory-mb <MB>` - Memory budget in MB (default: auto-detect)
- `-j <N>` or `--threads <N>` - Number of threads (default: auto-detect)
- `-o <DIR>` or `--output <DIR>` - Output directory (default: `./output`)
- `-c <FILE>` or `--config <FILE>` - Configuration file
- `-v` or `--verbose` - Verbose logging

### Subcommands
- `analyze` - Full analysis pipeline (default)
- `assemble` - Assembly only
- `features` - Feature extraction only
- `database` - Database operations
- `config` - Generate config template
- `resume` - Resume from checkpoint
- `list-runs` - List previous runs

---

## Usage Examples

### Standard Analysis
```bash
meta-forge -m 4096 -j 8 analyze reads.fastq --sample-name sample1
```

### High-Performance (Server)
```bash
meta-forge -m 16384 -j 32 analyze reads.fastq \
  --sample-name bigdata \
  --mode sensitive
```

### Low-Memory (Laptop)
```bash
meta-forge -m 2048 -j 4 analyze reads.fastq \
  --sample-name laptop_run \
  --mode fast
```

### Assembly Only
```bash
meta-forge assemble reads.fastq \
  --k-range 21-51 \
  --min-coverage 3 \
  -o assembly_output/
```

### Multiple Input Files
```bash
# Paired-end reads
meta-forge -m 4096 analyze sample_R1.fastq sample_R2.fastq
```

---

## Analysis Modes

- `standard` - Full analysis with balanced speed/accuracy (default)
- `fast` - Quick analysis, reduced accuracy (5-10x faster)
- `sensitive` - More sensitive taxonomic assignment (slower)
- `assembly-only` - Assembly without classification
- `classification-only` - Classify pre-assembled contigs

**Example:**
```bash
meta-forge -m 4096 analyze reads.fastq --mode fast
```

---

## Database Setup

### Initialize Database
```bash
meta-forge database init data/metadb.sqlite
```

### Import NCBI Taxonomy
```bash
meta-forge database import-taxonomy data/metadb.sqlite \
  --names nodes.dmp \
  --nodes names.dmp
```

### Build K-mer Index
```bash
meta-forge database build-index data/metadb.sqlite \
  --k 21 \
  --input-fasta references.fasta
```

### Query Database
```bash
meta-forge database query data/metadb.sqlite \
  --sequence ATCGATCGATCG \
  --k 21 \
  --output results.json
```

---

## Configuration

### Generate Config Template
```bash
meta-forge config standard --output my_config.toml
```

### Example Config (config.toml)
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

### Use Config
```bash
meta-forge -c my_config.toml analyze reads.fastq
```

---

## Output Format

### Classification TSV
Tab-separated values, compatible with MetaBAT2/MaxBin2/Anvi'o:
```
contig_id    bin_id    taxonomy_name           confidence    lineage                          method
1            100       Escherichia coli        0.9500       Bacteria;Escherichia coli        hybrid_kmer_taxonomy
2            101       Salmonella enterica     0.8700       Bacteria;Salmonella enterica     ml_kmer_clustering
```

### Classification Summary
Human-readable text format:
```
# Classification Summary
Total Contigs: 150

## Taxa Distribution
Escherichia coli: 75 contigs (50.0%)
Salmonella enterica: 45 contigs (30.0%)
Bin_3: 30 contigs (20.0%)

## Classification Methods
hybrid_kmer_taxonomy: 120 contigs (80.0%)
ml_kmer_clustering: 30 contigs (20.0%)
```

---

## Performance Tips

### Memory Guidelines
- **4GB laptop**: `-m 1024` (1GB for analysis)
- **8GB laptop**: `-m 2048` (2GB for analysis)
- **16GB+ laptop**: `-m 4096` (4GB for analysis)
- **Server**: `-m 8192+` (8GB+ for analysis)

### Speed Optimization
1. Use `--mode fast` for quick exploration
2. Increase threads on multi-core: `-j <cores>`
3. Enable auto-detection (default): remove `-m` and `-j` flags
4. Use SSD for output directory

### Resume Failed Runs
```bash
# List available runs
meta-forge list-runs

# Resume specific run
meta-forge resume --run-id 20240102_123456 --section classification
```

---

## Troubleshooting

### Error: "Invalid input: '-m' looks like a command-line flag"
**Cause**: Flags placed after subcommand
**Fix**: Move flags before subcommand
```bash
# Wrong: meta-forge analyze input.fastq -m 4096
# Right: meta-forge -m 4096 analyze input.fastq
```

### Error: "No input files provided"
**Cause**: Missing input file
**Fix**: Provide valid input file path
```bash
# Check file exists
ls -lh input.fastq

# Run with correct path
meta-forge analyze input.fastq
```

### Error: "Failed to parse k-range"
**Cause**: Invalid k-mer range format
**Fix**: Use `MIN-MAX` format
```bash
meta-forge assemble input.fastq --k-range 21-51
```

### Low Performance
1. Check system resources: `htop` (Linux) or Activity Monitor (macOS)
2. Increase memory: `-m 4096` or higher
3. Increase threads: `-j 8` or more
4. Use faster mode: `--mode fast`

---

## Help Commands

```bash
# General help
meta-forge --help

# Subcommand help
meta-forge analyze --help
meta-forge assemble --help
meta-forge database --help

# Version
meta-forge --version
```

---

## Development

### Build for Development
```bash
cargo build
```

### Run Tests
```bash
cargo test
```

### Run Benchmarks
```bash
cargo bench
```

### Profile Performance
```bash
cargo install flamegraph
cargo flamegraph --bin meta-forge -- analyze input.fastq
```

---

## Contributing

Issues, PRs, and test datasets welcome!

1. Fork the repository
2. Create feature branch: `git checkout -b feature-name`
3. Make changes and test: `cargo test`
4. Commit: `git commit -m "Description"`
5. Push: `git push origin feature-name`
6. Create pull request

GitHub: https://github.com/mladen5000/meta-forge

---

## Citation

If you use MetaForge in your research, please cite:

```
MetaForge: A High-Performance Rust Pipeline for Metagenomic Analysis
Mladen Rasic, 2024
https://github.com/mladen5000/meta-forge
```

---

## License

MIT License

© 2024 Mladen Rasic

---

**Quick Reference:**

```bash
# Minimal (auto-detect)
meta-forge analyze input.fastq

# Typical laptop
meta-forge -m 2048 -j 4 analyze input.fastq --sample-name my_sample

# High-performance server
meta-forge -m 16384 -j 32 analyze input.fastq --mode sensitive
```

**Remember: Global flags BEFORE subcommand!** ✅
