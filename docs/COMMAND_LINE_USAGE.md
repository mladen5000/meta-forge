# Command Line Usage - Meta-Forge Pipeline

## ⚠️ IMPORTANT: Flag Order Matters!

The CLI parser requires **global flags BEFORE the subcommand**.

---

## ✅ CORRECT Command Format

### Structure
```bash
meta-forge [GLOBAL_FLAGS] <SUBCOMMAND> [SUBCOMMAND_ARGS]
```

### Global Flags (MUST come first)
- `-m <MB>` or `--memory-mb <MB>` - Memory budget in MB
- `-j <N>` or `--threads <N>` - Number of threads
- `-c <FILE>` or `--config <FILE>` - Configuration file
- `-v` or `--verbose` - Verbose logging
- `-o <DIR>` or `--output <DIR>` - Output directory
- `--auto-detect` - Auto-detect system resources (default: true)

---

## Common Commands

### 1. Run Complete Analysis
```bash
# ✅ CORRECT - Flags BEFORE 'analyze'
meta-forge -m 4096 -j 8 analyze input.fastq

# ✅ CORRECT - With sample name
meta-forge -m 4096 -j 8 analyze input.fastq --sample-name my_sample

# ✅ CORRECT - With mode selection
meta-forge -m 4096 analyze input.fastq --mode fast

# ❌ WRONG - Flags AFTER 'analyze'
meta-forge analyze input.fastq -m 4096  # This will fail!
```

### 2. Assembly Only
```bash
# ✅ CORRECT
meta-forge -m 2048 assemble input.fastq

# ✅ With k-mer range
meta-forge -m 2048 assemble input.fastq --k-range 21-51

# ✅ With minimum coverage
meta-forge assemble input.fastq --min-coverage 3
```

### 3. Feature Extraction Only
```bash
# ✅ CORRECT
meta-forge features input.fasta --types composition,kmer

# ✅ With output format
meta-forge features input.fasta --types composition --format csv
```

### 4. Database Operations
```bash
# Initialize database
meta-forge database init --db-path ./database

# Add reference
meta-forge database add-reference genome.fasta --db-path ./database

# Query database
meta-forge database query sequence.fasta --db-path ./database
```

### 5. Configuration
```bash
# Generate config template
meta-forge config standard --output my_config.toml

# Use config file
meta-forge -c my_config.toml analyze input.fastq
```

---

## Analysis Modes

Available via `-M` or `--mode` flag in `analyze` command:

- `standard` - Full analysis (default)
- `fast` - Quick analysis with reduced accuracy
- `sensitive` - More sensitive taxonomic assignment
- `assembly-only` - Assembly without classification
- `classification-only` - Classify pre-assembled contigs

**Example**:
```bash
meta-forge -m 4096 analyze input.fastq --mode fast
```

---

## Memory & Thread Settings

### Auto-detection (Recommended)
```bash
# Uses --auto-detect (default: true)
meta-forge analyze input.fastq
```

### Manual Configuration
```bash
# 4GB memory, 8 threads
meta-forge -m 4096 -j 8 analyze input.fastq

# 2GB memory, 4 threads
meta-forge -m 2048 -j 4 analyze input.fastq
```

### Memory Budget Guidelines
- **4GB laptop**: `-m 1024` (1GB for analysis)
- **8GB laptop**: `-m 2048` (2GB for analysis)
- **16GB+ laptop**: `-m 4096` (4GB for analysis)
- **Server**: `-m 8192` or higher

---

## Output Options

### Specify Output Directory
```bash
meta-forge -o ./results analyze input.fastq
```

### Output Files Created
```
output/
├── classifications.tsv           ← Classification results (TSV)
├── classification_summary.txt    ← Human-readable summary
├── contigs.fasta                 ← Assembled contigs
├── corrected_reads.fastq         ← Error-corrected reads
├── final_report.json             ← Complete analysis report
└── assembly_graph.gfa            ← Assembly graph (GFA format)
```

---

## Common Error Fixes

### Error: "Unsupported file format: -m"
**Cause**: Flags placed AFTER subcommand
**Fix**: Move flags BEFORE subcommand
```bash
# ❌ WRONG
meta-forge analyze input.fastq -m 4096

# ✅ CORRECT
meta-forge -m 4096 analyze input.fastq
```

### Error: "No input files provided"
**Cause**: Missing input file argument
**Fix**: Provide input file path
```bash
meta-forge -m 4096 analyze input.fastq
```

### Error: "Failed to parse k-range"
**Cause**: Invalid k-mer range format
**Fix**: Use format `MIN-MAX`
```bash
meta-forge assemble input.fastq --k-range 21-51
```

---

## Quick Reference

### Minimal Command (Auto-detect Everything)
```bash
meta-forge analyze input.fastq
```

### Typical Laptop Command
```bash
meta-forge -m 2048 -j 4 analyze input.fastq --sample-name my_sample
```

### High-Performance Server
```bash
meta-forge -m 16384 -j 32 analyze input.fastq --mode sensitive
```

### Resume from Checkpoint
```bash
meta-forge resume --run-id 20240102_123456 --section classification
```

### List Previous Runs
```bash
meta-forge list-runs
```

---

## Help Commands

### General Help
```bash
meta-forge --help
```

### Subcommand Help
```bash
meta-forge analyze --help
meta-forge assemble --help
meta-forge features --help
meta-forge database --help
```

### Version
```bash
meta-forge --version
```

---

## Examples with Real Data

### Example 1: Simple Analysis
```bash
# Basic analysis with auto-detection
meta-forge analyze sample1.fastq
```

### Example 2: Multiple Input Files
```bash
# Paired-end reads
meta-forge -m 4096 analyze sample_R1.fastq sample_R2.fastq
```

### Example 3: With All Options
```bash
meta-forge \
  -m 4096 \
  -j 8 \
  -o ./results \
  -v \
  analyze sample.fastq \
  --sample-name metagenome_01 \
  --mode sensitive \
  --database ./refdb
```

### Example 4: Fast Mode for Quick Check
```bash
# Quick analysis (5-10x faster, slightly less accurate)
meta-forge -m 2048 analyze sample.fastq --mode fast
```

---

## Performance Tips

1. **Use auto-detect** for first run: `meta-forge analyze input.fastq`
2. **Adjust memory** based on dataset size:
   - Small (<1M reads): `-m 1024`
   - Medium (1-10M reads): `-m 2048`
   - Large (>10M reads): `-m 4096+`
3. **Use fast mode** for exploration: `--mode fast`
4. **Enable threading** on multi-core systems: `-j <cores>`
5. **Resume failed runs**: `meta-forge resume --run-id <id> --section <section>`

---

## Troubleshooting

### Check System Resources
```bash
# Memory
free -h  # Linux
vm_stat | grep free  # macOS

# CPU cores
nproc  # Linux
sysctl -n hw.ncpu  # macOS
```

### Monitor Pipeline Progress
```bash
# Verbose mode shows detailed progress
meta-forge -v analyze input.fastq
```

### Check Output Files
```bash
ls -lh output/
cat output/classification_summary.txt
column -t -s $'\t' output/classifications.tsv | less
```

---

## Flag Reference Table

| Flag | Short | Long | Description | Example |
|------|-------|------|-------------|---------|
| Memory | `-m` | `--memory-mb` | Memory budget (MB) | `-m 4096` |
| Threads | `-j` | `--threads` | CPU threads | `-j 8` |
| Config | `-c` | `--config` | Config file | `-c config.toml` |
| Verbose | `-v` | `--verbose` | Debug logging | `-v` |
| Output | `-o` | `--output` | Output directory | `-o ./results` |
| Auto | - | `--auto-detect` | Auto system detection | `--auto-detect` |

---

**Remember**: Global flags MUST come BEFORE the subcommand! ✅
