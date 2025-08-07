# MetaForge – Your Rust-Powered Metagenomics Sidekick

[![Rust](https://img.shields.io/badge/rust-1.70+-blue.svg)](https://www.rust-lang.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

Hey there—welcome to MetaForge, the no-fluff metagenomics toolkit written in Rust. Whether you’re piecing together genomes, classifying species, or squeezing out every drop of performance, MetaForge is built to keep up without demanding a PhD in pipeline assembly.

---

## Why You’ll Love MetaForge

* **Smart K-mer Sizes**: It adjusts k-mer lengths on the fly so you don’t have to guess what works best.
* **ML Where It Matters**: Machine-learning models tackle tricky classification and repeat regions, but we don’t use AI for the basics—those are fast enough.
* **Deep Feature Mining**: Get both classic sequence stats and graph-based insights in one go.
* **Speed That Scales**: Multithreaded from day one. Memory-hungry? Only if you ask for it.
* **Built-In Database**: An embedded SQLite store for k-mer lookups and taxonomy—no extra dependencies.
* **Live Feedback**: See progress bars, performance metrics, and logs as your job runs.
* **Config Your Way**: TOML files with sensible defaults, plus a few presets to get started.

---

## Getting Started

### Prerequisites

* Rust ≥ 1.70
* SQLite3
* OpenMP (optional, for parallel bits)

### Clone & Build

```bash
git clone https://github.com/mladen5000/meta-forge.git
cd meta-forge
cargo build --release
# You’ll find the `meta-forge` binary in target/release/
```

Or, if you prefer one-liners:

```bash
cargo install --path .
```

---

## First Run

Analyze a FASTQ in standard mode:

```bash
meta-forge analyze my_reads.fastq --sample-name project1 --mode standard
```

Need a config template? We got you:

```bash
meta-forge config standard --output my_project.toml
```

---

## Common Workflows

### Standard Analysis

```bash
meta-forge analyze reads.fastq \
  --sample-name sampleA \
  --output ./results \
  --threads 8
```

### Big Data, High Performance

```bash
meta-forge analyze big.fastq \
  --sample-name bigA \
  --mode high-performance \
  --memory 16 \
  --threads 16 \
  --output ./big_results
```

### Low-Memory Mode (Tiny Machines)

```bash
meta-forge analyze reads.fastq --sample-name lowmem \
  --mode low-memory --memory 4 --threads 4
```

### Paired-End Reads

```bash
meta-forge analyze R1.fastq R2.fastq --sample-name pairX --output ./pair_results
```

### Just Assembly

```bash
meta-forge assemble reads.fastq --k-range 21-31 --min-coverage 3 --output ./assembly
```

---

## Your Data, Your DB

Initialize a fresh database:

```bash
meta-forge database init ./data/metagenomics.db
```

Load NCBI taxonomy:

```bash
meta-forge database import-taxonomy ./data/metagenomics.db \
  --taxonomy-file ncbi_taxonomy.txt \
  --names-file names.dmp \
  --nodes-file nodes.dmp
```

Build your k-mer index:

```bash
meta-forge database build-index ./data/metagenomics.db \
  --k 21 \
  --input-fasta ref_genomes.fasta \
  --batch-size 10000
```

---

## Configuration at a Glance

Drop this in `config.toml` to get rolling:

```toml
[general]
working_directory = "./data"
output_directory  = "./results"
max_threads       = 0  # auto-detect cores

[assembly]
k_min               = 21
k_max               = 31
k_step              = 2
min_coverage        = 3
min_length          = 200
graph_simplification = true

[features]
composition_features = true
complexity_features  = true
topology_features    = true
feature_dimensions   = 100

[database]
connection_string = "sqlite:./data/metagenomics.db"
batch_size        = 1000
index_k_size      = 21

[performance]
memory_limit_gb   = 8
enable_monitoring = true
checkpoint_interval = 1000
```

---

## What You’ll See

**Results** land in your `--output` folder with clear names:

* `contigs.fasta` + `assembly_stats.json` + `assembly_graph.gfa`
* `features.json` / `features.csv` + interactive HTML summary
* `taxonomy.json` + `abundance_profile.csv` + HTML report
* `performance.json` + `memory_profile.json`

---

## Tips & Tricks

* **OOM?** Drop into low-memory mode:
  `meta-forge analyze reads.fastq --mode low-memory --memory 2`

* **Feeling slow?** Crank threads or hit high-performance:
  `meta-forge analyze reads.fastq --mode high-performance --threads 16`

* **Database hiccup?**
  `meta-forge database init ./data/metagenomics.db`

* **Need more logs?**
  `meta-forge analyze reads.fastq --verbose --debug`

---

## Contributing & Support

Spotted a bug or missing feature? Open an issue or PR:
[https://github.com/mladen5000/meta-forge](https://github.com/mladen5000/meta-forge)

Need help? Jump into Discussions or mail me: [youremail@example.com](mailto:youremail@example.com)

---

Licensed under MIT. Have fun forging genomes! 🚀
