# MetaForge – Rust Metagenomics Pipeline

&#x20;&#x20;

MetaForge began as a side project in a small lab that needed a faster way to assemble and classify metagenomic reads. Over time it grew into a compact Rust toolkit that handles everything from k‑mer assembly to basic machine‑learning–assisted taxonomy, all without drowning you in dependencies.

---

## What It Does

* **Adaptive k‑mers**: Picks k‑mer sizes based on local complexity, so you don’t run three assemblies to find the sweet spot.
* **Pragmatic ML**: A lightweight neural model helps resolve repeats and tricky taxa, but core tasks stay in Rust for speed and reproducibility.
* **Feature extraction**: Compute sequence stats and graph metrics in one go, then dump to JSON/CSV for downstream analysis.
* **SQLite backend**: Simple embedded database for k‑mer lookups and taxonomy tables—no need to install PostgreSQL.
* **Streamable & multithreaded**: Process large FASTQ files in chunks, use as many cores as you have, and keep memory in check.
* **Minimal setup**: One `cargo build --release`, one binary, no Python or R required.

---

## Quickstart

1. Clone and build:

   ```bash
   git clone https://github.com/mladen5000/meta-forge.git
   cd meta-forge
   cargo build --release
   ```

   The `meta-forge` binary will land in `target/release/`.

2. Point it at your FASTQ:

   ```bash
   ./target/release/meta-forge analyze reads.fastq --sample-name sample1 --mode standard --threads 4
   ```

3. Check `./results` (or whatever folder you set) for:

   * `contigs.fasta`
   * `assembly_stats.json`
   * `taxonomy.json`
   * `abundance_profile.csv`
   * performance logs in `performance.json`

---

## Modes & Examples

### Standard analysis

```bash
meta-forge analyze reads.fastq \
  --sample-name projectX \
  --output results/ \
  --threads 8
```

### High‑throughput (more RAM/cores)

```bash
meta-forge analyze big_dataset.fastq \
  --sample-name bigRun \
  --mode high-performance \
  --memory 32 \
  --threads 16 \
  --output big_results/
```

### Low‑memory (for laptops)

```bash
meta-forge analyze reads.fastq \
  --sample-name lightRun \
  --mode low-memory \
  --memory 4 \
  --threads 2 \
  --output lm_results/
```

### Assembly‑only

```bash
meta-forge assemble reads.fastq --k-range 21-31 --min-coverage 4 --output asm_only/
```

---

## Taxonomy Database

Initialize DB:

```bash
meta-forge database init data/metadb.sqlite
```

Import NCBI dump files:

```bash
meta-forge database import-taxonomy data/metadb.sqlite \
  --names nodes.dmp \
  --nodes names.dmp
```

Build k‑mer index:

```bash
meta-forge database build-index data/metadb.sqlite --k 21 --input-fasta refs.fasta
```

Query examples:

```bash
meta-forge database query data/metadb.sqlite --sequence ATCGATCG --k 21 --output dbhits.json
```

---

## Config template

A basic `config.toml` lives under `configs/`; tweak thread count, memory, k‑ranges, and feature flags there.

Example snippet:

```toml
[general]
output = "./results"
threads = 4

[assembly]
k_min = 21
k_max = 31
min_coverage = 3

[features]
use_ml = true
extract_graph_metrics = true
```

---

## Contributing

This is a small open‑source project maintained by a former grad student. Issues, PRs, and real‑world test datasets welcome: [https://github.com/mladen5000/meta-forge](https://github.com/mladen5000/meta-forge)

© 2024 Mladen Rasic. MIT License.
