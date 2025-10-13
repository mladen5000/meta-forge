# MetaForge Pipeline - Missing Components Analysis

## Executive Summary

MetaForge has a **functional pipeline** from raw reads to binned contigs, but it's missing the final critical component: **real taxonomy integration**. Currently, the pipeline outputs "Bin 0", "Bin 1", etc. instead of actual species names like "Escherichia coli", "Staphylococcus aureus".

## What Works ✅

### 1. **Complete Assembly Pipeline**
- ✅ QC & Error Correction (`src/qc/`)
- ✅ Laptop-optimized k-mer Assembly (`src/assembly/laptop_assembly.rs`)
- ✅ Contig Generation with coverage tracking
- ✅ Feature Extraction (k-mer profiles, GC content, coverage)

### 2. **ML-Based Binning**
- ✅ K-mer composition analysis (`src/ml/simple_classifier.rs`)
- ✅ Coverage-based clustering
- ✅ Adaptive bin count detection (prevents over-binning)
- ✅ Bin quality metrics (completeness, contamination estimates)

### 3. **Database Infrastructure**
- ✅ SQLite database with proper schema (`src/database/integration.rs`)
- ✅ Tables for: taxonomy, sequences, contigs, annotations, k-mer index
- ✅ Performance: WAL mode, indices, caching
- ✅ Methods for: `insert_taxonomy_entries()`, `get_taxonomy_name()`, `find_similar_sequences()`

### 4. **Output Generation**
- ✅ TSV format (MetaBAT2/MaxBin2 compatible)
- ✅ JSON reports with metrics
- ✅ FASTA contigs with coverage
- ✅ Kraken2-compatible reports (structure exists)

## What's Missing ❌

### **PRIMARY GAP: Real Taxonomy Assignment**

Currently, the pipeline outputs:
```
Bin 0: 45 contigs (1.2 MB)
Bin 1: 23 contigs (0.8 MB)
Bin 2: 12 contigs (0.3 MB)
```

**Should output:**
```
Escherichia coli: 45 contigs (1.2 MB, 32% abundance)
Staphylococcus aureus: 23 contigs (0.8 MB, 21% abundance)
Bacillus subtilis: 12 contigs (0.3 MB, 9% abundance)
```

### Missing Components Breakdown

#### 1. **Taxonomy Database Population** (CRITICAL)

**Location:** `src/database/integration.rs`

**What exists:**
- ✅ Database schema with `taxonomy` table
- ✅ Methods: `insert_taxonomy_entries()`, `get_taxonomy_name()`
- ✅ NCBI taxonomy import function: `DatabaseMigrator::import_taxonomy_from_ncbi()`

**What's missing:**
- ❌ **No real taxonomy data loaded** - database is empty
- ❌ No reference genome database (NCBI RefSeq, GTDB, etc.)
- ❌ No pre-built k-mer → taxon mapping

**Where it's used:**
- `src/ml/kmer_taxonomy.rs` - Has infrastructure but uses 4 hardcoded synthetic profiles:
  ```rust
  pub fn load_default_references(&mut self) -> Result<()> {
      // TODO: Load from real database
      self.add_example_reference("Escherichia coli", "species", 0.50);
      self.add_example_reference("Staphylococcus aureus", "species", 0.33);
      self.add_example_reference("Bacillus subtilis", "species", 0.44);
      self.add_example_reference("Pseudomonas aeruginosa", "species", 0.66);
      Ok(())
  }
  ```

#### 2. **Bin → Taxon Mapping** (CRITICAL)

**Current flow:**
```
Reads → Assembly → Contigs → Binning (ML) → Output "Bin 0, Bin 1..."
                                              ❌ STOPS HERE
```

**Needed flow:**
```
Reads → Assembly → Contigs → Binning (ML) → Taxonomy Assignment → Output "E. coli, S. aureus..."
                                              ✅ NEEDS THIS STEP
```

**Implementation needed:**
- Method to map bin_id → taxonomy_id → species name
- Options:
  1. **K-mer based** (fastest, already partially implemented)
  2. **BLAST/Diamond** (most accurate, slower)
  3. **Marker gene** (CheckM2 style, good balance)

**Where to add:**
- New module: `src/ml/taxonomy_assignment.rs`
- Or extend: `src/ml/kmer_taxonomy.rs::KmerTaxonomyClassifier`
- Integration point: `src/ml/simple_classifier.rs::classify_contigs()` → add taxonomy lookup

#### 3. **Abundance Profiling** (IMPORTANT)

**What exists:**
- ✅ Contig-level coverage tracking
- ✅ Bin-level size metrics
- ✅ `AbundanceProfile` struct in `complete_integration.rs` (unused)

**What's missing:**
- ❌ **Relative abundance calculation** (% of total reads/bases per taxon)
- ❌ Normalization by genome size
- ❌ Multi-sample comparison support

**Formula needed:**
```
Abundance(taxon) = (Σ contig_coverage * contig_length) / total_coverage
                   normalized by expected_genome_size(taxon)
```

**Where to add:**
- New module: `src/ml/abundance_estimator.rs`
- Or extend: `src/utils/streaming_abundance.rs` (partially implemented)

#### 4. **Reference Database Integration**

**Options (in order of complexity):**

**Option A: Minimal - GTDB Sketches** (Fastest to implement)
- Use pre-built k-mer sketches from GTDB (Genome Taxonomy Database)
- ~300K genomes, taxonomy-aware
- Implementation: Sourmash or similar
- Size: ~2-5 GB
- **Time to implement: 1-2 days**

**Option B: Medium - Custom K-mer DB** (Balance)
- Build k-mer index from RefSeq representative genomes
- ~200 common pathogens/commensals
- Implementation: Use existing `database/integration.rs::build_kmer_index()`
- Size: ~500 MB - 2 GB
- **Time to implement: 3-5 days**

**Option C: Full - NCBI RefSeq** (Most comprehensive)
- Complete bacterial/viral reference genomes
- ~250K genomes
- Implementation: Full BLAST/Diamond pipeline
- Size: 50+ GB
- **Time to implement: 1-2 weeks**

## Implementation Roadmap

### Phase 1: Minimal Viable Taxonomy (1-3 days)

**Goal:** Replace "Bin 0" with real species names

```rust
// 1. Add taxonomy lookup to classification
// File: src/ml/simple_classifier.rs

impl SimpleContigClassifier {
    pub fn classify_contigs_with_taxonomy(
        &self,
        contigs: &[Contig],
        taxonomy_db: &MetagenomicsDatabase,
    ) -> Result<Vec<TaxonomicContigClassification>> {
        // Existing binning
        let bin_classifications = self.classify_contigs(contigs)?;

        // NEW: Taxonomy assignment per bin
        let taxonomy_classifier = KmerTaxonomyClassifier::new(4);
        taxonomy_classifier.load_reference_database(taxonomy_db)?;

        let mut taxon_classifications = Vec::new();
        for classification in bin_classifications {
            let contig = &contigs[classification.contig_id];
            let taxon = taxonomy_classifier.classify_sequence(&contig.sequence)?;

            taxon_classifications.push(TaxonomicContigClassification {
                contig_id: classification.contig_id,
                bin_id: classification.bin_id,
                taxon: taxon.taxon,
                rank: taxon.rank,
                confidence: classification.confidence * taxon.confidence,
                features: classification.features,
            });
        }

        Ok(taxon_classifications)
    }
}
```

**Tasks:**
1. ✅ Implement `KmerTaxonomyClassifier::load_reference_database()`
2. ✅ Populate database with GTDB taxonomy (download & import)
3. ✅ Wire into main pipeline (`complete_integration.rs`)
4. ✅ Update output writers to use taxon names instead of bin_id

### Phase 2: Abundance Profiling (2-3 days)

```rust
// File: src/ml/abundance_estimator.rs (NEW)

pub struct AbundanceEstimator {
    genome_sizes: AHashMap<String, usize>, // taxon -> expected genome size
}

impl AbundanceEstimator {
    pub fn estimate_abundances(
        &self,
        taxon_classifications: &[TaxonomicContigClassification],
        contigs: &[Contig],
    ) -> Result<Vec<TaxonAbundance>> {
        let mut taxon_coverages: AHashMap<String, f64> = AHashMap::new();

        for classification in taxon_classifications {
            let contig = &contigs[classification.contig_id];
            let coverage_contribution = contig.coverage * contig.length as f64;

            *taxon_coverages.entry(classification.taxon.clone()).or_insert(0.0)
                += coverage_contribution;
        }

        // Normalize by genome sizes and calculate relative abundances
        let total_coverage: f64 = taxon_coverages.values().sum();

        let mut abundances = Vec::new();
        for (taxon, coverage) in taxon_coverages {
            let genome_size = self.genome_sizes.get(&taxon).copied().unwrap_or(3_000_000);
            let normalized_coverage = coverage / genome_size as f64;
            let relative_abundance = (normalized_coverage / total_coverage) * 100.0;

            abundances.push(TaxonAbundance {
                taxon,
                relative_abundance,
                read_count: coverage as usize,
                num_contigs: taxon_classifications.iter()
                    .filter(|c| c.taxon == taxon)
                    .count(),
            });
        }

        abundances.sort_by(|a, b| b.relative_abundance.partial_cmp(&a.relative_abundance).unwrap());
        Ok(abundances)
    }
}
```

**Tasks:**
1. ✅ Implement abundance calculation
2. ✅ Add genome size database (NCBI genome_sizes.tsv)
3. ✅ Integrate into pipeline
4. ✅ Update reports with abundance %

### Phase 3: Advanced Features (Optional, 1-2 weeks)

1. **Multi-sample comparison**
   - Beta diversity metrics
   - PCA/NMDS visualization
   - Differential abundance testing

2. **Strain-level resolution**
   - SNP calling within species
   - Strain tracking across samples

3. **Functional profiling**
   - Gene prediction (Prodigal)
   - KEGG/COG annotation
   - Metabolic pathway inference

## Data Requirements

### Minimal Setup (Phase 1)

1. **GTDB Taxonomy** (~100 MB)
   ```bash
   wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv
   ./target/release/meta-forge database import-taxonomy \
     --taxonomy bac120_taxonomy.tsv
   ```

2. **Representative Genomes** (2-5 GB)
   ```bash
   # Download GTDB representative genomes (R220)
   wget https://data.gtdb.ecogenomic.org/releases/latest/genomic_files_reps/gtdb_genomes_reps.tar.gz

   # Build k-mer index
   ./target/release/meta-forge database build-index \
     --k 21 --input-fasta gtdb_genomes_reps/*.fna
   ```

### Full Setup (Phase 3)

1. **NCBI RefSeq** (50+ GB)
2. **CheckM2 database** (3 GB)
3. **KEGG/COG databases** (5-10 GB)

## Example Output Comparison

### Current Output ❌
```
Classification Results
=====================
Total contigs: 180
Bins created: 10

Bin 0: 45 contigs (1.2 MB)
Bin 1: 23 contigs (0.8 MB)
Bin 2: 18 contigs (0.6 MB)
...
```

### Target Output ✅
```
Taxonomic Classification Results
================================
Total contigs: 180
Species identified: 6

Escherichia coli (GTDB: s__Escherichia_coli)
  - Abundance: 32.4%
  - Contigs: 45 (1.2 MB)
  - Completeness: 94.2%
  - Contamination: 1.3%
  - N50: 28,456 bp

Staphylococcus aureus (GTDB: s__Staphylococcus_aureus)
  - Abundance: 21.8%
  - Contigs: 23 (0.8 MB)
  - Completeness: 87.6%
  - Contamination: 2.1%
  - N50: 34,892 bp

Bacillus subtilis (GTDB: s__Bacillus_subtilis)
  - Abundance: 15.3%
  - Contigs: 18 (0.6 MB)
  - Completeness: 91.3%
  - Contamination: 0.8%
  - N50: 41,223 bp

[... 3 more species ...]

Unclassified: 38.5% (71 contigs, 2.1 MB)
```

## Priority Actions

### Immediate (This Week)
1. **Download GTDB taxonomy** → Populate database
2. **Implement bin → taxon mapping** → `src/ml/taxonomy_assignment.rs`
3. **Wire taxonomy into pipeline** → Update `complete_integration.rs`
4. **Update output writers** → Replace bin_id with taxon names

### Short Term (Next 2 Weeks)
1. **Implement abundance calculation**
2. **Add genome size normalization**
3. **Generate Kraken2-compatible reports**
4. **Add visualization (bar charts, krona plots)**

### Long Term (1-2 Months)
1. **Multi-sample support**
2. **Strain-level tracking**
3. **Functional profiling**
4. **Web dashboard**

## Code Locations Reference

### Key Files to Modify

1. **`src/ml/simple_classifier.rs`** (lines 230-428)
   - Add taxonomy lookup after binning
   - Return `TaxonomicContigClassification` instead of `ContigClassification`

2. **`src/ml/kmer_taxonomy.rs`** (lines 79-96)
   - Replace `load_default_references()` with database loading
   - Add method: `load_reference_database(db: &MetagenomicsDatabase)`

3. **`src/database/integration.rs`** (lines 375-405, 562-583)
   - Already has taxonomy methods - just need to populate data
   - Add method: `get_taxonomy_by_kmer_similarity()`

4. **`src/pipeline/complete_integration.rs`** (lines 2300-2500)
   - Update classification call to include taxonomy
   - Modify output generation to use taxon names

5. **`src/utils/kraken_reporter.rs`**
   - Already structured for Kraken2 output
   - Needs wiring to taxonomy results

## Conclusion

The pipeline is **90% complete**. The missing 10% (taxonomy integration) is the most visible part to users. With the database infrastructure already in place, implementing real taxonomy assignment is a **1-3 day task** for minimal functionality, extending to 1-2 weeks for production-quality with abundance profiling.

The hardest part (assembly, binning, database design) is done. The remaining work is mostly **data integration and plumbing**.
