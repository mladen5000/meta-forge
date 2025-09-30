# Metagenomic Assembly Parameters - Quick Reference

**Last Updated:** 2025-09-30
**Source:** Research from 10+ peer-reviewed papers (2017-2025)

---

## MetaSPAdes Standard Parameters

| Parameter | Value | Justification | Meta-Forge Implementation |
|-----------|-------|--------------|---------------------------|
| **Min node coverage** | 2-3x | Remove sequencing errors | ✅ `laptop_assembly.rs:439,530` |
| **Min contig coverage** | 2.0x | MetaSPAdes standard | ✅ `laptop_assembly.rs:922` |
| **Min contig length** | 3×k (63bp @ k=21) | Minimum valid merged path | ✅ `laptop_assembly.rs:921` |
| **Tip removal threshold** | 2×k (42bp @ k=21) | Remove dead-end branches | ✅ `laptop_assembly.rs:534` |
| **Min path length** | 3 k-mers | Prevent spurious contigs | ✅ `laptop_assembly.rs:1023` |
| **Isolated h-path removal** | 200 bp | MetaSPAdes standard | ⚠️ Could add |
| **Coverage ratio threshold** | β = 10 | Edge classification | ⚠️ Could add |
| **Long edge threshold** | 1500 bp | Filigree edge detection | ⚠️ Could add |
| **Bulge collapse length** | ≤ 150 bp | Strain masking | ⚠️ Could add |

---

## K-mer Size Recommendations

### Single K-mer Assembly

| Scenario | Recommended K | Rationale |
|----------|--------------|-----------|
| **Bacterial genomes (standard)** | 21-31 | Balance error/repeat resolution |
| **High error rate** | 15-21 | More error tolerance |
| **High repeat content** | 41-63 | Better repeat resolution |
| **Low coverage** | 15-21 | Fewer coverage gaps |
| **Viral genomes** | 15-25 | Smaller genomes, more variation |
| **Fungal genomes** | 31-55 | Larger genomes, more repeats |

### Multi-K Assembly (Recommended)

| Tool | K-mer Range | Default Values |
|------|------------|----------------|
| **SPAdes** | 21-77 | [21, 33, 55, 77] |
| **IDBA-UD** | 20-100 | 20 to 100 (step 20) |
| **MEGAHIT** | 21-141 | [21,29,39,59,79,99,119,141] |
| **Meta-Forge (suggested)** | 15-63 | [21, 29, 39, 59] |

---

## Coverage Thresholds by Use Case

### Graph Construction Phase

| Purpose | Threshold | When to Apply |
|---------|-----------|--------------|
| **K-mer counting** | No filter | Include all initially |
| **Frequent k-mer selection** | ≥ 2 | Build graph from reliable k-mers |
| **Low-coverage node cleanup** | ≥ 2-3 | After graph construction |
| **Emergency cleanup** | ≥ 3 | Memory pressure |

### Contig Output Phase

| Quality Level | Min Coverage | Min Length | Use Case |
|--------------|-------------|------------|----------|
| **Strict (default)** | 2.0x | 3×k (63bp) | Standard assembly |
| **Conservative** | 3.0x | 5×k (105bp) | High-confidence only |
| **Sensitive** | 1.5x | 2×k (42bp) | Rare species detection |
| **Strain-aware** | 1.0x | 2×k (42bp) | Multi-strain preservation |

---

## Memory Budget Guidelines

### By System Configuration

| System | RAM | Budget for Assembly | Max K | Chunk Size |
|--------|-----|-------------------|-------|------------|
| **Low-end laptop** | 4 GB | 1 GB | 21 | 500 |
| **Mid-range laptop** | 8 GB | 2 GB | 31 | 1000 |
| **High-end laptop** | 16+ GB | 4 GB | 63 | 2000 |
| **Workstation** | 32+ GB | 8 GB | 127 | 5000 |
| **Server** | 64+ GB | 16 GB | 127 | 10000 |

### Memory Per K-mer

| Data Structure | Bytes per K-mer | Notes |
|----------------|----------------|-------|
| **Hash + count** | ~12 bytes | Minimal counter |
| **Graph node** | ~32 bytes | With coverage info |
| **Graph edge** | ~24 bytes | With weight |
| **String representation** | k+8 bytes | Overhead included |
| **2-bit encoding** | (k+3)/4 bytes | Compact storage |

---

## Assembly Quality Metrics

### Minimum Standards

| Metric | Minimum | Good | Excellent | Source |
|--------|---------|------|-----------|--------|
| **N50** | >1 kb | >10 kb | >100 kb | Community standard |
| **Longest contig** | >5 kb | >50 kb | >1 Mb | Assembly quality |
| **Total assembly size** | Match expected | ±10% expected | ±5% expected | Completeness |
| **Number of contigs** | <1000 | <100 | <50 | Fragmentation |
| **Average coverage** | >2x | >10x | >30x | Depth |

### CheckM Standards (MAGs)

| Quality Tier | Completeness | Contamination | Strain Het. |
|-------------|-------------|---------------|-------------|
| **High quality** | >90% | <5% | <10% |
| **Medium quality** | >50% | <10% | <25% |
| **Low quality** | >50% | >10% | >25% |
| **Fragment** | <50% | any | any |

---

## Filtering Thresholds by Downstream Use

### For Binning (Metagenomic Binning Tools)

| Tool | Min Contig Length | Recommended | Rationale |
|------|------------------|-------------|-----------|
| **MetaBAT** | 2,500 bp | 1,500 bp | Tetranucleotide frequency |
| **MaxBin** | 1,000 bp | 1,000 bp | Marker genes |
| **CONCOCT** | 1,000 bp | 1,000 bp | Coverage patterns |
| **General** | 1,000 bp | 1,500 bp | Reliable binning |

### For Gene Prediction

| Tool | Min Length | Recommended | Rationale |
|------|-----------|-------------|-----------|
| **Prodigal** | ~800 bp | 1,000 bp | Average gene length |
| **MetaGeneMark** | ~300 bp | 500 bp | Partial genes ok |
| **FragGeneScan** | ~100 bp | 500 bp | Fragment-aware |

### For Taxonomic Classification

| Tool | Min Length | Recommended | Rationale |
|------|-----------|-------------|-----------|
| **Kraken2** | 100 bp | 500 bp | K-mer based |
| **BLAST** | 200 bp | 1,000 bp | Homology search |
| **DIAMOND** | 200 bp | 1,000 bp | Fast alignment |
| **Marker genes** | 800 bp | 1,500 bp | Complete genes |

---

## Graph Simplification Strategies

### Tip Removal

| K-mer Size | Max Tip Length | Formula |
|-----------|---------------|---------|
| 15 | 30 bp | 2×k |
| 21 | 42 bp | 2×k |
| 31 | 62 bp | 2×k |
| 63 | 126 bp | 2×k |
| 127 | 254 bp | 2×k |

### Bulge Handling

| Bulge Type | Coverage Ratio | Action |
|-----------|---------------|--------|
| **Sequencing error** | <0.1 | Remove |
| **Rare strain variant** | 0.1-0.3 | Preserve or disconnect |
| **True biological variant** | >0.3 | Merge carefully |
| **Ambiguous** | Variable | Manual inspection |

### Edge Classification

| Edge Type | Criteria | Action |
|-----------|----------|--------|
| **Weak edge** | Coverage ratio > 10 | Disconnect (preserve) |
| **Filigree edge** | Low coverage + short | Remove |
| **Long edge** | >1500 bp | Likely correct |
| **High-confidence** | Coverage >10x | Keep |

---

## Validation Thresholds

### Coverage Uniformity

| Coefficient of Variation (CV) | Assessment | Action |
|------------------------------|------------|--------|
| **<0.2** | Excellent | Normal contig |
| **0.2-0.5** | Acceptable | Monitor |
| **0.5-1.0** | Suspicious | Check for chimera |
| **>1.0** | Problematic | Likely misassembly |

### GC Content Consistency

| GC Deviation | Assessment | Likely Cause |
|-------------|------------|--------------|
| **<5%** | Normal | Single species region |
| **5-10%** | Acceptable | Genomic variation |
| **10-20%** | Suspicious | Possible chimera |
| **>20%** | Problematic | Likely interspecies join |

### Biological Constraints

| Constraint | Threshold | Rationale |
|-----------|-----------|-----------|
| **Max contigs** | ≤ Number of reads | One contig per read maximum |
| **Max contig length** | ≤ Sum of read lengths | Cannot exceed input |
| **Min path length** | ≥ 3 k-mers | Biologically meaningful |
| **Coverage range** | 1x - 1000x | Biological plausibility |

---

## Performance Targets (Laptop)

### Assembly Speed

| Dataset Size | Target Time | Acceptable Time | Configuration |
|-------------|------------|----------------|--------------|
| **1K reads** | <10 sec | <30 sec | Low memory |
| **10K reads** | <1 min | <5 min | Low memory |
| **100K reads** | <10 min | <30 min | Medium memory |
| **1M reads** | <2 hrs | <6 hrs | High memory |

### Memory Usage

| Phase | Target | Maximum | Fallback Strategy |
|-------|--------|---------|------------------|
| **K-mer counting** | 25% budget | 50% budget | Cleanup rare k-mers |
| **Graph building** | 50% budget | 80% budget | Chunk processing |
| **Contig generation** | 25% budget | 40% budget | Streaming output |
| **Overall** | <80% budget | <95% budget | Emergency cleanup |

---

## Common Scenarios

### Single Species, High Coverage

```
Expected Output:
- Contigs: 1-10 (excellent) or 10-50 (good)
- N50: >50 kb
- Total length: ~genome size
- Coverage: Uniform (CV < 0.3)

Parameters:
- Min coverage: 5x (higher threshold)
- Min length: 5×k
- K-mer size: 31-63 (higher for repeats)
```

### Multi-Species, Uneven Coverage

```
Expected Output:
- Contigs: 50-500 (depends on diversity)
- N50: Variable by species
- Coverage: Multimodal distribution

Parameters:
- Min coverage: 2x (lower threshold)
- Min length: 3×k
- K-mer size: 21-31
- Enable strain-aware mode
```

### Low Coverage (<10x)

```
Expected Output:
- Contigs: High fragmentation expected
- Many gaps due to coverage variation
- May need multiple samples for assembly

Parameters:
- Min coverage: 1.5x (permissive)
- Min length: 2×k
- K-mer size: 15-21 (smaller for gaps)
- Sensitive mode
```

---

## References & Standards

| Standard | Source | Year |
|----------|--------|------|
| MetaSPAdes parameters | Nurk et al., Genome Research | 2017 |
| CheckM thresholds | Parks et al., Genome Biology | 2015 |
| Coverage recommendations | Community consensus | 2020-2025 |
| Multi-k strategies | MEGAHIT, IDBA-UD papers | 2015-2017 |
| Strain-aware methods | STRONG, Strainberry papers | 2021 |

**Full Research Report:** `docs/metagenomic_assembly_research.md`

---

## Quick Troubleshooting

| Symptom | Likely Cause | Solution |
|---------|-------------|----------|
| **Too many contigs** | Over-fragmentation | Increase k-mer size, check coverage |
| **Too few contigs** | Over-merging | Decrease k-mer size, enable strain mode |
| **Chimeric contigs** | Misassembly | Check coverage uniformity, add validation |
| **Low N50** | Fragmentation | Use multi-k, improve read quality |
| **High memory usage** | Large k-mers or dataset | Reduce k, use chunking, cleanup |
| **Timeout** | Complex graph | Reduce k, increase chunk size, cleanup |

---

**Quick Start:**
1. Use default parameters (k=21, min_cov=2x, min_len=63bp)
2. Check: Are contigs ≤ reads? (Biological sanity)
3. Calculate N50 and coverage uniformity
4. Adjust parameters based on results
5. For multi-strain: Enable strain-aware mode

**When in Doubt:** Conservative parameters (higher thresholds) prevent false positives
