# MetaForge Pipeline Refactoring Plan

**Date**: October 10, 2025
**Status**: Phase 1 Complete ✅

## Overview

Refactoring the monolithic 4,319-line `complete_integration.rs` file into modular, maintainable components.

**Goal**: Reduce pipeline orchestrator to ~450 lines by extracting domain-specific logic to appropriate modules.

---

## Progress

### ✅ Phase 1: Core Pipeline Types (COMPLETED)

**Created**: `src/core/pipeline_types.rs` (170 lines)

**Extracted data structures:**
- `AbundanceProfile` - K-mer abundance profile
- `AnalysisResults` - Complete pipeline results
- `AssemblyResults` - Assembly phase results
- `TaxonomicClassification` - Classification data
- `FeatureCollection` - Extracted features
- `AnalysisReport` - Final report structure
- `PerformanceMetrics` - Performance data
- `ReportSummary` - Summary statistics
- `QualityMetrics` - Quality assessment
- `FileFormat` - File type detection

**Changes:**
- Created new module at [src/core/pipeline_types.rs](src/core/pipeline_types.rs)
- Updated [src/core/mod.rs](src/core/mod.rs) to export new types
- All types are now shared across pipeline coordinators

**Build Status**: ✅ Compiles successfully

---

## Remaining Phases

### 📋 Phase 2: QC File Reader Module

**Target file**: `src/qc/file_reader.rs` (~600 lines)

**Functions to extract from complete_integration.rs:**
- `preprocess_inputs()` (lines 1724-1894)
- `process_fastq_file()` (lines 1896-2088)
- `process_fasta_file()` (lines 2090-2120)
- `preprocess_inputs_with_progress()` (lines 2790-2885)
- `process_fastq_file_with_progress()` (lines 2887-3169)
- `detect_file_format()` (lines 2749-2786)

**New structure:**
```rust
pub struct FileReader {
    qc_config: QCPipelineConfig,
}

impl FileReader {
    pub fn new(config: QCPipelineConfig) -> Result<Self>;
    pub async fn preprocess_files(&self, inputs: &[PathBuf]) -> Result<(Vec<CorrectedRead>, QCStats)>;
    pub async fn preprocess_with_progress(&self, inputs: &[PathBuf], progress: &mut MultiProgress, line: usize) -> Result<(Vec<CorrectedRead>, QCStats)>;
    fn process_fastq(&self, path: &Path) -> Result<(Vec<CorrectedRead>, QCStats)>;
    fn process_fasta(&self, path: &Path) -> Result<Vec<CorrectedRead>>;
}
```

---

### 🧬 Phase 3: Assembly Coordinator Module

**Target file**: `src/assembly/coordinator.rs` (~300 lines)

**Functions to extract:**
- `run_assembly()` (lines 2122-2196)
- `run_assembly_with_progress()` (lines 3171-3258)
- `calculate_n50()` (lines 3574-3595)
- `calculate_n90()` (lines 3628-3649)
- `calculate_average_gc()` (lines 3597-3626)
- `calculate_coverage_std()` (lines 3651-3665)

**New structure:**
```rust
pub struct AssemblyCoordinator {
    config: AssemblyConfig,
    laptop_config: LaptopConfig,
}

impl AssemblyCoordinator {
    pub fn new(config: AssemblyConfig) -> Result<Self>;
    pub async fn assemble(&self, reads: &[CorrectedRead]) -> Result<AssemblyResults>;
    pub async fn assemble_with_progress(&self, reads: &[CorrectedRead], progress: &mut MultiProgress, line: usize) -> Result<AssemblyResults>;
    fn calculate_stats(&self, contigs: &[Contig]) -> AssemblyStats;
}
```

---

### 🔍 Phase 4: Feature Coordinator Module

**Target file**: `src/features/coordinator.rs` (~200 lines)

**Functions to extract:**
- `extract_features()` (lines 2198-2239)
- `extract_features_with_progress()` (lines 3260-3365)

**New structure:**
```rust
pub struct FeatureCoordinator {
    config: FeatureConfig,
    extractor: Arc<AdvancedFeatureExtractor>,
}

impl FeatureCoordinator {
    pub fn new(config: FeatureConfig) -> Result<Self>;
    pub async fn extract_all(&self, assembly: &AssemblyResults) -> Result<FeatureCollection>;
    pub async fn extract_with_progress(&self, assembly: &AssemblyResults, progress: &mut MultiProgress, line: usize) -> Result<FeatureCollection>;
}
```

---

### 🏷️ Phase 5: Classification Coordinator Module

**Target file**: `src/ml/coordinator.rs` (~300 lines)

**Functions to extract:**
- `classify_sequences()` (lines 2241-2314)
- `classify_sequences_with_progress()` (lines 3367-3444)
- `bin_contigs_by_coverage()` (lines 2316-2369)

**New structure:**
```rust
pub struct ClassificationCoordinator {
    config: ClassificationConfig,
    classifier: SimpleContigClassifier,
}

impl ClassificationCoordinator {
    pub fn new(config: ClassificationConfig) -> Result<Self>;
    pub async fn classify(&self, assembly: &AssemblyResults, features: &FeatureCollection) -> Result<Vec<TaxonomicClassification>>;
    pub async fn classify_with_progress(&self, assembly: &AssemblyResults, features: &FeatureCollection, progress: &mut MultiProgress, line: usize) -> Result<Vec<TaxonomicClassification>>;
}
```

---

### 📊 Phase 6: Abundance Estimator Module

**Target file**: `src/ml/abundance_estimator.rs` (~200 lines)

**Functions to extract:**
- `estimate_abundance()` (lines 2371-2385)
- `estimate_abundance_with_progress()` (lines 3446-3504)
- `hash_kmer()` (lines 3564-3572)

**New structure:**
```rust
pub struct AbundanceEstimator {
    kmer_size: usize,
}

impl AbundanceEstimator {
    pub fn new(kmer_size: usize) -> Self;
    pub async fn estimate(&self, reads: &[CorrectedRead]) -> Result<AbundanceProfile>;
    pub async fn estimate_with_progress(&self, reads: &[CorrectedRead], progress: &mut MultiProgress, line: usize) -> Result<AbundanceProfile>;
    fn hash_kmer(&self, kmer: &[u8]) -> u64;
}
```

---

### 📝 Phase 7: Report Coordinator Module

**Target file**: `src/reporting/report_coordinator.rs` (~500 lines)

**Functions to extract:**
- `generate_report()` (lines 2387-2439)
- `generate_report_with_progress()` (lines 3506-3562)
- `write_report_files()` (lines 2441-2547)
- `generate_kraken_reports()` (lines 2549-2591)
- `generate_html_report()` (lines 2593-2704)
- `generate_tsv_summary()` (lines 2706-2729)
- `calculate_assembly_completeness()` (lines 3667-3681)
- `calculate_coverage_uniformity()` (lines 3683-3700)

**New structure:**
```rust
pub struct ReportCoordinator {
    config: IOConfig,
    output_manager: IntermediateOutputManager,
}

impl ReportCoordinator {
    pub fn new(config: IOConfig, output_manager: IntermediateOutputManager) -> Result<Self>;
    pub async fn generate(&self, sample_name: &str, assembly: &AssemblyResults, classifications: &[TaxonomicClassification], abundance: &AbundanceProfile) -> Result<AnalysisReport>;
    pub async fn generate_with_progress(&self, sample_name: &str, assembly: &AssemblyResults, classifications: &[TaxonomicClassification], abundance: &AbundanceProfile, progress: &mut MultiProgress, line: usize) -> Result<AnalysisReport>;
    async fn write_all_formats(&self, report: &AnalysisReport) -> Result<()>;
}
```

---

### 🎯 Phase 8: Simplified Pipeline Orchestrator

**Target**: Reduce `complete_integration.rs` from 4,319 to ~450 lines

**Final structure:**
```rust
pub struct MetagenomicsPipeline {
    config: PipelineConfiguration,
    database: Option<MetagenomicsDatabase>,
    resource_monitor: ResourceMonitor,
    output_manager: IntermediateOutputManager,
    performance_profiler: PerformanceProfiler,

    // Coordinators (composition over inheritance)
    file_reader: FileReader,
    assembly_coordinator: AssemblyCoordinator,
    feature_coordinator: FeatureCoordinator,
    classification_coordinator: ClassificationCoordinator,
    abundance_estimator: AbundanceEstimator,
    report_coordinator: ReportCoordinator,
}

impl MetagenomicsPipeline {
    // High-level orchestration only
    pub async fn run_analysis(&mut self, inputs: &[PathBuf], sample_name: &str, mode: AnalysisMode) -> Result<AnalysisResults> {
        // Delegate to coordinators
        let (reads, qc_stats) = self.file_reader.preprocess_files(inputs).await?;
        let assembly = self.assembly_coordinator.assemble(&reads).await?;
        let features = self.feature_coordinator.extract_all(&assembly).await?;
        let classifications = self.classification_coordinator.classify(&assembly, &features).await?;
        let abundance = self.abundance_estimator.estimate(&reads).await?;
        let report = self.report_coordinator.generate(sample_name, &assembly, &classifications, &abundance).await?;

        Ok(AnalysisResults { /* ... */ })
    }
}
```

---

## Benefits of Refactoring

### ✅ Modularity
- Each phase (QC, assembly, classification, etc.) is self-contained
- Clear separation of concerns
- Easy to test individual components

### ✅ Maintainability
- Changes to QC don't affect classification
- Each coordinator has single responsibility
- Easier code review and debugging

### ✅ Reusability
- Coordinators can be used independently
- Other tools can import specific coordinators
- No need to use entire pipeline for one operation

### ✅ Testability
- Unit tests for each coordinator
- Integration tests for pipeline orchestration
- Mock coordinators for testing

### ✅ Code Organization
- Files under 500 lines each
- Logical grouping by domain
- Clear API boundaries

---

## File Structure (After Refactoring)

```
src/
├── core/
│   ├── data_structures.rs      (existing)
│   ├── paired_reads.rs          (existing)
│   └── pipeline_types.rs        (NEW - 170 lines) ✅
├── qc/
│   ├── quality_filter.rs        (existing)
│   ├── adapter_trimmer.rs       (existing)
│   ├── qc_stats.rs              (existing)
│   └── file_reader.rs           (NEW - ~600 lines)
├── assembly/
│   ├── laptop_assembly.rs       (existing)
│   ├── adaptive_k.rs            (existing)
│   ├── orchestrator.rs          (existing)
│   └── coordinator.rs           (NEW - ~300 lines)
├── features/
│   ├── extraction.rs            (existing)
│   └── coordinator.rs           (NEW - ~200 lines)
├── ml/
│   ├── simple_classifier.rs     (existing)
│   ├── kmer_taxonomy.rs         (existing)
│   ├── coordinator.rs           (NEW - ~300 lines)
│   └── abundance_estimator.rs   (NEW - ~200 lines)
├── reporting/
│   ├── report_generator.rs      (existing)
│   └── report_coordinator.rs    (NEW - ~500 lines)
└── pipeline/
    └── complete_integration.rs  (4319 → ~450 lines)
```

---

## Next Steps

### Immediate (This Session)
1. ✅ Create `core/pipeline_types.rs`
2. ✅ Update imports in `core/mod.rs`
3. ✅ Verify compilation
4. 📝 Create this refactoring plan document

### Near Term (Next Session)
1. Create `qc/file_reader.rs` - Extract preprocessing logic
2. Update `complete_integration.rs` to use `FileReader`
3. Create `assembly/coordinator.rs` - Extract assembly orchestration
4. Create `ml/abundance_estimator.rs` - Extract abundance logic

### Medium Term
1. Create `ml/coordinator.rs` - Extract classification logic
2. Create `features/coordinator.rs` - Extract feature extraction
3. Create `reporting/report_coordinator.rs` - Extract reporting logic

### Final Steps
1. Simplify `complete_integration.rs` to pure orchestration
2. Update all imports across the codebase
3. Run full test suite
4. Update documentation
5. Performance benchmarks (ensure no regression)

---

## Migration Strategy

### Incremental Approach
- Extract one module at a time
- Compile and test after each extraction
- Use git commits per phase for easy rollback
- Keep `complete_integration.rs` functional throughout

### Testing Strategy
- Existing tests should continue to pass
- Add unit tests for each new coordinator
- Integration tests for full pipeline

### Risk Mitigation
- Git branch for refactoring (`refactor/modular-pipeline`)
- Benchmark before/after each phase
- Keep old code commented until all tests pass

---

## Success Criteria

### Code Metrics
- ✅ Pipeline orchestrator < 500 lines
- ✅ No file > 600 lines
- ✅ All coordinators < 400 lines
- ✅ Clear module boundaries

### Functionality
- ✅ All existing tests pass
- ✅ No performance regression (< 5% slower acceptable)
- ✅ Same output files generated
- ✅ CLI interface unchanged

### Code Quality
- ✅ No circular dependencies
- ✅ Clear API documentation
- ✅ Consistent error handling
- ✅ Proper logging throughout

---

## Notes

### Dependencies
- Data structures in `core/pipeline_types.rs` prevent circular dependencies
- All coordinators depend only on `core` types, not each other
- Pipeline orchestrator depends on all coordinators

### Backward Compatibility
- Public API remains the same
- CLI unchanged
- Output format identical
- Configuration format preserved

---

**Last Updated**: October 10, 2025
**Author**: Claude Code
**Status**: Phase 1 Complete - Ready for Phase 2
