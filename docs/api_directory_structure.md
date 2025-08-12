# REST API Directory Structure and Organization

## Complete Directory Structure

```
├── benches
│   └──  performance_profiler.rs
├── config
│   └──  default.toml
├── data
│   ├──  metagenomics.db
│   ├──  metagenomics.db-shm
│   └──  metagenomics.db-wal
├── docs
│   ├──  additional_test_coverage_recommendations.md
│   ├──  agent_orchestration_strategy.md
│   ├──  api_directory_structure.md
│   ├──  assembly_performance_analysis_report.md
│   ├──  bioinformatics_test_fixes.md
│   ├──  bioinformatics_test_validation_report.md
│   ├──  CLAUDE.md
│   ├──  metagenomic_assembly_algorithm_analysis.md
│   ├──  optimization_implementation_examples.md
│   ├── 󰂺 README.md
│   ├──  rest_api_architecture.md
│   ├──  tdd_remediation_plan.md
│   ├──  tui_architecture.md
│   └──  tui_progress_integration_analysis.md
├── examples
│   └──  test_kmer_ambiguous.rs
├── logs
│   ├──  pipeline.log.2025-08-08
│   ├──  pipeline.log.2025-08-09
│   ├──  pipeline.log.2025-08-10
│   └──  pipeline.log.2025-08-11
├── output
│   ├──  Sample_001_report.html
│   ├──  Sample_001_report.json
│   ├──  Sample_001_summary.tsv
│   ├──  short_SRR390728_1_report.html
│   ├──  short_SRR390728_1_report.json
│   ├──  short_SRR390728_1_summary.tsv
│   ├──  SRR390728_1_report.html
│   ├──  SRR390728_1_report.json
│   └──  SRR390728_1_summary.tsv
├── reads
│   ├──  2017.12.04_18.45.54_sample_0
│   ├──  output
│   ├──  SRR390728
│   ├──  tmp
│   ├──  work
│   ├──  index.html
│   ├──  sample_0.tar.gz
│   ├──  short_SRR390728_1.fastq
│   ├──  short_SRR390728_2.fastq
│   ├──  SRR390728_1.fastq
│   └──  SRR390728_2.fastq
├── scripts
│   ├──  batch_check.sh
│   ├──  cc.sh
│   ├──  cc_diff.sh
│   ├──  cc_fix.sh
│   ├──  coordinate_tdd_fixes.sh
│   ├──  diagnostic_script.sh
│   ├──  direct_cargo_check.sh
│   ├──  rs_preflight.sh
│   └──  run_checks.sh
├── src
│   ├──  assembly
│   ├──  bin
│   ├──  core
│   ├──  database
│   ├──  features
│   ├──  ml
│   ├──  pipeline
│   ├──  tui
│   ├──  utils
│   ├──  lib.rs
│   └──  main.rs
├── target
│   ├──  debug
│   ├──  tmp
│   └──  CACHEDIR.TAG
├── tests
│   ├──  tests_db
│   ├──  assembly_performance_test.rs
│   ├──  comprehensive_test_suite.rs
│   ├──  mod.rs
│   ├──  test_ambiguous_bases.rs
│   ├──  test_assembly_graph.rs
│   ├──  test_assembly_pipeline.rs
│   ├──  test_core_data_structures.rs
│   ├──  test_edge_cases_bioinformatics.rs
│   ├──  test_graph_algorithms.rs
│   ├──  test_integration_end_to_end.rs
│   ├──  test_property_based_genomics.rs
│   ├──  test_statistical_accuracy.rs
│   └──  test_synthetic_biological_data.rs
├── tmp
├── work
├── Cargo.lock
├──  Cargo.toml
├── check_compilation.sh
├── CLAUDE.md
├── imessage_recover.sh
├── kraken2_benchmark.sh
├── Makefile
├── names.dmp
├── OPTIMIZATION_REPORT.md
├── quick_check.sh
├── README.md
├── run_analysis.sh
├── run_diagnostics.sh
├── rust_toolchain.toml
├── simple_check.sh
└── TESTING_FIXES_ANALYSIS.md
```

## File Organization Principles

### 1. Layer Separation
Each major architectural layer has its own directory with clear responsibilities:
- **API Layer**: HTTP handling, routing, middleware
- **Service Layer**: Business logic and orchestration  
- **Data Layer**: Database access and data models
- **Job Layer**: Background processing and scheduling

### 2. Domain-Driven Structure
Within each layer, files are organized by domain (analysis, samples, users, etc.) rather than by technical concerns.

### 3. Trait-Based Testing
Service traits in `src/services/traits/` enable easy mocking and testing by providing interfaces that can be implemented by both real services and test doubles.

### 4. Configuration Management
Environment-specific configuration files allow easy deployment across different environments while maintaining security through environment variables for secrets.

### 5. Documentation Co-location
API documentation lives alongside code to ensure it stays current, with examples and guides organized by audience (developers, operators, users).

## Key Design Decisions

### Module Organization
- **Flat hierarchy within domains**: Avoid deep nesting that makes navigation difficult
- **Clear module boundaries**: Each module has well-defined responsibilities
- **Consistent naming**: Use domain-specific naming consistently across layers

### Testing Strategy
- **Unit tests**: Co-located with source code using `#[cfg(test)]` modules
- **Integration tests**: Separate `tests/` directory for cross-service testing
- **Test utilities**: Shared fixtures and utilities in `tests/common/`

### Configuration Approach
- **Layered configuration**: Default -> Environment -> CLI arguments
- **Type-safe configuration**: Rust structs with serde for configuration parsing
- **Environment-specific overrides**: Different files for dev/staging/production

### Error Handling
- **Structured errors**: Different error types for different layers
- **Error context**: Rich error information for debugging
- **API error mapping**: Clean conversion from internal errors to HTTP responses

This structure provides a solid foundation for the REST API while maintaining clear separation of concerns and enabling easy testing and maintenance.