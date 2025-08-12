# Core Module Guidelines

## Purpose
Fundamental data structures and utilities for metagenomic analysis.

## Key Components
- `data_structures.rs` - Core genomic data types
- `paired_reads.rs` - Paired-end read handling

## Development Rules
- This is foundational code - extra care with API design
- All public APIs must be well-documented
- Prefer zero-copy operations where possible
- Use `#[inline]` for hot path functions

## Data Structure Principles
- Memory-efficient representations
- Cache-friendly layouts
- Generic over sequence types when reasonable
- Implement standard traits (Clone, Debug, etc.)

## Safety & Performance
- Use `unsafe` only when necessary with clear justification
- Profile memory usage with large datasets
- Consider SIMD operations for bulk operations
- Validate inputs at API boundaries

## Testing Strategy
- Property-based tests for core algorithms
- Memory leak detection in CI
- Fuzz testing for parsing functions
- Benchmark critical data structure operations

## API Stability
- This module forms the foundation - maintain compatibility
- Use deprecation warnings before breaking changes
- Document performance characteristics in rustdoc