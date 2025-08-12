# Features Module Guidelines

## Purpose
Feature extraction from genomic sequences and assembly graphs.

## Key Components
- `extraction.rs` - Sequence and graph feature computation

## Development Rules
- Ensure feature calculations are numerically stable
- Document the biological meaning of each feature
- Use appropriate data types for precision requirements
- Validate feature ranges and handle edge cases

## Feature Categories
- Sequence composition (GC content, k-mer frequencies)
- Graph topology metrics (connectivity, centrality)
- Assembly quality metrics (N50, coverage distribution)
- Taxonomic classification features

## Performance Guidelines
- Use iterators for memory efficiency
- Consider parallel computation for independent features
- Cache expensive calculations when appropriate
- Profile feature extraction on large datasets

## Numerical Considerations
- Use appropriate floating-point precision
- Handle zero-division and overflow cases
- Consider numerical stability for ratios and logarithms
- Document precision requirements for each feature

## Testing Strategy
- Unit tests with known expected values
- Test edge cases (empty sequences, single nucleotides)
- Compare against reference implementations
- Validate feature ranges and distributions