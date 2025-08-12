# Assembly Module Guidelines

## Purpose
Metagenomic sequence assembly with adaptive k-mer strategies and graph-based approaches.

## Key Components
- `adaptive_k.rs` - Dynamic k-mer size selection
- `graph_construction.rs` - De Bruijn graph assembly
- `bioinformatics_optimizations.rs` - Performance optimizations
- `optimized_structures.rs` - Memory-efficient data structures

## Development Rules
- Follow TDD: write tests before implementation
- Use `cargo clippy` and `cargo fmt` before commits
- Prefer `Result` over panics for error handling
- Keep functions under 50 lines, files under 500 lines

## Performance Considerations
- Optimize for large genomic datasets
- Use SIMD when possible for k-mer operations
- Consider memory layout for cache efficiency
- Profile with `cargo bench` for critical paths

## Testing
- Unit tests for each algorithm
- Integration tests with real FASTQ data
- Benchmark comparisons with existing tools

## Architecture
- Modular design allows swapping algorithms
- Clear separation between data structures and algorithms
- Streaming processing for large files
- Error propagation via `Result<T, AssemblyError>`

## Dependencies
- `rayon` for parallelization
- `bio` crate for sequence handling
- Custom data structures in `core::data_structures`