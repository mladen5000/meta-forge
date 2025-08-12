# Utils Module Guidelines

## Purpose
Shared utilities and helper functions for the metagenomic pipeline.

## Key Components
- `configuration.rs` - Configuration parsing and validation
- `genomic_validator.rs` - Sequence validation and quality checks
- `progress_display.rs` - Progress reporting utilities
- `streaming_abundance.rs` - Abundance calculation streaming

## Development Rules
- Write utilities to be reusable across modules
- Provide comprehensive error handling
- Document performance characteristics
- Keep utilities focused and single-purpose

## Design Principles
- Prefer composition over inheritance
- Use trait objects for extensible behavior
- Minimize dependencies between utilities
- Design for testability and mockability

## Configuration Handling
- Support multiple configuration formats (TOML, JSON, YAML)
- Validate configuration at startup
- Provide clear error messages for invalid config
- Support environment variable overrides

## Validation Functions
- Use appropriate data types for genomic data
- Handle edge cases in sequence validation
- Provide detailed validation error messages
- Consider performance for large-scale validation

## Progress Reporting
- Use consistent progress reporting interfaces
- Support both CLI and programmatic progress tracking
- Handle cancellation and interruption gracefully
- Provide ETA estimation when possible

## Testing Strategy
- Unit tests for each utility function
- Property-based tests for validation functions
- Mock external dependencies in tests
- Test error conditions and edge cases