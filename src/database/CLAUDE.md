# Database Module Guidelines

## Purpose
SQLite-based storage for k-mer indices, taxonomy data, and analysis results.

## Key Components
- `integration.rs` - Database connection and query interfaces

## Development Rules
- Always use prepared statements to prevent SQL injection
- Handle connection errors gracefully with retries
- Use transactions for batch operations
- Test with both small and large datasets

## Database Design Principles
- Normalize taxonomy data to reduce storage
- Index k-mer tables for fast lookups
- Use appropriate SQLite pragmas for performance
- Handle concurrent access safely

## Performance Optimization
- Batch inserts for k-mer data
- Use memory-mapped files for read-only operations
- Consider VACUUM operations for maintenance
- Profile query performance with EXPLAIN QUERY PLAN

## Error Handling
- Wrap SQLite errors in domain-specific types
- Provide clear error messages for database issues
- Handle schema migration failures gracefully
- Log database operations for debugging

## Testing Approach
- Test with real-world database sizes
- Verify data integrity after operations
- Test concurrent access scenarios
- Benchmark query performance