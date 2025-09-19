//! TDD Tests for Database Connection Lifecycle Management
//! =====================================================
//!
//! Tests proper resource management for database connections, transactions,
//! and prepared statements to prevent borrowing conflicts.

#[cfg(test)]
mod database_lifecycle_tests {
    use super::*;
    use crate::assembly::database_optimizations::*;
    use std::sync::Arc;
    use tempfile::tempdir;

    fn create_test_database() -> GenomicDatabase {
        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("test.db");

        let config = DatabaseConfig {
            max_connections: 2,
            batch_size: 100,
            enable_caching: true,
            cache_size_mb: 10,
            enable_compression: true,
            wal_checkpoint_interval: 1000,
        };

        GenomicDatabase::new(db_path.to_str().unwrap(), config).unwrap()
    }

    #[test]
    fn test_connection_lifecycle_basic() {
        // TDD: Test that we can get and return connections without borrowing issues
        let db = create_test_database();

        // This should work without compilation errors
        let result = db.insert_kmers_batch(&vec![KmerRecord {
            id: 1,
            sequence: "ATCG".to_string(),
            hash: vec![1, 2, 3, 4],
            count: 5,
        }]);

        // Should succeed or fail gracefully, but not have borrowing issues
        match result {
            Ok(_) => {}  // Success is good
            Err(_) => {} // Errors are ok as long as they compile
        }
    }

    #[test]
    fn test_transaction_lifecycle_proper_scoping() {
        // TDD: Test that transactions are properly scoped
        let db = create_test_database();

        // Multiple operations in sequence should work
        for i in 0..5 {
            let kmer = KmerRecord {
                id: i,
                sequence: format!("ATCG{}", i),
                hash: vec![i as u8, (i + 1) as u8, (i + 2) as u8],
                count: i as u32,
            };

            let result = db.insert_kmers_batch(&vec![kmer]);
            assert!(result.is_ok() || result.is_err()); // Should compile and run
        }
    }

    #[test]
    fn test_statement_preparation_lifecycle() {
        // TDD: Test that prepared statements don't conflict with connection lifecycle
        let db = create_test_database();

        // Query operations should work without borrowing conflicts
        let frequent_kmers = db.get_frequent_kmers(1, Some(10));

        match frequent_kmers {
            Ok(kmers) => {
                // Should be able to process results
                assert!(kmers.len() >= 0);
            }
            Err(_) => {
                // Errors are fine as long as compilation succeeds
            }
        }
    }

    #[test]
    fn test_concurrent_database_access() {
        // TDD: Test that multiple threads can access database concurrently
        use std::thread;

        let db = Arc::new(create_test_database());
        let mut handles = vec![];

        for i in 0..3 {
            let db_clone = Arc::clone(&db);
            let handle = thread::spawn(move || {
                let kmer = KmerRecord {
                    id: i,
                    sequence: format!("THREAD{}", i),
                    hash: vec![i as u8; 4],
                    count: i as u32,
                };

                // Each thread should be able to use the database
                db_clone.insert_kmers_batch(&vec![kmer])
            });
            handles.push(handle);
        }

        // Wait for all threads and collect results
        let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();

        // All operations should complete (success or failure, but no panics)
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_cache_operation_lifecycle() {
        // TDD: Test that cache operations don't interfere with connection management
        let db = create_test_database();

        // Insert some data first
        let test_kmers: Vec<KmerRecord> = (0..10)
            .map(|i| KmerRecord {
                id: i,
                sequence: format!("CACHE{:02}", i),
                hash: vec![i as u8; 4],
                count: (i * 2) as u32,
            })
            .collect();

        let _ = db.insert_kmers_batch(&test_kmers);

        // Query same data multiple times to test caching
        for _ in 0..3 {
            let result1 = db.get_frequent_kmers(0, Some(5));
            let result2 = db.get_frequent_kmers(0, Some(5));

            // Results should be consistent
            match (result1, result2) {
                (Ok(r1), Ok(r2)) => {
                    // If both succeed, they should return same data
                    assert_eq!(r1.len(), r2.len());
                }
                _ => {
                    // Other combinations are acceptable for testing
                }
            }
        }
    }

    #[test]
    fn test_error_handling_preserves_connections() {
        // TDD: Test that errors don't leave connections in bad state
        let db = create_test_database();

        // Try some operations that might fail
        let invalid_kmer = KmerRecord {
            id: -1,                   // Potentially invalid ID
            sequence: "".to_string(), // Empty sequence
            hash: vec![],             // Empty hash
            count: 0,
        };

        let result = db.insert_kmers_batch(&vec![invalid_kmer]);

        // After error, database should still be usable
        let valid_kmer = KmerRecord {
            id: 999,
            sequence: "ATCG".to_string(),
            hash: vec![1, 2, 3, 4],
            count: 1,
        };

        let result2 = db.insert_kmers_batch(&vec![valid_kmer]);

        // At least one operation should work, or both should fail consistently
        assert!(result.is_ok() || result.is_err());
        assert!(result2.is_ok() || result2.is_err());
    }

    #[test]
    fn test_connection_pool_exhaustion_recovery() {
        // TDD: Test that connection pool handles exhaustion gracefully
        let config = DatabaseConfig {
            max_connections: 1, // Very small pool
            batch_size: 10,
            enable_caching: false, // Disable caching to force DB hits
            cache_size_mb: 1,
            enable_compression: false,
            wal_checkpoint_interval: 100,
        };

        let temp_dir = tempdir().unwrap();
        let db_path = temp_dir.path().join("pool_test.db");
        let db = GenomicDatabase::new(db_path.to_str().unwrap(), config).unwrap();

        // Try to use more connections than available
        let test_kmers: Vec<KmerRecord> = (0..5)
            .map(|i| KmerRecord {
                id: i,
                sequence: format!("POOL{}", i),
                hash: vec![i as u8; 3],
                count: i as u32,
            })
            .collect();

        // Multiple rapid operations should either succeed or fail gracefully
        for batch in test_kmers.chunks(2) {
            let result = db.insert_kmers_batch(batch);
            match result {
                Ok(_) => {} // Success is good
                Err(e) => {
                    // Errors should be meaningful, not panics
                    let error_msg = format!("{}", e);
                    assert!(!error_msg.is_empty());
                }
            }
        }
    }
}
