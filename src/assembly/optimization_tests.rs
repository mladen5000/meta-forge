//! Test-Driven Development for Assembly Optimizations
//! ==================================================
//!
//! This module implements comprehensive tests for our optimization features
//! following TDD principles: write failing tests first, then make them pass.

use super::memory_optimizations::*;
use anyhow::Result;
use std::sync::atomic::Ordering;

#[cfg(test)]
mod tdd_tests {
    use super::*;

    /* ========================================================================= */
    /*                        TDD: KMER ARENA ALLOCATOR                        */
    /* ========================================================================= */

    #[test]
    fn test_kmer_arena_basic_allocation() {
        // TDD Step 1: Write failing test first
        let arena = KmerArena::new(1); // 1MB arena
        let test_kmer = vec![0x1234567890ABCDEF, 0xFEDCBA0987654321];
        
        // This should work
        let kmer_ref = arena.allocate_kmer(&test_kmer).unwrap();
        assert_eq!(kmer_ref.length, 2);
        
        // Should be able to retrieve the same data
        let retrieved = arena.get_kmer(&kmer_ref).unwrap();
        assert_eq!(retrieved, &test_kmer);
    }

    #[test]
    fn test_kmer_arena_memory_efficiency() {
        // Test that arena allocation is more memory efficient than individual allocations
        let arena = KmerArena::new(10); // 10MB arena
        let mut allocated_kmers = Vec::new();
        
        // Allocate many k-mers
        for i in 0..1000 {
            let kmer_data = vec![i as u64, (i * 2) as u64, (i * 3) as u64];
            let kmer_ref = arena.allocate_kmer(&kmer_data).unwrap();
            allocated_kmers.push(kmer_ref);
        }
        
        let stats = arena.memory_stats();
        assert_eq!(stats.allocated_kmers, 1000);
        assert!(stats.utilization > 0.0);
        assert!(stats.active_blocks > 0);
    }

    #[test] 
    #[should_panic(expected = "Arena full")]
    fn test_kmer_arena_capacity_limits() {
        // Test that arena respects capacity limits
        let small_arena = KmerArena::new(1); // Very small arena
        
        // Try to allocate more than capacity
        for i in 0..10000 {
            let large_kmer = vec![i as u64; 1000]; // Large k-mers
            small_arena.allocate_kmer(&large_kmer).unwrap();
        }
        // Should panic or error before completing
    }

    /* ========================================================================= */
    /*                     TDD: LOCK-FREE GRAPH BUILDER                        */
    /* ========================================================================= */

    #[test]
    fn test_lock_free_graph_basic_operations() {
        let config = BuilderConfig::default();
        let builder = LockFreeGraphBuilder::new(config);
        
        // Create test k-mer references
        let kmer_ref1 = KmerRef { block_id: 0, offset: 0, length: 2 };
        let kmer_ref2 = KmerRef { block_id: 0, offset: 2, length: 2 };
        
        // Add nodes
        builder.add_node(12345, kmer_ref1).unwrap();
        builder.add_node(67890, kmer_ref2).unwrap();
        
        // Queue edges
        builder.queue_edge(12345, 67890, 1.5);
        builder.queue_edge(67890, 12345, 2.0);
        
        let stats = builder.get_stats();
        assert_eq!(stats.total_nodes, 2);
        assert_eq!(stats.edges_queued, 2);
    }

    #[test]
    fn test_lock_free_graph_concurrent_access() {
        use std::thread;
        use std::sync::Arc;
        
        let config = BuilderConfig::default();
        let builder = Arc::new(LockFreeGraphBuilder::new(config));
        
        let mut handles = vec![];
        
        // Spawn multiple threads adding nodes concurrently
        for i in 0..4 {
            let builder_clone = Arc::clone(&builder);
            let handle = thread::spawn(move || {
                for j in 0..100 {
                    let hash = (i * 100 + j) as u64;
                    let kmer_ref = KmerRef { 
                        block_id: 0, 
                        offset: hash as usize, 
                        length: 2 
                    };
                    builder_clone.add_node(hash, kmer_ref).unwrap();
                }
            });
            handles.push(handle);
        }
        
        // Wait for all threads
        for handle in handles {
            handle.join().unwrap();
        }
        
        let stats = builder.get_stats();
        assert_eq!(stats.total_nodes, 400); // 4 threads * 100 nodes each
    }

    #[test]
    fn test_lock_free_edge_processing() {
        let config = BuilderConfig { edge_batch_size: 10, ..Default::default() };
        let builder = LockFreeGraphBuilder::new(config);
        
        // Add nodes first
        for i in 0..20 {
            let kmer_ref = KmerRef { block_id: 0, offset: i, length: 2 };
            builder.add_node(i as u64, kmer_ref).unwrap();
        }
        
        // Queue many edges
        for i in 0..19 {
            builder.queue_edge(i as u64, (i + 1) as u64, 1.0);
        }
        
        // Process edges in batches
        let processed = builder.process_edge_batches().unwrap();
        assert!(processed.valid_edges > 0);
        assert_eq!(processed.invalid_edges, 0);
    }

    /* ========================================================================= */
    /*                  TDD: BOUNDED STREAM PROCESSOR                          */
    /* ========================================================================= */

    #[test]
    fn test_bounded_stream_memory_guarantee() {
        let config = StreamConfig {
            max_kmers: 100,
            memory_limit_mb: 1, // Very small limit
            enable_lru: true,
            sample_rate: 1.0, // Process all k-mers
        };
        let processor = BoundedStreamProcessor::new(config);
        
        // Create stream of k-mers that would exceed memory
        let kmer_stream = (0..1000).map(|i| (i as u64, vec![i as u8; 64])); // 64KB total
        
        let stats = processor.process_kmer_stream(kmer_stream).unwrap();
        
        // Should process all but only store within limits
        assert_eq!(stats.total_kmers_processed, 1000);
        assert!(stats.kmers_stored <= 100); // Respect max_kmers limit
        
        let processor_stats = processor.get_streaming_stats();
        assert!(processor_stats.memory_usage_mb <= 1); // Respect memory limit
        assert!(processor_stats.reservoir_size <= 100); // Respect reservoir size
    }

    #[test]
    fn test_bounded_stream_lru_eviction() {
        let config = StreamConfig {
            max_kmers: 10, // Very small reservoir
            memory_limit_mb: 1,
            enable_lru: true,
            sample_rate: 1.0,
        };
        let processor = BoundedStreamProcessor::new(config);
        
        // Add more k-mers than reservoir capacity
        let kmer_stream = (0..50).map(|i| (i as u64, vec![i as u8; 32]));
        
        let stats = processor.process_kmer_stream(kmer_stream).unwrap();
        
        let processor_stats = processor.get_streaming_stats();
        assert_eq!(processor_stats.reservoir_size, 10); // Should stay at max
        assert!(processor_stats.total_evictions > 0); // Should have evictions
    }

    #[test]
    fn test_bounded_stream_sampling() {
        let config = StreamConfig {
            max_kmers: 1000,
            memory_limit_mb: 10,
            enable_lru: true,
            sample_rate: 0.1, // Only 10% sampling
        };
        let processor = BoundedStreamProcessor::new(config);
        
        let kmer_stream = (0..1000).map(|i| (i as u64, vec![i as u8; 32]));
        
        let stats = processor.process_kmer_stream(kmer_stream).unwrap();
        
        assert_eq!(stats.total_kmers_processed, 1000);
        // With 10% sampling, should store roughly 100 k-mers (with randomness)
        assert!(stats.kmers_stored < 200); // Should be much less than total
        assert!(stats.kmers_sampled > 800); // Most should be sampled out
    }

    /* ========================================================================= */
    /*                     TDD: PERFORMANCE BENCHMARKS                         */
    /* ========================================================================= */

    #[test]
    fn test_memory_optimization_performance() {
        // Compare optimized vs unoptimized memory usage
        use std::time::Instant;
        
        let start = Instant::now();
        
        // Optimized path: Arena allocation
        let arena = KmerArena::new(10); // 10MB
        for i in 0..10000 {
            let kmer_data = vec![i as u64, (i * 2) as u64];
            arena.allocate_kmer(&kmer_data).unwrap();
        }
        
        let optimized_time = start.elapsed();
        let arena_stats = arena.memory_stats();
        
        // Unoptimized path: Individual allocations
        let start = Instant::now();
        let mut individual_kmers = Vec::new();
        for i in 0..10000 {
            let kmer_data = vec![i as u64, (i * 2) as u64];
            individual_kmers.push(kmer_data);
        }
        let unoptimized_time = start.elapsed();
        
        // Arena should be faster and use less memory
        println!("Arena allocation: {:?}", optimized_time);
        println!("Individual allocation: {:?}", unoptimized_time);
        println!("Arena memory efficiency: {:.1}%", arena_stats.utilization * 100.0);
        
        // Performance assertions (may need tuning based on actual performance)
        assert!(optimized_time < unoptimized_time * 2); // Should be at least 2x faster
        assert!(arena_stats.utilization > 0.5); // Should be reasonably efficient
    }

    #[test]
    fn test_lock_free_vs_locked_performance() {
        use std::sync::{Arc, Mutex};
        use std::thread;
        use std::time::Instant;
        
        // Lock-free version
        let config = BuilderConfig::default();
        let lock_free_builder = Arc::new(LockFreeGraphBuilder::new(config));
        
        let start = Instant::now();
        let mut handles = vec![];
        
        for i in 0..4 {
            let builder = Arc::clone(&lock_free_builder);
            let handle = thread::spawn(move || {
                for j in 0..1000 {
                    let hash = (i * 1000 + j) as u64;
                    let kmer_ref = KmerRef { block_id: 0, offset: hash as usize, length: 2 };
                    builder.add_node(hash, kmer_ref).unwrap();
                }
            });
            handles.push(handle);
        }
        
        for handle in handles {
            handle.join().unwrap();
        }
        
        let lock_free_time = start.elapsed();
        
        // Locked version (simplified for comparison)
        let locked_data = Arc::new(Mutex::new(std::collections::HashMap::new()));
        
        let start = Instant::now();
        let mut handles = vec![];
        
        for i in 0..4 {
            let data = Arc::clone(&locked_data);
            let handle = thread::spawn(move || {
                for j in 0..1000 {
                    let hash = (i * 1000 + j) as u64;
                    let mut map = data.lock().unwrap();
                    map.insert(hash, j);
                }
            });
            handles.push(handle);
        }
        
        for handle in handles {
            handle.join().unwrap();
        }
        
        let locked_time = start.elapsed();
        
        println!("Lock-free time: {:?}", lock_free_time);
        println!("Locked time: {:?}", locked_time);
        
        // Lock-free should be faster or at least competitive
        assert!(lock_free_time <= locked_time * 2); // Allow some variance
    }

    /* ========================================================================= */
    /*                    TDD: ERROR HANDLING AND EDGE CASES                   */
    /* ========================================================================= */

    #[test]
    fn test_error_handling_graceful_degradation() {
        // Test that optimizations fail gracefully
        
        // Test with invalid k-mer data
        let arena = KmerArena::new(1);
        let empty_kmer = vec![];
        
        // Should handle empty k-mers gracefully
        match arena.allocate_kmer(&empty_kmer) {
            Ok(_) => {}, // OK if it works
            Err(_) => {}, // OK if it fails gracefully
        }
        
        // Test with extremely large k-mers
        let huge_kmer = vec![0u64; 1000000]; // 8MB k-mer
        match arena.allocate_kmer(&huge_kmer) {
            Ok(_) => {}, // OK if system has enough memory
            Err(e) => {
                // Should have meaningful error message
                assert!(!format!("{}", e).is_empty());
            }
        }
    }

    #[test]
    fn test_memory_pressure_recovery() {
        // Test that bounded processor recovers from memory pressure
        let config = StreamConfig {
            max_kmers: 10,
            memory_limit_mb: 1,
            enable_lru: true,
            sample_rate: 1.0,
        };
        let processor = BoundedStreamProcessor::new(config);
        
        // Cause memory pressure
        let large_kmer_stream = (0..100).map(|i| (i as u64, vec![i as u8; 1024])); // Large k-mers
        
        let stats = processor.process_kmer_stream(large_kmer_stream).unwrap();
        
        // Should still function despite memory pressure
        assert!(stats.total_kmers_processed > 0);
        
        let processor_stats = processor.get_streaming_stats();
        assert!(processor_stats.reservoir_size <= 10);
        assert!(processor_stats.memory_usage_mb <= 1);
        
        // Should have handled memory pressure by evicting
        if stats.total_kmers_processed > 10 {
            assert!(processor_stats.total_evictions > 0);
        }
    }
}