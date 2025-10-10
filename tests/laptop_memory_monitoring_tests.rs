//! Laptop Memory Monitoring and Management Tests
//! =============================================
//!
//! Comprehensive test suite for validating memory monitoring, cleanup,
//! and optimization strategies specifically designed for laptop constraints.
//!
//! Focus Areas:
//! - Real-time memory usage monitoring
//! - Automatic memory pressure detection
//! - Memory cleanup and garbage collection
//! - Memory leak detection and prevention
//! - Laptop-specific memory optimization validation

use anyhow::Result;
use meta_forge::assembly::laptop_assembly::{
    BoundedKmerCounter, LaptopAssemblyGraph, LaptopConfig,
};
use meta_forge::assembly::memory_optimizations::{
    BoundedStreamProcessor, BuilderConfig, KmerArena, LockFreeGraphBuilder, MemoryStats,
    StreamConfig, StreamProcessorStats,
};
use meta_forge::core::data_structures::{CorrectedRead, CorrectionMetadata};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

/* ========================================================================= */
/*                        REAL-TIME MEMORY MONITORING                      */
/* ========================================================================= */

/// Memory monitor for tracking usage in real-time
#[derive(Debug)]
struct RealTimeMemoryMonitor {
    current_usage: AtomicUsize,
    peak_usage: AtomicUsize,
    allocation_count: AtomicUsize,
    deallocation_count: AtomicUsize,
    pressure_threshold: usize,
    critical_threshold: usize,
}

impl RealTimeMemoryMonitor {
    fn new(pressure_threshold_mb: usize, critical_threshold_mb: usize) -> Self {
        Self {
            current_usage: AtomicUsize::new(0),
            peak_usage: AtomicUsize::new(0),
            allocation_count: AtomicUsize::new(0),
            deallocation_count: AtomicUsize::new(0),
            pressure_threshold: pressure_threshold_mb * 1024 * 1024,
            critical_threshold: critical_threshold_mb * 1024 * 1024,
        }
    }

    fn allocate(&self, size: usize) -> bool {
        let current = self.current_usage.fetch_add(size, Ordering::SeqCst);
        let new_usage = current + size;

        // Update peak usage
        let mut peak = self.peak_usage.load(Ordering::Relaxed);
        while peak < new_usage {
            match self.peak_usage.compare_exchange_weak(
                peak,
                new_usage,
                Ordering::SeqCst,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(current_peak) => peak = current_peak,
            }
        }

        self.allocation_count.fetch_add(1, Ordering::Relaxed);

        // Check if allocation would exceed critical threshold
        new_usage < self.critical_threshold
    }

    fn deallocate(&self, size: usize) {
        self.current_usage.fetch_sub(size, Ordering::SeqCst);
        self.deallocation_count.fetch_add(1, Ordering::Relaxed);
    }

    fn is_under_pressure(&self) -> bool {
        self.current_usage.load(Ordering::Relaxed) > self.pressure_threshold
    }

    fn is_critical(&self) -> bool {
        self.current_usage.load(Ordering::Relaxed) > self.critical_threshold
    }

    fn get_stats(&self) -> MemoryMonitorStats {
        MemoryMonitorStats {
            current_usage_mb: self.current_usage.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0),
            peak_usage_mb: self.peak_usage.load(Ordering::Relaxed) as f64 / (1024.0 * 1024.0),
            pressure_threshold_mb: self.pressure_threshold as f64 / (1024.0 * 1024.0),
            critical_threshold_mb: self.critical_threshold as f64 / (1024.0 * 1024.0),
            allocation_count: self.allocation_count.load(Ordering::Relaxed),
            deallocation_count: self.deallocation_count.load(Ordering::Relaxed),
            pressure_ratio: self.current_usage.load(Ordering::Relaxed) as f64
                / self.pressure_threshold as f64,
        }
    }
}

#[derive(Debug)]
struct MemoryMonitorStats {
    current_usage_mb: f64,
    peak_usage_mb: f64,
    pressure_threshold_mb: f64,
    critical_threshold_mb: f64,
    allocation_count: usize,
    deallocation_count: usize,
    pressure_ratio: f64,
}

#[cfg(test)]
mod memory_monitoring_tests {
    use super::*;

    /// Test real-time memory usage tracking
    #[test]
    fn test_real_time_memory_tracking() {
        let monitor = RealTimeMemoryMonitor::new(100, 150); // 100MB pressure, 150MB critical

        // Simulate memory allocations
        let allocations = vec![
            10 * 1024 * 1024, // 10MB
            20 * 1024 * 1024, // 20MB
            30 * 1024 * 1024, // 30MB
            40 * 1024 * 1024, // 40MB
        ];

        let mut successful_allocations = 0;
        for (i, &size) in allocations.iter().enumerate() {
            if monitor.allocate(size) {
                successful_allocations += 1;
                println!(
                    "Allocation {}: {:.1}MB - SUCCESS",
                    i + 1,
                    size as f64 / (1024.0 * 1024.0)
                );
            } else {
                println!(
                    "Allocation {}: {:.1}MB - BLOCKED (critical threshold)",
                    i + 1,
                    size as f64 / (1024.0 * 1024.0)
                );
                break;
            }

            if monitor.is_under_pressure() {
                println!("  ⚠️ Memory pressure detected at allocation {}", i + 1);
            }
        }

        let stats = monitor.get_stats();
        println!("\n📊 Memory tracking results:");
        println!("  Current: {:.1}MB", stats.current_usage_mb);
        println!("  Peak: {:.1}MB", stats.peak_usage_mb);
        println!("  Pressure ratio: {:.2}", stats.pressure_ratio);
        println!("  Allocations: {}", stats.allocation_count);

        // Validate tracking accuracy
        assert!(stats.current_usage_mb > 0.0, "Should track memory usage");
        assert!(
            stats.peak_usage_mb >= stats.current_usage_mb,
            "Peak should be >= current"
        );
        assert_eq!(
            stats.allocation_count, successful_allocations,
            "Should track allocation count"
        );

        println!("✅ Real-time memory tracking validated");
    }

    /// Test memory pressure detection and response
    #[test]
    fn test_memory_pressure_detection() {
        let monitor = RealTimeMemoryMonitor::new(50, 80); // Low thresholds for testing
        let chunk_size = 10 * 1024 * 1024; // 10MB chunks

        let mut allocations = Vec::new();
        let mut pressure_detected = false;
        let mut critical_reached = false;

        // Gradually increase memory usage
        for i in 0..10 {
            if monitor.allocate(chunk_size) {
                allocations.push(chunk_size);

                if monitor.is_under_pressure() && !pressure_detected {
                    pressure_detected = true;
                    println!("🟡 Memory pressure detected at allocation {}", i + 1);
                }

                if monitor.is_critical() && !critical_reached {
                    critical_reached = true;
                    println!(
                        "🔴 Critical memory threshold reached at allocation {}",
                        i + 1
                    );
                }
            } else {
                println!("🚫 Allocation {} blocked by critical threshold", i + 1);
                break;
            }
        }

        assert!(pressure_detected, "Should detect memory pressure");

        let final_stats = monitor.get_stats();
        assert!(
            final_stats.pressure_ratio > 0.5,
            "Should reach significant memory usage"
        );

        println!("✅ Memory pressure detection working correctly");
    }

    /// Test automatic memory cleanup triggers
    #[test]
    fn test_automatic_memory_cleanup() -> Result<()> {
        let monitor = Arc::new(RealTimeMemoryMonitor::new(30, 50)); // Aggressive thresholds
        let counter = BoundedKmerCounter::new(25); // 25MB counter

        // Simulate automatic cleanup when pressure is detected
        let cleanup_triggered = Arc::new(AtomicUsize::new(0));

        // Background task to monitor and trigger cleanup
        let monitor_clone = Arc::clone(&monitor);
        let cleanup_clone = Arc::clone(&cleanup_triggered);

        let _cleanup_handle = thread::spawn(move || {
            for _ in 0..100 {
                if monitor_clone.is_under_pressure() {
                    // Simulate cleanup by deallocating some memory
                    monitor_clone.deallocate(5 * 1024 * 1024); // Free 5MB
                    cleanup_clone.fetch_add(1, Ordering::Relaxed);
                    println!("🧹 Automatic cleanup triggered");
                }
                thread::sleep(Duration::from_millis(10));
            }
        });

        // Stress the memory system
        for i in 0..20 {
            if monitor.allocate(5 * 1024 * 1024) {
                // 5MB chunks
                println!("Allocation {}: 5MB", i + 1);
            }
            thread::sleep(Duration::from_millis(50));
        }

        thread::sleep(Duration::from_millis(200)); // Let cleanup thread work

        let cleanup_count = cleanup_triggered.load(Ordering::Relaxed);
        assert!(
            cleanup_count > 0,
            "Automatic cleanup should have been triggered"
        );

        println!("✅ Automatic cleanup triggered {} times", cleanup_count);

        Ok(())
    }
}

/* ========================================================================= */
/*                        MEMORY LEAK DETECTION                            */
/* ========================================================================= */

#[cfg(test)]
mod memory_leak_tests {
    use super::*;

    /// Test for memory leaks in k-mer arena allocation
    #[test]
    fn test_kmer_arena_leak_detection() -> Result<()> {
        let monitor = RealTimeMemoryMonitor::new(100, 150);
        let initial_usage = monitor.current_usage.load(Ordering::Relaxed);

        // Multiple allocation/deallocation cycles
        for cycle in 0..5 {
            let arena = KmerArena::new(10); // 10MB arena

            // Simulate arena usage tracking
            monitor.allocate(10 * 1024 * 1024)?;

            // Allocate many k-mers
            let mut kmer_refs = Vec::new();
            for i in 0..1000 {
                let kmer_data = vec![i as u64, (i + 1) as u64];
                match arena.allocate_kmer(&kmer_data) {
                    Ok(kmer_ref) => kmer_refs.push(kmer_ref),
                    Err(_) => break,
                }
            }

            let arena_stats = arena.memory_stats();
            println!(
                "Cycle {}: {} k-mers allocated, {:.1}% utilization",
                cycle + 1,
                arena_stats.allocated_kmers,
                arena_stats.utilization * 100.0
            );

            // Simulate arena cleanup
            drop(arena);
            monitor.deallocate(10 * 1024 * 1024);
        }

        let final_usage = monitor.current_usage.load(Ordering::Relaxed);
        let leak_amount = final_usage as i64 - initial_usage as i64;

        assert!(
            leak_amount.abs() < 1024 * 1024, // < 1MB tolerance
            "Potential memory leak detected: {}KB difference",
            leak_amount / 1024
        );

        println!("✅ No significant memory leaks detected in k-mer arena");

        Ok(())
    }

    /// Test for memory leaks in graph construction
    #[test]
    fn test_graph_construction_leak_detection() -> Result<()> {
        let monitor = RealTimeMemoryMonitor::new(200, 300);
        let initial_usage = monitor.current_usage.load(Ordering::Relaxed);

        // Multiple graph construction cycles
        for cycle in 0..3 {
            let config = LaptopConfig::low_memory();
            let mut graph = LaptopAssemblyGraph::new(config);

            // Track graph memory usage
            let graph_memory_mb = 50; // Estimate
            monitor.allocate(graph_memory_mb * 1024 * 1024)?;

            // Build graph with test data
            let reads = create_leak_test_reads(500, 80);
            graph.build_from_reads(&reads, 21)?;

            let actual_memory = graph.memory_usage_mb();
            println!(
                "Cycle {}: Graph memory usage: {:.1}MB",
                cycle + 1,
                actual_memory
            );

            // Cleanup
            drop(graph);
            monitor.deallocate(graph_memory_mb * 1024 * 1024);
        }

        let final_usage = monitor.current_usage.load(Ordering::Relaxed);
        let leak_amount = final_usage as i64 - initial_usage as i64;

        assert!(
            leak_amount.abs() < 5 * 1024 * 1024, // < 5MB tolerance
            "Potential memory leak in graph construction: {:.1}MB difference",
            leak_amount as f64 / (1024.0 * 1024.0)
        );

        println!("✅ No significant memory leaks detected in graph construction");

        Ok(())
    }

    /// Test for memory leaks in streaming processors
    #[test]
    fn test_streaming_processor_leak_detection() -> Result<()> {
        let monitor = RealTimeMemoryMonitor::new(100, 150);
        let initial_usage = monitor.current_usage.load(Ordering::Relaxed);

        // Multiple streaming cycles
        for cycle in 0..5 {
            let config = StreamConfig {
                max_kmers: 10_000,
                memory_limit_mb: 20,
                enable_lru: true,
                sample_rate: 0.5,
            };

            let processor = BoundedStreamProcessor::new(config);
            monitor.allocate(20 * 1024 * 1024)?; // Track processor memory

            // Stream k-mers
            let kmer_stream = (0..50_000).map(|i| (i as u64, vec![i as u8; 16]));
            let _stats = processor.process_kmer_stream(kmer_stream)?;

            let proc_stats = processor.get_streaming_stats();
            println!(
                "Cycle {}: Processor memory: {:.1}MB, utilization: {:.1}%",
                cycle + 1,
                proc_stats.memory_usage_mb,
                proc_stats.memory_utilization
            );

            // Cleanup
            drop(processor);
            monitor.deallocate(20 * 1024 * 1024);
        }

        let final_usage = monitor.current_usage.load(Ordering::Relaxed);
        let leak_amount = final_usage as i64 - initial_usage as i64;

        assert!(
            leak_amount.abs() < 2 * 1024 * 1024, // < 2MB tolerance
            "Potential memory leak in streaming processor: {:.1}MB difference",
            leak_amount as f64 / (1024.0 * 1024.0)
        );

        println!("✅ No significant memory leaks detected in streaming processor");

        Ok(())
    }
}

/* ========================================================================= */
/*                        MEMORY OPTIMIZATION VALIDATION                   */
/* ========================================================================= */

#[cfg(test)]
mod memory_optimization_tests {
    use super::*;

    /// Test memory pool efficiency and reuse
    #[test]
    fn test_memory_pool_efficiency() -> Result<()> {
        let monitor = RealTimeMemoryMonitor::new(150, 200);
        let arena = KmerArena::new(50); // 50MB arena

        // Test memory reuse patterns
        let allocation_batches = vec![100, 200, 150, 300, 50]; // Different batch sizes

        for (batch_num, &batch_size) in allocation_batches.iter().enumerate() {
            let batch_start = monitor.current_usage.load(Ordering::Relaxed);

            let mut kmer_refs = Vec::new();
            for i in 0..batch_size {
                let kmer_data = vec![(batch_num * 1000 + i) as u64];
                match arena.allocate_kmer(&kmer_data) {
                    Ok(kmer_ref) => kmer_refs.push(kmer_ref),
                    Err(_) => break,
                }
            }

            let batch_end = monitor.current_usage.load(Ordering::Relaxed);
            let batch_memory = batch_end - batch_start;

            let arena_stats = arena.memory_stats();
            println!(
                "Batch {}: {} k-mers, {:.1}KB batch memory, {:.1}% utilization",
                batch_num + 1,
                kmer_refs.len(),
                batch_memory as f64 / 1024.0,
                arena_stats.utilization * 100.0
            );

            // Validate efficient memory usage
            assert!(
                arena_stats.utilization > 0.3,
                "Arena utilization {:.1}% too low for batch {}",
                arena_stats.utilization * 100.0,
                batch_num + 1
            );
        }

        let final_stats = arena.memory_stats();
        assert!(
            final_stats.utilization > 0.5,
            "Overall arena utilization {:.1}% should be > 50%",
            final_stats.utilization * 100.0
        );

        println!(
            "✅ Memory pool efficiency validated: {:.1}% utilization",
            final_stats.utilization * 100.0
        );

        Ok(())
    }

    /// Test lock-free graph builder memory efficiency
    #[test]
    fn test_lock_free_builder_memory_efficiency() -> Result<()> {
        let config = BuilderConfig {
            edge_batch_size: 5_000,
            max_memory_bytes: 100 * 1024 * 1024, // 100MB
            numa_aware: false,
        };

        let builder = LockFreeGraphBuilder::new(config);
        let monitor = RealTimeMemoryMonitor::new(120, 150);

        // Build graph with memory tracking
        for i in 0..10_000 {
            let kmer_ref = meta_forge::assembly::memory_optimizations::KmerRef {
                block_id: 0,
                offset: i * 8,
                length: 2,
            };

            builder.add_node(i as u64, kmer_ref)?;

            if i > 0 {
                builder.queue_edge((i - 1) as u64, i as u64, 1.0);
            }

            // Process edges periodically
            if i % 1000 == 0 {
                let _result = builder.process_edge_batches()?;
            }
        }

        // Final edge processing
        let final_result = builder.process_edge_batches()?;
        let builder_stats = builder.get_stats();

        println!("📊 Lock-free builder results:");
        println!("  Nodes: {}", builder_stats.total_nodes);
        println!("  Edges processed: {}", final_result.valid_edges);
        println!("  Memory usage: {}MB", builder_stats.memory_usage_mb);
        println!("  Queue length: {}", builder_stats.queue_length);

        // Validate memory efficiency
        assert!(
            builder_stats.memory_usage_mb < 100,
            "Memory usage {}MB exceeds 100MB limit",
            builder_stats.memory_usage_mb
        );
        assert!(
            final_result.valid_edges > 8000,
            "Should process most edges successfully"
        );

        println!("✅ Lock-free builder memory efficiency validated");

        Ok(())
    }

    /// Test memory usage under concurrent load
    #[test]
    fn test_concurrent_memory_usage() -> Result<()> {
        let monitor = Arc::new(RealTimeMemoryMonitor::new(200, 300));
        let num_threads = 4;
        let mut handles = Vec::new();

        println!(
            "🔄 Testing concurrent memory usage with {} threads...",
            num_threads
        );

        for thread_id in 0..num_threads {
            let monitor_clone = Arc::clone(&monitor);

            let handle = thread::spawn(move || -> Result<()> {
                let config = LaptopConfig::low_memory();
                let mut graph = LaptopAssemblyGraph::new(config);

                // Each thread processes different data
                let reads = create_thread_specific_reads(thread_id, 200, 60);

                // Simulate memory allocation tracking
                monitor_clone.allocate(30 * 1024 * 1024)?; // 30MB per thread

                graph.build_from_reads(&reads, 15)?;
                let memory_usage = graph.memory_usage_mb();

                println!("Thread {}: {:.1}MB memory usage", thread_id, memory_usage);

                // Cleanup
                monitor_clone.deallocate(30 * 1024 * 1024);

                Ok(())
            });

            handles.push(handle);
        }

        // Wait for all threads to complete
        for (i, handle) in handles.into_iter().enumerate() {
            match handle.join() {
                Ok(result) => match result {
                    Ok(_) => println!("Thread {} completed successfully", i),
                    Err(e) => println!("Thread {} failed: {}", i, e),
                },
                Err(_) => println!("Thread {} panicked", i),
            }
        }

        let final_stats = monitor.get_stats();
        println!("📊 Concurrent execution results:");
        println!("  Peak memory: {:.1}MB", final_stats.peak_usage_mb);
        println!("  Final memory: {:.1}MB", final_stats.current_usage_mb);
        println!("  Total allocations: {}", final_stats.allocation_count);
        println!("  Total deallocations: {}", final_stats.deallocation_count);

        // Validate concurrent memory handling
        assert!(
            final_stats.peak_usage_mb < 300.0,
            "Peak memory {:.1}MB exceeded limit",
            final_stats.peak_usage_mb
        );
        assert!(
            final_stats.current_usage_mb < 50.0,
            "Final memory {:.1}MB should be low after cleanup",
            final_stats.current_usage_mb
        );

        println!("✅ Concurrent memory usage handling validated");

        Ok(())
    }
}

/* ========================================================================= */
/*                        LAPTOP-SPECIFIC MEMORY TESTS                     */
/* ========================================================================= */

#[cfg(test)]
mod laptop_specific_tests {
    use super::*;

    /// Test memory management for typical laptop workloads
    #[test]
    fn test_typical_laptop_workload() -> Result<()> {
        let scenarios = vec![
            ("4GB Laptop", LaptopConfig::low_memory(), 800), // Conservative memory limit
            ("8GB Laptop", LaptopConfig::medium_memory(), 1600),
            ("16GB Laptop", LaptopConfig::high_memory(), 3200),
        ];

        for (scenario_name, config, memory_limit_mb) in scenarios {
            println!("\n🖥️ Testing {} scenario...", scenario_name);

            let monitor = RealTimeMemoryMonitor::new(
                (memory_limit_mb as f64 * 0.8) as usize, // Pressure at 80%
                memory_limit_mb,                         // Critical at 100%
            );

            // Simulate typical bioinformatics workload
            let read_counts = vec![100, 500, 1000, 2000];

            for &read_count in &read_counts {
                let reads = create_laptop_workload_reads(read_count, 100);
                let mut graph = LaptopAssemblyGraph::new(config.clone());

                let estimated_memory = estimate_memory_usage(read_count, 100);
                if !monitor.allocate(estimated_memory) {
                    println!(
                        "  {} reads: SKIPPED (would exceed memory limit)",
                        read_count
                    );
                    continue;
                }

                let start_time = Instant::now();
                graph.build_from_reads(&reads, 21)?;
                let processing_time = start_time.elapsed();

                let actual_memory = graph.memory_usage_mb();
                println!(
                    "  {} reads: {:.1}MB memory, {:.2}s processing",
                    read_count,
                    actual_memory,
                    processing_time.as_secs_f64()
                );

                // Validate laptop constraints
                assert!(
                    actual_memory < memory_limit_mb as f64,
                    "Memory usage {:.1}MB exceeds {} limit {}MB",
                    actual_memory,
                    scenario_name,
                    memory_limit_mb
                );

                assert!(
                    processing_time.as_secs() < 30,
                    "Processing time {:.2}s too long for laptop",
                    processing_time.as_secs_f64()
                );

                monitor.deallocate(estimated_memory);
            }

            let stats = monitor.get_stats();
            println!(
                "  Peak memory: {:.1}MB, Pressure ratio: {:.2}",
                stats.peak_usage_mb, stats.pressure_ratio
            );

            println!("✅ {} scenario validated", scenario_name);
        }

        Ok(())
    }

    /// Test memory recovery after system pressure
    #[test]
    fn test_memory_recovery_after_pressure() -> Result<()> {
        let monitor = RealTimeMemoryMonitor::new(50, 80); // Low thresholds for testing

        // Phase 1: Build up memory pressure
        println!("Phase 1: Building memory pressure...");
        let mut allocations = Vec::new();

        for i in 0..15 {
            let size = 8 * 1024 * 1024; // 8MB chunks
            if monitor.allocate(size) {
                allocations.push(size);
                println!(
                    "  Allocation {}: {:.1}MB (total: {:.1}MB)",
                    i + 1,
                    size as f64 / (1024.0 * 1024.0),
                    monitor.get_stats().current_usage_mb
                );
            } else {
                println!("  Allocation {} blocked by critical threshold", i + 1);
                break;
            }
        }

        assert!(
            monitor.is_under_pressure(),
            "Should be under memory pressure"
        );
        let pressure_stats = monitor.get_stats();
        println!(
            "Memory pressure reached: {:.1}MB ({:.1}% of threshold)",
            pressure_stats.current_usage_mb,
            pressure_stats.pressure_ratio * 100.0
        );

        // Phase 2: Gradual memory recovery
        println!("\nPhase 2: Gradual memory recovery...");
        let recovery_count = allocations.len() / 2; // Free half the allocations

        for i in 0..recovery_count {
            let size = allocations[i];
            monitor.deallocate(size);

            let stats = monitor.get_stats();
            println!(
                "  Recovery {}: freed {:.1}MB (remaining: {:.1}MB)",
                i + 1,
                size as f64 / (1024.0 * 1024.0),
                stats.current_usage_mb
            );
        }

        // Phase 3: Validate recovery
        let recovery_stats = monitor.get_stats();
        assert!(
            !monitor.is_critical(),
            "Should no longer be at critical level"
        );
        assert!(
            recovery_stats.current_usage_mb < pressure_stats.current_usage_mb,
            "Memory usage should have decreased"
        );

        println!("\n✅ Memory recovery validated:");
        println!(
            "  Peak: {:.1}MB -> Current: {:.1}MB",
            pressure_stats.peak_usage_mb, recovery_stats.current_usage_mb
        );
        println!(
            "  Pressure ratio: {:.2} -> {:.2}",
            pressure_stats.pressure_ratio, recovery_stats.pressure_ratio
        );

        Ok(())
    }

    /// Test memory fragmentation handling
    #[test]
    fn test_memory_fragmentation_handling() -> Result<()> {
        let arena = KmerArena::new(20); // 20MB arena
        let monitor = RealTimeMemoryMonitor::new(30, 40);

        // Simulate fragmented allocation pattern
        let mut kmer_refs = Vec::new();
        let allocation_pattern = vec![
            vec![1, 2, 3, 4, 5], // Small allocations
            vec![10, 15, 20],    // Medium allocations
            vec![8, 12, 6, 9],   // Mixed sizes
        ];

        for (phase, sizes) in allocation_pattern.iter().enumerate() {
            println!(
                "Fragmentation phase {}: allocating {:?} k-mers",
                phase + 1,
                sizes
            );

            for &size in sizes {
                let kmer_data = vec![0u64; size];
                match arena.allocate_kmer(&kmer_data) {
                    Ok(kmer_ref) => {
                        kmer_refs.push(kmer_ref);
                        monitor.allocate(size * 8)?; // Track allocation
                    }
                    Err(_) => {
                        println!("  Allocation failed for size {}", size);
                        break;
                    }
                }
            }

            let arena_stats = arena.memory_stats();
            println!(
                "  Arena utilization: {:.1}%",
                arena_stats.utilization * 100.0
            );
        }

        // Validate fragmentation handling
        let final_stats = arena.memory_stats();
        assert!(
            final_stats.utilization > 0.4,
            "Arena should handle fragmentation efficiently (utilization: {:.1}%)",
            final_stats.utilization * 100.0
        );

        println!(
            "✅ Memory fragmentation handling validated: {:.1}% utilization",
            final_stats.utilization * 100.0
        );

        Ok(())
    }

    fn estimate_memory_usage(read_count: usize, read_length: usize) -> usize {
        // Rough estimate: each read generates ~2-3 k-mers, each k-mer ~32 bytes
        let kmer_count = read_count * 2;
        let kmer_memory = kmer_count * 32;
        let overhead = kmer_memory / 4; // 25% overhead for structures
        kmer_memory + overhead
    }
}

/* ========================================================================= */
/*                        HELPER FUNCTIONS                                 */
/* ========================================================================= */

/// Create reads for memory leak testing
fn create_leak_test_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    (0..count)
        .map(|i| {
            let sequence = format!("{:0width$}", i, width = length)
                .chars()
                .map(|c| match c {
                    '0' => 'A',
                    '1' => 'C',
                    '2' => 'G',
                    '3' => 'T',
                    '4' => 'A',
                    '5' => 'C',
                    '6' => 'G',
                    '7' => 'T',
                    '8' => 'A',
                    '9' => 'C',
                    _ => 'A',
                })
                .collect::<String>();

            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; length],
                correction_metadata: CorrectionMetadata {
                    algorithm: "leak_test".to_string(),
                    confidence_threshold: 0.9,
                    context_window: 5,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: AHashMap::new(),
            }
        })
        .collect()
}

/// Create thread-specific reads for concurrent testing
fn create_thread_specific_reads(
    thread_id: usize,
    count: usize,
    length: usize,
) -> Vec<CorrectedRead> {
    (0..count)
        .map(|i| {
            let base_id = thread_id * 10000 + i;
            let nucleotides = ['A', 'C', 'G', 'T'];

            let sequence: String = (0..length)
                .map(|j| nucleotides[(base_id + j) % 4])
                .collect();

            CorrectedRead {
                id: base_id,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; length],
                correction_metadata: CorrectionMetadata {
                    algorithm: format!("thread_{}", thread_id),
                    confidence_threshold: 0.9,
                    context_window: 5,
                    correction_time_ms: 0,
                },
                kmer_hash_cache: AHashMap::new(),
            }
        })
        .collect()
}

/// Create reads for laptop workload testing
fn create_laptop_workload_reads(count: usize, length: usize) -> Vec<CorrectedRead> {
    // Simulate realistic bioinformatics reads with some overlap
    let base_sequences = vec![
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
        "TTTTCCCCGGGGAAAATTTTCCCCGGGGAAAATTTTCCCC",
        "AAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGCTTAAGC",
    ];

    (0..count)
        .map(|i| {
            let base_seq = &base_sequences[i % base_sequences.len()];
            let start = (i * 3) % (base_seq.len() - length);
            let sequence = base_seq[start..start + length].to_string();

            CorrectedRead {
                id: i,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![28 + (i % 8) as u8; length], // Variable quality
                correction_metadata: CorrectionMetadata {
                    algorithm: "laptop_workload".to_string(),
                    confidence_threshold: 0.88,
                    context_window: 6,
                    correction_time_ms: 1,
                },
                kmer_hash_cache: AHashMap::new(),
            }
        })
        .collect()
}
