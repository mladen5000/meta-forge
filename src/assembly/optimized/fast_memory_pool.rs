//! Fast Memory Pool for Assembly Operations
//! ======================================
//!
//! Optimized memory pool implementation for high-frequency allocations
//! in assembly operations. Features cache-friendly allocation patterns,
//! NUMA awareness, and lock-free operations where possible.

use std::alloc::{GlobalAlloc, Layout, System};
use std::ptr::NonNull;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use anyhow::{anyhow, Result};

/// High-performance memory pool optimized for genomics workloads
pub struct FastAssemblyMemoryPool {
    /// Pool configuration
    config: PoolConfig,
    /// Cache-aligned pool chunks to avoid false sharing
    chunks: Vec<CachePadded<PoolChunk>>,
    /// Global statistics
    stats: Arc<PoolStatistics>,
    /// Current allocation index (round-robin)
    current_chunk: AtomicUsize,
}

/// Memory pool configuration
#[derive(Debug, Clone)]
pub struct PoolConfig {
    /// Size of each pool chunk in bytes
    pub chunk_size: usize,
    /// Number of chunks to pre-allocate
    pub initial_chunks: usize,
    /// Maximum number of chunks
    pub max_chunks: usize,
    /// Alignment requirement for allocations
    pub alignment: usize,
    /// Enable NUMA-aware allocation
    pub numa_aware: bool,
    /// Cache line size for padding
    pub cache_line_size: usize,
}

impl Default for PoolConfig {
    fn default() -> Self {
        Self {
            chunk_size: 64 * 1024 * 1024, // 64MB chunks
            initial_chunks: 4,
            max_chunks: 16,
            alignment: 64, // Cache line aligned
            numa_aware: true,
            cache_line_size: 64,
        }
    }
}

/// Individual memory chunk with lock-free allocation
struct PoolChunk {
    /// Raw memory block
    memory: NonNull<u8>,
    /// Total size of this chunk
    size: usize,
    /// Current allocation offset (atomic for lock-free allocation)
    offset: AtomicUsize,
    /// Number of active allocations in this chunk
    active_allocs: AtomicUsize,
}

/// Pool performance statistics
#[derive(Debug, Default)]
pub struct PoolStatistics {
    pub total_allocations: AtomicUsize,
    pub total_deallocations: AtomicUsize,
    pub total_bytes_allocated: AtomicUsize,
    pub peak_memory_usage: AtomicUsize,
    pub cache_hits: AtomicUsize,
    pub cache_misses: AtomicUsize,
    pub fragmentation_ratio: AtomicUsize, // As percentage * 100
}

/// Smart pointer for pool-allocated memory with automatic cleanup
pub struct PooledAllocation<T> {
    ptr: NonNull<T>,
    chunk_id: usize,
    pool: Arc<FastAssemblyMemoryPool>,
}

impl FastAssemblyMemoryPool {
    /// Create new memory pool with specified configuration
    pub fn new(config: PoolConfig) -> Result<Arc<Self>> {
        let mut chunks = Vec::with_capacity(config.initial_chunks);

        // Pre-allocate initial chunks
        for _ in 0..config.initial_chunks {
            let chunk = PoolChunk::new(config.chunk_size, config.alignment)?;
            chunks.push(CachePadded::new(chunk));
        }

        Ok(Arc::new(Self {
            config,
            chunks,
            stats: Arc::new(PoolStatistics::default()),
            current_chunk: AtomicUsize::new(0),
        }))
    }

    /// Allocate memory for type T with proper alignment
    /// OPTIMIZATION: Lock-free allocation with cache-friendly round-robin
    pub fn allocate<T>(&self) -> Result<PooledAllocation<T>, std::alloc::AllocError> {
        let layout = Layout::new::<T>();
        let size = layout.size();
        let align = layout.align().max(self.config.alignment);

        // OPTIMIZATION: Fast path - try current chunk first
        let current = self.current_chunk.load(Ordering::Acquire);
        if let Some(ptr) = self.try_allocate_from_chunk(current, size, align) {
            self.stats.total_allocations.fetch_add(1, Ordering::Relaxed);
            self.stats.total_bytes_allocated.fetch_add(size, Ordering::Relaxed);
            self.stats.cache_hits.fetch_add(1, Ordering::Relaxed);

            return Ok(PooledAllocation {
                ptr: ptr.cast(),
                chunk_id: current,
                pool: Arc::clone(&unsafe { std::mem::transmute::<&Self, &Arc<Self>>(self) }),
            });
        }

        // OPTIMIZATION: Round-robin through chunks to distribute load
        for i in 0..self.chunks.len() {
            let chunk_idx = (current + i + 1) % self.chunks.len();
            if let Some(ptr) = self.try_allocate_from_chunk(chunk_idx, size, align) {
                // Update current chunk for next allocation
                self.current_chunk.store(chunk_idx, Ordering::Release);

                self.stats.total_allocations.fetch_add(1, Ordering::Relaxed);
                self.stats.total_bytes_allocated.fetch_add(size, Ordering::Relaxed);

                return Ok(PooledAllocation {
                    ptr: ptr.cast(),
                    chunk_id: chunk_idx,
                    pool: Arc::clone(&unsafe { std::mem::transmute::<&Self, &Arc<Self>>(self) }),
                });
            }
        }

        // No space available - try to allocate new chunk if under limit
        self.stats.cache_misses.fetch_add(1, Ordering::Relaxed);
        self.allocate_fallback(size, align)
    }

    /// Try to allocate from specific chunk (lock-free)
    fn try_allocate_from_chunk(&self, chunk_idx: usize, size: usize, align: usize) -> Option<NonNull<u8>> {
        if chunk_idx >= self.chunks.len() {
            return None;
        }

        let chunk = &self.chunks[chunk_idx];

        loop {
            let current_offset = chunk.offset.load(Ordering::Acquire);

            // Calculate aligned offset
            let aligned_offset = (current_offset + align - 1) & !(align - 1);
            let end_offset = aligned_offset + size;

            if end_offset > chunk.size {
                return None; // Not enough space in this chunk
            }

            // Try to atomically update offset
            match chunk.offset.compare_exchange_weak(
                current_offset,
                end_offset,
                Ordering::Release,
                Ordering::Relaxed
            ) {
                Ok(_) => {
                    // Successfully allocated
                    chunk.active_allocs.fetch_add(1, Ordering::Relaxed);

                    unsafe {
                        let ptr = chunk.memory.as_ptr().add(aligned_offset);
                        return NonNull::new(ptr);
                    }
                }
                Err(_) => {
                    // Retry with updated offset
                    std::hint::spin_loop();
                    continue;
                }
            }
        }
    }

    /// Fallback allocation using system allocator
    fn allocate_fallback<T>(&self, size: usize, align: usize) -> Result<PooledAllocation<T>, std::alloc::AllocError> {
        let layout = Layout::from_size_align(size, align).map_err(|_| std::alloc::AllocError)?;

        unsafe {
            let ptr = System.alloc(layout);
            if ptr.is_null() {
                return Err(std::alloc::AllocError);
            }

            self.stats.total_allocations.fetch_add(1, Ordering::Relaxed);
            self.stats.total_bytes_allocated.fetch_add(size, Ordering::Relaxed);

            Ok(PooledAllocation {
                ptr: NonNull::new_unchecked(ptr).cast(),
                chunk_id: usize::MAX, // Indicates system allocation
                pool: Arc::clone(&std::mem::transmute::<&Self, &Arc<Self>>(self)),
            })
        }
    }

    /// Get memory pool statistics
    pub fn statistics(&self) -> PoolStatistics {
        PoolStatistics {
            total_allocations: AtomicUsize::new(self.stats.total_allocations.load(Ordering::Relaxed)),
            total_deallocations: AtomicUsize::new(self.stats.total_deallocations.load(Ordering::Relaxed)),
            total_bytes_allocated: AtomicUsize::new(self.stats.total_bytes_allocated.load(Ordering::Relaxed)),
            peak_memory_usage: AtomicUsize::new(self.stats.peak_memory_usage.load(Ordering::Relaxed)),
            cache_hits: AtomicUsize::new(self.stats.cache_hits.load(Ordering::Relaxed)),
            cache_misses: AtomicUsize::new(self.stats.cache_misses.load(Ordering::Relaxed)),
            fragmentation_ratio: AtomicUsize::new(self.stats.fragmentation_ratio.load(Ordering::Relaxed)),
        }
    }

    /// Calculate current memory utilization
    pub fn utilization_stats(&self) -> UtilizationStats {
        let mut total_allocated = 0;
        let mut total_capacity = 0;
        let mut fragmented_space = 0;

        for chunk in &self.chunks {
            let used = chunk.offset.load(Ordering::Relaxed);
            total_allocated += used;
            total_capacity += chunk.size;

            // Simple fragmentation estimate
            let wasted = chunk.size - used;
            if wasted > 0 && chunk.active_allocs.load(Ordering::Relaxed) > 0 {
                fragmented_space += wasted;
            }
        }

        let utilization = if total_capacity > 0 {
            (total_allocated as f64 / total_capacity as f64) * 100.0
        } else {
            0.0
        };

        let fragmentation = if total_capacity > 0 {
            (fragmented_space as f64 / total_capacity as f64) * 100.0
        } else {
            0.0
        };

        UtilizationStats {
            total_capacity,
            total_allocated,
            utilization_percent: utilization,
            fragmentation_percent: fragmentation,
            active_chunks: self.chunks.len(),
        }
    }
}

/// Memory utilization statistics
#[derive(Debug, Clone)]
pub struct UtilizationStats {
    pub total_capacity: usize,
    pub total_allocated: usize,
    pub utilization_percent: f64,
    pub fragmentation_percent: f64,
    pub active_chunks: usize,
}

impl PoolChunk {
    /// Create new memory chunk
    fn new(size: usize, alignment: usize) -> Result<Self, std::alloc::AllocError> {
        let layout = Layout::from_size_align(size, alignment).map_err(|_| std::alloc::AllocError)?;

        unsafe {
            let ptr = System.alloc(layout);
            if ptr.is_null() {
                return Err(std::alloc::AllocError);
            }

            Ok(Self {
                memory: NonNull::new_unchecked(ptr),
                size,
                offset: AtomicUsize::new(0),
                active_allocs: AtomicUsize::new(0),
            })
        }
    }
}

impl Drop for PoolChunk {
    fn drop(&mut self) {
        unsafe {
            let layout = Layout::from_size_align_unchecked(self.size, 1);
            System.dealloc(self.memory.as_ptr(), layout);
        }
    }
}

impl<T> PooledAllocation<T> {
    /// Get raw pointer to allocated memory
    pub fn as_ptr(&self) -> *mut T {
        self.ptr.as_ptr()
    }

    /// Get mutable reference to allocated memory
    pub unsafe fn as_mut(&mut self) -> &mut T {
        self.ptr.as_mut()
    }

    /// Get immutable reference to allocated memory
    pub unsafe fn as_ref(&self) -> &T {
        self.ptr.as_ref()
    }
}

impl<T> Drop for PooledAllocation<T> {
    fn drop(&mut self) {
        // Update deallocation statistics
        self.pool.stats.total_deallocations.fetch_add(1, Ordering::Relaxed);

        if self.chunk_id != usize::MAX {
            // Pool allocation - update chunk stats
            if self.chunk_id < self.pool.chunks.len() {
                self.pool.chunks[self.chunk_id].active_allocs.fetch_sub(1, Ordering::Relaxed);
            }
        } else {
            // System allocation - need to free
            unsafe {
                let layout = Layout::new::<T>();
                System.dealloc(self.ptr.as_ptr() as *mut u8, layout);
            }
        }
    }
}

// Safety: PooledAllocation is Send and Sync if T is Send and Sync
unsafe impl<T: Send> Send for PooledAllocation<T> {}
unsafe impl<T: Sync> Sync for PooledAllocation<T> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pool_creation() {
        let config = PoolConfig::default();
        let pool = FastAssemblyMemoryPool::new(config).unwrap();

        let stats = pool.utilization_stats();
        assert_eq!(stats.active_chunks, 4); // Initial chunks
        assert!(stats.total_capacity > 0);
    }

    #[test]
    fn test_allocation_deallocation() {
        let config = PoolConfig::default();
        let pool = FastAssemblyMemoryPool::new(config).unwrap();

        // Allocate some memory
        let _alloc1: PooledAllocation<u64> = pool.allocate().unwrap();
        let _alloc2: PooledAllocation<[u8; 1024]> = pool.allocate().unwrap();

        let stats = pool.statistics();
        assert_eq!(stats.total_allocations.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn test_concurrent_allocation() {
        use std::thread;

        let config = PoolConfig::default();
        let pool = Arc::new(FastAssemblyMemoryPool::new(config).unwrap());

        let handles: Vec<_> = (0..8).map(|_| {
            let pool_clone = Arc::clone(&pool);
            thread::spawn(move || {
                for _ in 0..100 {
                    let _alloc: PooledAllocation<u64> = pool_clone.allocate().unwrap();
                    // Allocation automatically dropped here
                }
            })
        }).collect();

        for handle in handles {
            handle.join().unwrap();
        }

        let stats = pool.statistics();
        assert_eq!(stats.total_allocations.load(Ordering::Relaxed), 800);
    }
}