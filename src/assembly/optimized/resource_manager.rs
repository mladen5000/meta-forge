//! Adaptive Resource Manager
//! ========================
//!
//! Dynamic system resource management with memory pressure detection,
//! CPU utilization monitoring, and adaptive workload balancing.

use crate::assembly::laptop_assembly::LaptopConfig as OptimizedConfig;
use anyhow::Result;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

/// Adaptive resource manager for dynamic system monitoring and optimization
pub struct AdaptiveResourceManager {
    config: OptimizedConfig,
    system_monitor: Arc<SystemMonitor>,
    memory_manager: Arc<MemoryManager>,
    cpu_manager: Arc<CpuManager>,
    metrics: ResourceMetrics,
    monitoring_active: Arc<AtomicBool>,
    _monitor_thread: Option<thread::JoinHandle<()>>,
}

impl AdaptiveResourceManager {
    /// Create new adaptive resource manager
    pub fn new(config: &OptimizedConfig) -> Self {
        let system_monitor = Arc::new(SystemMonitor::new());
        let memory_manager = Arc::new(MemoryManager::new(config.memory_budget_mb));
        let cpu_manager = Arc::new(CpuManager::new(config.cpu_cores));

        let monitoring_active = Arc::new(AtomicBool::new(true));

        // Start background monitoring thread
        let monitor_thread = Self::start_monitoring_thread(
            system_monitor.clone(),
            memory_manager.clone(),
            cpu_manager.clone(),
            monitoring_active.clone(),
        );

        Self {
            config: config.clone(),
            system_monitor,
            memory_manager,
            cpu_manager,
            metrics: ResourceMetrics::default(),
            monitoring_active,
            _monitor_thread: Some(monitor_thread),
        }
    }

    /// Get current memory pressure (0.0-1.0)
    pub fn memory_pressure(&self) -> f64 {
        self.memory_manager.memory_pressure()
    }

    /// Get available memory in bytes
    pub fn available_memory(&self) -> usize {
        self.memory_manager.available_memory()
    }

    /// Get current CPU utilization (0.0-1.0)
    pub fn cpu_utilization(&self) -> f64 {
        self.cpu_manager.utilization()
    }

    /// Check if system is under resource pressure
    pub fn is_under_pressure(&self) -> bool {
        self.memory_pressure() > 0.8 || self.cpu_utilization() > 0.9
    }

    /// Get optimal chunk size based on current resources
    pub fn optimal_chunk_size(&self, base_size: usize) -> usize {
        let memory_factor = 1.0 - self.memory_pressure().max(0.0).min(1.0);
        let cpu_factor = 1.0 - self.cpu_utilization().max(0.0).min(1.0);

        let resource_factor = (memory_factor + cpu_factor) / 2.0;
        let adjusted_size = base_size as f64 * (0.5 + resource_factor * 0.5);

        adjusted_size.max(100.0).min(base_size as f64 * 2.0) as usize
    }

    /// Get optimal number of parallel workers
    pub fn optimal_workers(&self) -> usize {
        let cpu_usage = self.cpu_utilization();
        let memory_pressure = self.memory_pressure();

        let available_cpu_ratio = (1.0 - cpu_usage).max(0.1);
        let memory_safety_ratio = (1.0 - memory_pressure).max(0.2);

        let optimal_ratio = available_cpu_ratio.min(memory_safety_ratio);
        let optimal_workers = (self.config.cpu_cores as f64 * optimal_ratio) as usize;

        optimal_workers.max(1).min(self.config.cpu_cores)
    }

    /// Request memory allocation with tracking
    pub fn allocate_memory(&self, size: usize, purpose: &str) -> Result<MemoryAllocation> {
        self.memory_manager.allocate(size, purpose)
    }

    /// Get comprehensive resource metrics
    pub fn metrics(&self) -> ResourceMetricsSnapshot {
        ResourceMetricsSnapshot {
            memory_pressure: self.memory_pressure(),
            cpu_utilization: self.cpu_utilization(),
            available_memory: self.available_memory(),
            total_memory: self.memory_manager.total_memory(),
            active_allocations: self.memory_manager.active_allocations(),
            peak_memory_usage: self.memory_manager.peak_usage(),
            system_load: self.system_monitor.load_average(),
            uptime: self.system_monitor.uptime(),
        }
    }

    /// Start background monitoring thread
    fn start_monitoring_thread(
        system_monitor: Arc<SystemMonitor>,
        memory_manager: Arc<MemoryManager>,
        cpu_manager: Arc<CpuManager>,
        active: Arc<AtomicBool>,
    ) -> thread::JoinHandle<()> {
        thread::spawn(move || {
            while active.load(Ordering::Relaxed) {
                // Update system metrics
                system_monitor.update();
                memory_manager.update();
                cpu_manager.update();

                // Monitor every 100ms
                thread::sleep(Duration::from_millis(100));
            }
        })
    }
}

impl Drop for AdaptiveResourceManager {
    fn drop(&mut self) {
        self.monitoring_active.store(false, Ordering::Relaxed);
    }
}

/// System monitoring and metrics collection
pub struct SystemMonitor {
    start_time: Instant,
    load_average: AtomicU64, // Stored as fixed-point * 1000
    last_update: Mutex<Instant>,
}

impl SystemMonitor {
    pub fn new() -> Self {
        Self {
            start_time: Instant::now(),
            load_average: AtomicU64::new(0),
            last_update: Mutex::new(Instant::now()),
        }
    }

    /// Update system metrics
    pub fn update(&self) {
        let mut last_update = self.last_update.lock().unwrap();
        let now = Instant::now();

        // Update every second at most
        if now.duration_since(*last_update) < Duration::from_secs(1) {
            return;
        }

        // Get system load (simplified - in production would use system calls)
        let load = self.get_system_load();
        self.load_average.store((load * 1000.0) as u64, Ordering::Relaxed);

        *last_update = now;
    }

    /// Get system load average
    pub fn load_average(&self) -> f64 {
        self.load_average.load(Ordering::Relaxed) as f64 / 1000.0
    }

    /// Get system uptime
    pub fn uptime(&self) -> Duration {
        self.start_time.elapsed()
    }

    /// Get system load (platform-specific implementation)
    fn get_system_load(&self) -> f64 {
        // Simplified implementation - in production would use:
        // - Linux: /proc/loadavg
        // - macOS: sysctl
        // - Windows: Performance counters

        // For now, return a simulated load based on thread activity
        #[cfg(unix)]
        {
            self.get_unix_load()
        }
        #[cfg(not(unix))]
        {
            0.5 // Default value
        }
    }

    #[cfg(unix)]
    fn get_unix_load(&self) -> f64 {
        // Try to read /proc/loadavg on Linux
        if let Ok(content) = std::fs::read_to_string("/proc/loadavg") {
            if let Some(first_field) = content.split_whitespace().next() {
                return first_field.parse().unwrap_or(0.0);
            }
        }

        // Fallback for macOS or other Unix systems
        0.3 // Simulated value
    }
}

/// Memory management with allocation tracking
pub struct MemoryManager {
    total_memory_mb: usize,
    allocated_memory: AtomicUsize,
    peak_usage: AtomicUsize,
    active_allocations: AtomicUsize,
    allocation_map: Mutex<std::collections::HashMap<usize, AllocationInfo>>,
    next_id: AtomicUsize,
}

impl MemoryManager {
    pub fn new(total_memory_mb: usize) -> Self {
        Self {
            total_memory_mb,
            allocated_memory: AtomicUsize::new(0),
            peak_usage: AtomicUsize::new(0),
            active_allocations: AtomicUsize::new(0),
            allocation_map: Mutex::new(std::collections::HashMap::new()),
            next_id: AtomicUsize::new(1),
        }
    }

    /// Get current memory pressure (0.0-1.0)
    pub fn memory_pressure(&self) -> f64 {
        let allocated = self.allocated_memory.load(Ordering::Relaxed);
        let total = self.total_memory_mb * 1024 * 1024; // Convert to bytes

        if total == 0 {
            0.0
        } else {
            (allocated as f64) / (total as f64)
        }
    }

    /// Get available memory in bytes
    pub fn available_memory(&self) -> usize {
        let allocated = self.allocated_memory.load(Ordering::Relaxed);
        let total = self.total_memory_mb * 1024 * 1024;
        total.saturating_sub(allocated)
    }

    /// Get total memory in bytes
    pub fn total_memory(&self) -> usize {
        self.total_memory_mb * 1024 * 1024
    }

    /// Get number of active allocations
    pub fn active_allocations(&self) -> usize {
        self.active_allocations.load(Ordering::Relaxed)
    }

    /// Get peak memory usage
    pub fn peak_usage(&self) -> usize {
        self.peak_usage.load(Ordering::Relaxed)
    }

    /// Allocate memory with tracking
    pub fn allocate(&self, size: usize, purpose: &str) -> Result<MemoryAllocation> {
        let current_allocated = self.allocated_memory.load(Ordering::Relaxed);
        let total_memory = self.total_memory();

        // Check if allocation would exceed budget
        if current_allocated + size > total_memory {
            return Err(anyhow::anyhow!(
                "Memory allocation would exceed budget: {} + {} > {}",
                current_allocated, size, total_memory
            ));
        }

        let id = self.next_id.fetch_add(1, Ordering::Relaxed);
        let allocation_info = AllocationInfo {
            size,
            purpose: purpose.to_string(),
            timestamp: Instant::now(),
        };

        // Track allocation
        {
            let mut map = self.allocation_map.lock().unwrap();
            map.insert(id, allocation_info);
        }

        let new_allocated = self.allocated_memory.fetch_add(size, Ordering::Relaxed) + size;
        self.active_allocations.fetch_add(1, Ordering::Relaxed);

        // Update peak usage
        loop {
            let current_peak = self.peak_usage.load(Ordering::Relaxed);
            if new_allocated <= current_peak {
                break;
            }
            if self.peak_usage.compare_exchange_weak(
                current_peak,
                new_allocated,
                Ordering::Relaxed,
                Ordering::Relaxed
            ).is_ok() {
                break;
            }
        }

        Ok(MemoryAllocation {
            id,
            size,
            manager: self as *const Self,
        })
    }

    /// Deallocate memory
    fn deallocate(&self, id: usize, size: usize) {
        self.allocated_memory.fetch_sub(size, Ordering::Relaxed);
        self.active_allocations.fetch_sub(1, Ordering::Relaxed);

        let mut map = self.allocation_map.lock().unwrap();
        map.remove(&id);
    }

    /// Update memory statistics
    pub fn update(&self) {
        // Could collect system memory stats here
        // For now, we rely on our own allocation tracking
    }
}

/// CPU utilization monitoring
pub struct CpuManager {
    cores: usize,
    utilization: AtomicU64, // Stored as fixed-point * 1000
    last_update: Mutex<Instant>,
}

impl CpuManager {
    pub fn new(cores: usize) -> Self {
        Self {
            cores,
            utilization: AtomicU64::new(0),
            last_update: Mutex::new(Instant::now()),
        }
    }

    /// Get current CPU utilization (0.0-1.0)
    pub fn utilization(&self) -> f64 {
        self.utilization.load(Ordering::Relaxed) as f64 / 1000.0
    }

    /// Update CPU metrics
    pub fn update(&self) {
        let mut last_update = self.last_update.lock().unwrap();
        let now = Instant::now();

        // Update every second at most
        if now.duration_since(*last_update) < Duration::from_secs(1) {
            return;
        }

        // Get CPU utilization (simplified)
        let usage = self.get_cpu_usage();
        self.utilization.store((usage * 1000.0) as u64, Ordering::Relaxed);

        *last_update = now;
    }

    /// Get CPU usage (platform-specific implementation)
    fn get_cpu_usage(&self) -> f64 {
        // Simplified implementation - in production would use:
        // - Linux: /proc/stat
        // - macOS: host_processor_info
        // - Windows: GetSystemTimes

        // For now, return a simulated usage based on thread count
        let active_threads = thread::available_parallelism()
            .map(|p| p.get())
            .unwrap_or(self.cores);

        if active_threads > self.cores {
            0.8 // High utilization if more threads than cores
        } else {
            active_threads as f64 / self.cores as f64 * 0.6 // Moderate utilization
        }
    }
}

/// Memory allocation with automatic cleanup
pub struct MemoryAllocation {
    id: usize,
    size: usize,
    manager: *const MemoryManager,
}

impl Drop for MemoryAllocation {
    fn drop(&mut self) {
        unsafe {
            (*self.manager).deallocate(self.id, self.size);
        }
    }
}

/// Information about a memory allocation
#[derive(Debug, Clone)]
struct AllocationInfo {
    size: usize,
    purpose: String,
    timestamp: Instant,
}

/// Resource metrics
#[derive(Debug, Default)]
pub struct ResourceMetrics {
    pub memory_pressure_events: AtomicUsize,
    pub cpu_overload_events: AtomicUsize,
    pub allocation_failures: AtomicUsize,
    pub optimization_triggers: AtomicUsize,
}

/// Snapshot of resource metrics
#[derive(Debug, Clone)]
pub struct ResourceMetricsSnapshot {
    pub memory_pressure: f64,
    pub cpu_utilization: f64,
    pub available_memory: usize,
    pub total_memory: usize,
    pub active_allocations: usize,
    pub peak_memory_usage: usize,
    pub system_load: f64,
    pub uptime: Duration,
}

impl ResourceMetricsSnapshot {
    /// Format memory values for display
    pub fn format_memory(&self, bytes: usize) -> String {
        if bytes >= 1024 * 1024 * 1024 {
            format!("{:.1} GB", bytes as f64 / (1024.0 * 1024.0 * 1024.0))
        } else if bytes >= 1024 * 1024 {
            format!("{:.1} MB", bytes as f64 / (1024.0 * 1024.0))
        } else if bytes >= 1024 {
            format!("{:.1} KB", bytes as f64 / 1024.0)
        } else {
            format!("{} B", bytes)
        }
    }

    /// Get memory usage percentage
    pub fn memory_usage_percent(&self) -> f64 {
        if self.total_memory == 0 {
            0.0
        } else {
            ((self.total_memory - self.available_memory) as f64 / self.total_memory as f64) * 100.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::laptop_assembly::LaptopConfig;

    #[test]
    fn test_resource_manager_creation() {
        let config = OptimizedConfig::from_laptop_config(LaptopConfig::auto_detect());
        let manager = AdaptiveResourceManager::new(&config);

        assert!(manager.memory_pressure() >= 0.0);
        assert!(manager.memory_pressure() <= 1.0);
        assert!(manager.cpu_utilization() >= 0.0);
        assert!(manager.cpu_utilization() <= 1.0);
    }

    #[test]
    fn test_memory_allocation() {
        let manager = MemoryManager::new(1024); // 1GB

        let allocation = manager.allocate(1024 * 1024, "test").unwrap(); // 1MB
        assert_eq!(allocation.size, 1024 * 1024);
        assert_eq!(manager.active_allocations(), 1);

        drop(allocation);

        // Give a moment for cleanup
        std::thread::sleep(std::time::Duration::from_millis(1));
        assert_eq!(manager.active_allocations(), 0);
    }

    #[test]
    fn test_memory_pressure_calculation() {
        let manager = MemoryManager::new(100); // 100MB

        // No allocations - no pressure
        assert_eq!(manager.memory_pressure(), 0.0);

        // Allocate 50MB - moderate pressure
        let _allocation = manager.allocate(50 * 1024 * 1024, "test").unwrap();
        assert!(manager.memory_pressure() > 0.4);
        assert!(manager.memory_pressure() < 0.6);
    }

    #[test]
    fn test_memory_budget_enforcement() {
        let manager = MemoryManager::new(1); // 1MB total

        // Try to allocate more than budget
        let result = manager.allocate(2 * 1024 * 1024, "too big");
        assert!(result.is_err());
    }

    #[test]
    fn test_optimal_chunk_size() {
        let config = OptimizedConfig::from_laptop_config(LaptopConfig::auto_detect());
        let manager = AdaptiveResourceManager::new(&config);

        let base_size = 1000;
        let optimal = manager.optimal_chunk_size(base_size);

        // Should be within reasonable bounds
        assert!(optimal >= 100);
        assert!(optimal <= base_size * 2);
    }

    #[test]
    fn test_metrics_snapshot() {
        let config = OptimizedConfig::from_laptop_config(LaptopConfig::auto_detect());
        let manager = AdaptiveResourceManager::new(&config);

        let metrics = manager.metrics();

        assert!(metrics.memory_pressure >= 0.0);
        assert!(metrics.memory_pressure <= 1.0);
        assert!(metrics.cpu_utilization >= 0.0);
        assert!(metrics.cpu_utilization <= 1.0);
        assert!(metrics.total_memory > 0);
        assert!(metrics.available_memory <= metrics.total_memory);
    }
}