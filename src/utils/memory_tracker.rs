//! Memory usage tracking for performance monitoring
//!
//! Provides cross-platform memory tracking capabilities

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

/// Memory tracker that records peak and current usage
#[derive(Clone)]
pub struct MemoryTracker {
    current: Arc<AtomicU64>,
    peak: Arc<AtomicU64>,
}

impl MemoryTracker {
    /// Create a new memory tracker
    pub fn new() -> Self {
        Self {
            current: Arc::new(AtomicU64::new(0)),
            peak: Arc::new(AtomicU64::new(0)),
        }
    }

    /// Record memory allocation
    pub fn allocate(&self, bytes: u64) {
        let new_current = self.current.fetch_add(bytes, Ordering::SeqCst) + bytes;

        // Update peak if needed
        let mut current_peak = self.peak.load(Ordering::SeqCst);
        while new_current > current_peak {
            match self.peak.compare_exchange_weak(
                current_peak,
                new_current,
                Ordering::SeqCst,
                Ordering::SeqCst,
            ) {
                Ok(_) => break,
                Err(x) => current_peak = x,
            }
        }
    }

    /// Record memory deallocation
    pub fn deallocate(&self, bytes: u64) {
        self.current.fetch_sub(bytes, Ordering::SeqCst);
    }

    /// Get current memory usage in bytes
    pub fn current_usage(&self) -> u64 {
        self.current.load(Ordering::SeqCst)
    }

    /// Get peak memory usage in bytes
    pub fn peak_usage(&self) -> u64 {
        self.peak.load(Ordering::SeqCst)
    }

    /// Get system memory usage (process RSS)
    pub fn system_memory_usage() -> Result<u64, String> {
        #[cfg(target_os = "macos")]
        {
            Self::get_macos_memory()
        }
        #[cfg(target_os = "linux")]
        {
            Self::get_linux_memory()
        }
        #[cfg(not(any(target_os = "macos", target_os = "linux")))]
        {
            Ok(0)
        }
    }

    #[cfg(target_os = "macos")]
    fn get_macos_memory() -> Result<u64, String> {
        use std::process::Command;

        let output = Command::new("ps")
            .args(["-o", "rss=", "-p", &std::process::id().to_string()])
            .output()
            .map_err(|e| format!("Failed to get memory usage: {}", e))?;

        let rss_kb = String::from_utf8_lossy(&output.stdout)
            .trim()
            .parse::<u64>()
            .map_err(|e| format!("Failed to parse memory: {}", e))?;

        Ok(rss_kb * 1024) // Convert KB to bytes
    }

    #[cfg(target_os = "linux")]
    fn get_linux_memory() -> Result<u64, String> {
        let status = std::fs::read_to_string("/proc/self/status")
            .map_err(|e| format!("Failed to read /proc/self/status: {}", e))?;

        for line in status.lines() {
            if line.starts_with("VmRSS:") {
                let parts: Vec<&str> = line.split_whitespace().collect();
                if parts.len() >= 2 {
                    let rss_kb = parts[1]
                        .parse::<u64>()
                        .map_err(|e| format!("Failed to parse RSS: {}", e))?;
                    return Ok(rss_kb * 1024); // Convert KB to bytes
                }
            }
        }

        Err("VmRSS not found in /proc/self/status".to_string())
    }

    /// Sample current system memory and update tracker
    pub fn sample_system_memory(&self) {
        if let Ok(mem) = Self::system_memory_usage() {
            self.current.store(mem, Ordering::SeqCst);

            // Update peak
            let mut current_peak = self.peak.load(Ordering::SeqCst);
            while mem > current_peak {
                match self.peak.compare_exchange_weak(
                    current_peak,
                    mem,
                    Ordering::SeqCst,
                    Ordering::SeqCst,
                ) {
                    Ok(_) => break,
                    Err(x) => current_peak = x,
                }
            }
        }
    }

    /// Reset the tracker
    pub fn reset(&self) {
        self.current.store(0, Ordering::SeqCst);
        self.peak.store(0, Ordering::SeqCst);
    }
}

impl Default for MemoryTracker {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_memory_tracker() {
        let tracker = MemoryTracker::new();

        tracker.allocate(1000);
        assert_eq!(tracker.current_usage(), 1000);
        assert_eq!(tracker.peak_usage(), 1000);

        tracker.allocate(500);
        assert_eq!(tracker.current_usage(), 1500);
        assert_eq!(tracker.peak_usage(), 1500);

        tracker.deallocate(800);
        assert_eq!(tracker.current_usage(), 700);
        assert_eq!(tracker.peak_usage(), 1500); // Peak remains
    }

    #[test]
    fn test_system_memory() {
        let result = MemoryTracker::system_memory_usage();
        #[cfg(any(target_os = "macos", target_os = "linux"))]
        {
            assert!(result.is_ok());
            let mem = result.unwrap();
            assert!(mem > 0);
            println!("Current process memory: {} MB", mem / 1024 / 1024);
        }
    }
}
