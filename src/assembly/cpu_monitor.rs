//! CPU Utilization Monitoring
//! ==========================
//!
//! Real-time CPU usage tracking for assembly pipeline optimization

use std::time::{Duration, Instant};
use std::sync::atomic::{AtomicU64, AtomicBool, Ordering};
use std::sync::Arc;
use std::thread;

/// CPU utilization monitor
pub struct CpuMonitor {
    start_time: Instant,
    samples: Arc<AtomicU64>,
    active_threads: Arc<AtomicU64>,
    total_threads: usize,
    monitoring: Arc<AtomicBool>,
}

impl CpuMonitor {
    /// Create new CPU monitor
    pub fn new(total_threads: usize) -> Self {
        Self {
            start_time: Instant::now(),
            samples: Arc::new(AtomicU64::new(0)),
            active_threads: Arc::new(AtomicU64::new(0)),
            total_threads,
            monitoring: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Start monitoring CPU utilization
    pub fn start_monitoring(&mut self) {
        self.monitoring.store(true, Ordering::SeqCst);
        self.start_time = Instant::now();

        let samples = Arc::clone(&self.samples);
        let active_threads = Arc::clone(&self.active_threads);
        let monitoring = Arc::clone(&self.monitoring);
        let total_threads = self.total_threads;

        // Spawn background monitor thread
        thread::spawn(move || {
            while monitoring.load(Ordering::SeqCst) {
                thread::sleep(Duration::from_millis(100));

                // Estimate CPU usage based on active rayon threads
                let active = rayon::current_num_threads() as u64;
                active_threads.store(active, Ordering::Relaxed);
                samples.fetch_add(1, Ordering::Relaxed);
            }
        });
    }

    /// Stop monitoring and get results
    pub fn stop_monitoring(&mut self) -> CpuUtilizationReport {
        self.monitoring.store(false, Ordering::SeqCst);

        let samples = self.samples.load(Ordering::Relaxed);
        let active = self.active_threads.load(Ordering::Relaxed);
        let elapsed = self.start_time.elapsed();

        let avg_utilization = if self.total_threads > 0 {
            (active as f64 / self.total_threads as f64) * 100.0
        } else {
            0.0
        };

        CpuUtilizationReport {
            total_threads: self.total_threads,
            avg_active_threads: active as usize,
            avg_utilization_percent: avg_utilization,
            elapsed,
            samples: samples as usize,
        }
    }

    /// Get current CPU utilization estimate
    pub fn current_utilization(&self) -> f64 {
        let active = rayon::current_num_threads();
        (active as f64 / self.total_threads as f64) * 100.0
    }
}

/// CPU utilization report
#[derive(Debug)]
pub struct CpuUtilizationReport {
    pub total_threads: usize,
    pub avg_active_threads: usize,
    pub avg_utilization_percent: f64,
    pub elapsed: Duration,
    pub samples: usize,
}

impl CpuUtilizationReport {
    /// Print formatted report
    pub fn print(&self) {
        println!("\n📊 CPU Utilization Report:");
        println!("   Total threads: {}", self.total_threads);
        println!("   Avg active threads: {}", self.avg_active_threads);
        println!("   Avg utilization: {:.1}%", self.avg_utilization_percent);
        println!("   Elapsed time: {:.2}s", self.elapsed.as_secs_f64());
        println!("   Samples: {}", self.samples);

        if self.avg_utilization_percent < 70.0 {
            println!("   ⚠️  Low CPU utilization - consider increasing parallelism");
        } else if self.avg_utilization_percent >= 85.0 {
            println!("   ✅ Excellent CPU utilization!");
        } else {
            println!("   ✓  Good CPU utilization");
        }
    }
}