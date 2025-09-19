/// Progress display utilities for metagenomic pipeline
///
/// Provides clean terminal progress indicators that overwrite the same line
/// instead of creating new lines for better user experience
use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

/// Assembly-specific progress information
#[derive(Debug, Clone)]
pub struct AssemblyProgress {
    pub total_reads: usize,
    pub processed_reads: usize,
    pub kmers_extracted: usize,
    pub nodes_created: usize,
    pub edges_created: usize,
    pub memory_used_mb: f64,
    pub stage: AssemblyStage,
    pub stage_progress: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AssemblyStage {
    Initialization,
    KmerExtraction,
    GraphConstruction,
    ErrorCorrection,
    TransitiveReduction,
    ContigGeneration,
    Finalization,
    Complete,
}

impl AssemblyStage {
    pub fn description(&self) -> &str {
        match self {
            AssemblyStage::Initialization => "Initializing assembly parameters",
            AssemblyStage::KmerExtraction => "Extracting k-mers from sequences",
            AssemblyStage::GraphConstruction => "Building De Bruijn graph",
            AssemblyStage::ErrorCorrection => "Correcting sequencing errors",
            AssemblyStage::TransitiveReduction => "Simplifying graph structure",
            AssemblyStage::ContigGeneration => "Generating contiguous sequences",
            AssemblyStage::Finalization => "Finalizing assembly results",
            AssemblyStage::Complete => "Assembly completed",
        }
    }

    pub fn icon(&self) -> &str {
        match self {
            AssemblyStage::Initialization => "🔧",
            AssemblyStage::KmerExtraction => "✂️",
            AssemblyStage::GraphConstruction => "🌐",
            AssemblyStage::ErrorCorrection => "🔍",
            AssemblyStage::TransitiveReduction => "🔗",
            AssemblyStage::ContigGeneration => "📝",
            AssemblyStage::Finalization => "📊",
            AssemblyStage::Complete => "✅",
        }
    }
}

/// Progress bar for genomic processing operations
pub struct ProgressBar {
    total: u64,
    current: u64,
    width: usize,
    message: String,
    start_time: Instant,
    last_update: Instant,
    update_interval: Duration,
    finished: Arc<AtomicBool>,
    assembly_info: Option<AssemblyProgress>,
}

impl ProgressBar {
    /// Create a new progress bar
    pub fn new(total: u64, message: &str) -> Self {
        Self {
            total,
            current: 0,
            width: 50,
            message: message.to_string(),
            start_time: Instant::now(),
            last_update: Instant::now(),
            update_interval: Duration::from_millis(100), // Update every 100ms
            finished: Arc::new(AtomicBool::new(false)),
            assembly_info: None,
        }
    }

    /// Create a new assembly-specific progress bar
    pub fn new_assembly(total: u64, message: &str, total_reads: usize) -> Self {
        Self {
            total,
            current: 0,
            width: 50,
            message: message.to_string(),
            start_time: Instant::now(),
            last_update: Instant::now(),
            update_interval: Duration::from_millis(100),
            finished: Arc::new(AtomicBool::new(false)),
            assembly_info: Some(AssemblyProgress {
                total_reads,
                processed_reads: 0,
                kmers_extracted: 0,
                nodes_created: 0,
                edges_created: 0,
                memory_used_mb: 0.0,
                stage: AssemblyStage::Initialization,
                stage_progress: 0.0,
            }),
        }
    }

    /// Update progress with current count
    pub fn update(&mut self, current: u64) {
        self.current = current;

        // Only update display if enough time has passed
        let now = Instant::now();
        if now.duration_since(self.last_update) >= self.update_interval {
            self.display();
            self.last_update = now;
        }
    }

    /// Increment progress by 1
    pub fn tick(&mut self) {
        self.update(self.current + 1);
    }

    /// Increment progress by n
    pub fn tick_by(&mut self, n: u64) {
        self.update(self.current + n);
    }

    /// Display the progress bar  
    fn display(&self) {
        if self.finished.load(Ordering::Relaxed) {
            return;
        }

        // Use enhanced assembly display if available
        if self.assembly_info.is_some() {
            self.display_assembly_progress();
            return;
        }

        let percentage = if self.total > 0 {
            (self.current as f64 / self.total as f64) * 100.0
        } else {
            0.0
        };

        let filled = (self.width as f64 * percentage / 100.0) as usize;
        let empty = self.width - filled;

        let bar = format!("{}{}", "█".repeat(filled), "░".repeat(empty));

        let elapsed = self.start_time.elapsed();
        let rate = if elapsed.as_secs() > 0 {
            self.current as f64 / elapsed.as_secs_f64()
        } else {
            0.0
        };

        // Format rate based on the number
        let rate_str = if rate >= 1000.0 {
            format!("{:.1}K/s", rate / 1000.0)
        } else {
            format!("{rate:.0}/s")
        };

        // Clear line and move cursor to beginning
        print!("\r\x1b[2K");

        if self.total > 0 {
            print!(
                "{} [{}] {:.1}% ({}/{}) {} - {:?}",
                self.message,
                bar,
                percentage,
                format_number(self.current),
                format_number(self.total),
                rate_str,
                elapsed
            );
        } else {
            // Spinner for indeterminate progress
            let spinner_chars = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏'];
            let spinner_idx = (self.current % spinner_chars.len() as u64) as usize;
            print!(
                "{} {} {} {} - {:?}",
                self.message,
                spinner_chars[spinner_idx],
                format_number(self.current),
                rate_str,
                elapsed
            );
        }

        io::stdout().flush().unwrap();
    }

    /// Update assembly-specific progress information
    pub fn update_assembly_info(
        &mut self,
        stage: AssemblyStage,
        processed_reads: usize,
        kmers: usize,
        nodes: usize,
        edges: usize,
        memory_mb: f64,
    ) {
        if let Some(info) = &mut self.assembly_info {
            info.stage = stage;
            info.processed_reads = processed_reads;
            info.kmers_extracted = kmers;
            info.nodes_created = nodes;
            info.edges_created = edges;
            info.memory_used_mb = memory_mb;
            info.stage_progress = if info.total_reads > 0 {
                (processed_reads as f64 / info.total_reads as f64) * 100.0
            } else {
                0.0
            };
        }
    }

    /// Display enhanced assembly progress with detailed metrics
    fn display_assembly_progress(&self) {
        if let Some(info) = &self.assembly_info {
            let stage = &info.stage;
            let elapsed = self.start_time.elapsed();

            // Clear line and show comprehensive progress
            print!("\r\x1b[2K");

            // Stage header with icon and description
            println!(
                "\n{} {} - {:.1}% complete",
                stage.icon(),
                stage.description(),
                info.stage_progress
            );

            // Main progress bar
            let percentage = if self.total > 0 {
                (self.current as f64 / self.total as f64) * 100.0
            } else {
                info.stage_progress
            };
            let filled = (self.width as f64 * percentage / 100.0) as usize;
            let empty = self.width - filled;
            let bar = format!("{}{}", "█".repeat(filled), "░".repeat(empty));

            // Detailed metrics line
            print!(
                "├─ [{}] {:.1}% | Reads: {}/{} | K-mers: {} | Nodes: {} | Edges: {}",
                bar,
                percentage,
                format_number(info.processed_reads as u64),
                format_number(info.total_reads as u64),
                format_number(info.kmers_extracted as u64),
                format_number(info.nodes_created as u64),
                format_number(info.edges_created as u64)
            );

            // Memory and performance info
            let memory_bar = if info.memory_used_mb > 4000.0 {
                "🔴" // Red for high memory usage
            } else if info.memory_used_mb > 2000.0 {
                "🟡" // Yellow for medium memory usage
            } else {
                "🟢" // Green for low memory usage
            };

            let rate = if elapsed.as_secs() > 0 {
                info.processed_reads as f64 / elapsed.as_secs_f64()
            } else {
                0.0
            };

            println!(
                "\n├─ {} Memory: {:.1}MB | Rate: {:.0} reads/sec | Elapsed: {:?}",
                memory_bar, info.memory_used_mb, rate, elapsed
            );

            // ETA if available
            if let Some(eta) = self.eta() {
                println!("└─ ⏱️  ETA: {:?}", eta);
            } else {
                println!("└─ ⏱️  ETA: calculating...");
            }
        }
    }

    /// Get estimated time remaining
    pub fn eta(&self) -> Option<Duration> {
        if self.current == 0 || self.total == 0 {
            return None;
        }

        let elapsed = self.start_time.elapsed();
        let rate = self.current as f64 / elapsed.as_secs_f64();
        let remaining = self.total - self.current;

        if rate > 0.0 {
            Some(Duration::from_secs_f64(remaining as f64 / rate))
        } else {
            None
        }
    }

    /// Set a custom message while preserving progress
    pub fn set_message(&mut self, message: &str) {
        self.message = message.to_string();
    }

    /// Get current progress percentage
    pub fn percentage(&self) -> f64 {
        if self.total > 0 {
            (self.current as f64 / self.total as f64) * 100.0
        } else {
            0.0
        }
    }

    /// Finish the progress bar and move to next line
    pub fn finish(&mut self) {
        self.finished.store(true, Ordering::Relaxed);

        // Display final state
        self.display();

        // Move to next line
        println!();
        io::stdout().flush().unwrap();
    }

    /// Finish with a custom message
    pub fn finish_with_message(&mut self, message: &str) {
        self.finished.store(true, Ordering::Relaxed);

        // Clear line and print final message
        print!("\r\x1b[2K");
        let elapsed = self.start_time.elapsed();
        println!(
            "{} ✅ {} reads processed in {:?}",
            self.message,
            format_number(self.current),
            elapsed
        );

        if !message.is_empty() {
            println!("{message}");
        }

        io::stdout().flush().unwrap();
    }
}

/// Simple progress counter that overwrites the same line
pub struct ProgressCounter {
    message: String,
    count: u64,
    last_update: Instant,
    update_interval: Duration,
    start_time: Instant,
}

impl ProgressCounter {
    /// Create a new progress counter
    pub fn new(message: &str) -> Self {
        Self {
            message: message.to_string(),
            count: 0,
            last_update: Instant::now(),
            update_interval: Duration::from_millis(250), // Update every 250ms
            start_time: Instant::now(),
        }
    }

    /// Update the counter
    pub fn update(&mut self, count: u64) {
        self.count = count;

        let now = Instant::now();
        if now.duration_since(self.last_update) >= self.update_interval {
            self.display();
            self.last_update = now;
        }
    }

    /// Increment counter by 1
    pub fn tick(&mut self) {
        self.update(self.count + 1);
    }

    /// Increment counter by n
    pub fn tick_by(&mut self, n: u64) {
        self.update(self.count + n);
    }

    /// Display the current count
    fn display(&self) {
        let elapsed = self.start_time.elapsed();
        let rate = if elapsed.as_secs() > 0 {
            self.count as f64 / elapsed.as_secs_f64()
        } else {
            0.0
        };

        let rate_str = if rate >= 1000.0 {
            format!("{:.1}K/s", rate / 1000.0)
        } else {
            format!("{rate:.0}/s")
        };

        // Clear line and update
        print!("\r\x1b[2K");
        print!(
            "{} {} {} - {:?}",
            self.message,
            format_number(self.count),
            rate_str,
            elapsed
        );

        io::stdout().flush().unwrap();
    }

    /// Finish and move to next line
    pub fn finish(&mut self) {
        self.display();
        println!();
        io::stdout().flush().unwrap();
    }

    /// Finish with custom message
    pub fn finish_with_message(&mut self, message: &str) {
        print!("\r\x1b[2K");
        let elapsed = self.start_time.elapsed();
        println!(
            "{} ✅ {} processed in {:?}",
            self.message,
            format_number(self.count),
            elapsed
        );

        if !message.is_empty() {
            println!("{message}");
        }

        io::stdout().flush().unwrap();
    }
}

/// Format large numbers with appropriate suffixes
fn format_number(num: u64) -> String {
    if num >= 1_000_000_000 {
        format!("{:.1}G", num as f64 / 1_000_000_000.0)
    } else if num >= 1_000_000 {
        format!("{:.1}M", num as f64 / 1_000_000.0)
    } else if num >= 1_000 {
        format!("{:.1}K", num as f64 / 1_000.0)
    } else {
        num.to_string()
    }
}

/// Convenience function for quick progress updates
/// Uses carriage return to overwrite the same line
pub fn update_progress_line(message: &str) {
    print!("\r\x1b[2K{message}");
    io::stdout().flush().unwrap();
}

/// Finish progress line and move to next line
pub fn finish_progress_line(message: &str) {
    print!("\r\x1b[2K{message}\n");
    io::stdout().flush().unwrap();
}

/// Multi-line progress display for complex operations
pub struct MultiProgress {
    bars: Vec<String>,
    active_bar: usize,
}

impl MultiProgress {
    /// Create a new multi-line progress display
    pub fn new() -> Self {
        Self {
            bars: Vec::new(),
            active_bar: 0,
        }
    }

    /// Add a progress line
    pub fn add_line(&mut self, message: String) -> usize {
        self.bars.push(message);
        self.bars.len() - 1
    }

    /// Update a specific line
    pub fn update_line(&mut self, index: usize, message: String) {
        if index < self.bars.len() {
            self.bars[index] = message;
            self.display();
        }
    }

    /// Display all progress lines
    fn display(&self) {
        // Move cursor up to overwrite previous lines
        if !self.bars.is_empty() {
            print!("\x1b[{}A", self.bars.len());
        }

        // Clear and redraw each line
        for bar in &self.bars {
            print!("\r\x1b[2K{bar}\n");
        }

        io::stdout().flush().unwrap();
    }

    /// Finish all progress lines
    pub fn finish(&mut self) {
        println!();
        io::stdout().flush().unwrap();
    }
}

impl Default for MultiProgress {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::thread;
    use std::time::Duration;

    #[test]
    fn test_progress_bar() {
        let mut pb = ProgressBar::new(100, "Testing");

        for i in 0..=100 {
            pb.update(i);
            thread::sleep(Duration::from_millis(10));
        }

        pb.finish_with_message("Test completed successfully");
    }

    #[test]
    fn test_progress_counter() {
        let mut counter = ProgressCounter::new("Processing reads");

        for i in 0..1000 {
            counter.tick();
            if i % 100 == 0 {
                thread::sleep(Duration::from_millis(1));
            }
        }

        counter.finish_with_message("All reads processed");
    }

    #[test]
    fn test_format_number() {
        assert_eq!(format_number(123), "123");
        assert_eq!(format_number(1234), "1.2K");
        assert_eq!(format_number(1234567), "1.2M");
        assert_eq!(format_number(1234567890), "1.2G");
    }
}
