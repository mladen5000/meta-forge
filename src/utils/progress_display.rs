/// Progress display utilities for metagenomic pipeline
/// 
/// Provides clean terminal progress indicators that overwrite the same line
/// instead of creating new lines for better user experience
use std::io::{self, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

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

        let percentage = if self.total > 0 {
            (self.current as f64 / self.total as f64) * 100.0
        } else {
            0.0
        };

        let filled = (self.width as f64 * percentage / 100.0) as usize;
        let empty = self.width - filled;

        let bar = format!("{}{}",
            "█".repeat(filled),
            "░".repeat(empty)
        );

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
            format!("{:.0}/s", rate)
        };

        // Clear line and move cursor to beginning
        print!("\r\x1b[2K");
        
        if self.total > 0 {
            print!("{} [{}] {:.1}% ({}/{}) {} - {:?}",
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
            print!("{} {} {} {} - {:?}",
                self.message,
                spinner_chars[spinner_idx],
                format_number(self.current),
                rate_str,
                elapsed
            );
        }
        
        io::stdout().flush().unwrap();
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
        println!("{} ✅ {} reads processed in {:?}", 
            self.message, 
            format_number(self.current),
            elapsed
        );
        
        if !message.is_empty() {
            println!("{}", message);
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
            format!("{:.0}/s", rate)
        };

        // Clear line and update
        print!("\r\x1b[2K");
        print!("{} {} {} - {:?}", 
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
        println!("{} ✅ {} processed in {:?}", 
            self.message,
            format_number(self.count),
            elapsed
        );
        
        if !message.is_empty() {
            println!("{}", message);
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
    print!("\r\x1b[2K{}", message);
    io::stdout().flush().unwrap();
}

/// Finish progress line and move to next line
pub fn finish_progress_line(message: &str) {
    print!("\r\x1b[2K{}\n", message);
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
            print!("\r\x1b[2K{}\n", bar);
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