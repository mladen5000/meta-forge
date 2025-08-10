use anyhow::Result;
use meta_forge::utils::progress_display::{
    finish_progress_line, update_progress_line, MultiProgress, ProgressBar, ProgressCounter,
};
/// Demo of progress display utilities for metagenomic pipeline operations
use std::thread;
use std::time::Duration;

fn main() -> Result<()> {
    println!("🎯 Progress Display Utilities Demo");
    println!("==================================");

    // Demo 1: Progress Bar with Known Total
    println!("\n1️⃣  Progress Bar - Simulating Read Processing");
    println!("--------------------------------------------");
    let total_reads = 50000;
    let mut pb = ProgressBar::new(total_reads, "Processing reads");

    for i in 0..=total_reads {
        pb.update(i);
        if i % 1000 == 0 {
            thread::sleep(Duration::from_millis(10));
        }
    }
    pb.finish_with_message("Read processing completed!");

    // Demo 2: Progress Counter for Unknown Total
    println!("\n2️⃣  Progress Counter - Simulating K-mer Counting");
    println!("-----------------------------------------------");
    let mut counter = ProgressCounter::new("Counting k-mers");

    for i in 0..25000 {
        counter.update(i);
        if i % 500 == 0 {
            thread::sleep(Duration::from_millis(10));
        }
    }
    counter.finish_with_message("K-mer counting finished!");

    // Demo 3: Simple Line Updates
    println!("\n3️⃣  Simple Line Updates - File Operations");
    println!("----------------------------------------");

    let files = [
        "input.fastq",
        "output.fasta",
        "contigs.fa",
        "assembly_graph.gfa",
    ];
    for (i, file) in files.iter().enumerate() {
        update_progress_line(&format!(
            "Processing file {} ({}/{})",
            file,
            i + 1,
            files.len()
        ));
        thread::sleep(Duration::from_millis(500));
    }
    finish_progress_line("All files processed successfully");

    // Demo 4: Multi-line Progress for Complex Operations
    println!("\n4️⃣  Multi-line Progress - Pipeline Operations");
    println!("--------------------------------------------");

    let mut multi = MultiProgress::new();
    let read_line = multi.add_line("Quality filtering: Starting...".to_string());
    let assembly_line = multi.add_line("Assembly: Waiting...".to_string());
    let analysis_line = multi.add_line("Analysis: Waiting...".to_string());

    // Simulate complex pipeline
    for step in 0..20 {
        match step {
            0..=6 => {
                multi.update_line(
                    read_line,
                    format!("Quality filtering: Processing {} reads", step * 1000),
                );
            }
            7..=13 => {
                multi.update_line(read_line, "Quality filtering: ✅ Completed".to_string());
                multi.update_line(
                    assembly_line,
                    format!("Assembly: Building graph (k={})", 21 + (step - 7) * 2),
                );
            }
            14..=19 => {
                multi.update_line(assembly_line, "Assembly: ✅ Completed".to_string());
                multi.update_line(
                    analysis_line,
                    format!("Analysis: Analyzing contigs ({}/100)", (step - 13) * 20),
                );
            }
            _ => {}
        }
        thread::sleep(Duration::from_millis(200));
    }
    multi.update_line(analysis_line, "Analysis: ✅ Completed".to_string());
    thread::sleep(Duration::from_millis(500));
    multi.finish();

    // Demo 5: Different Progress Bar Styles
    println!("\n5️⃣  Progress Bar Variations");
    println!("---------------------------");

    // Fast processing simulation
    println!("Fast operation (high throughput):");
    let mut fast_pb = ProgressBar::new(100000, "High-speed processing");
    for i in 0..=100000 {
        fast_pb.update(i);
        if i % 5000 == 0 {
            thread::sleep(Duration::from_millis(10));
        }
    }
    fast_pb.finish();

    // Slow processing simulation
    println!("\nSlow operation (complex analysis):");
    let mut slow_pb = ProgressBar::new(20, "Complex analysis");
    for i in 0..=20 {
        slow_pb.update(i);
        thread::sleep(Duration::from_millis(150));
    }
    slow_pb.finish();

    // Demo 6: Spinner for Indeterminate Progress
    println!("\n6️⃣  Indeterminate Progress - Database Operations");
    println!("-----------------------------------------------");
    let mut spinner = ProgressBar::new(0, "Connecting to database"); // 0 = indeterminate

    for i in 0..50 {
        spinner.update(i);
        thread::sleep(Duration::from_millis(100));
    }
    spinner.finish_with_message("Database connection established!");

    println!("\n✨ Progress Display Demo Complete!");
    println!("=================================");
    println!("All progress utilities demonstrated:");
    println!("✅ Progress bars with known totals");
    println!("✅ Progress counters for unknown totals");
    println!("✅ Simple line updates");
    println!("✅ Multi-line progress displays");
    println!("✅ Different processing speeds");
    println!("✅ Indeterminate progress (spinners)");
    println!("\n🧬 Ready for integration with metagenomic pipelines!");

    Ok(())
}
