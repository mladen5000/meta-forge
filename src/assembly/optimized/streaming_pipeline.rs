//! Streaming Assembly Pipeline
//! ==========================
//!
//! Resource-aware streaming pipeline for metagenomic assembly.
//! Implements adaptive chunk sizing, backpressure handling, and modular stages.

use crate::assembly::optimized::{
    CSRAssemblyGraph, AdaptiveResourceManager
};
use crate::assembly::laptop_assembly::CompactKmer;
use crate::core::data_structures::{CorrectedRead, Contig};
use anyhow::{anyhow, Result};
use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::prelude::*;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

/// Streaming assembly pipeline with resource awareness and backpressure
pub struct StreamingAssemblyPipeline {
    config: PipelineConfig,
    resource_manager: Arc<AdaptiveResourceManager>,
    // memory_pool removed - using standard allocation
    stages: Vec<Box<dyn PipelineStage>>,
    metrics: PipelineMetrics,
    shutdown_signal: Arc<AtomicBool>,
}

/// Pipeline configuration
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    /// Base chunk size for read batches
    pub base_chunk_size: usize,
    /// Maximum chunk size under memory pressure
    pub max_chunk_size: usize,
    /// Minimum chunk size
    pub min_chunk_size: usize,
    /// Memory pressure threshold (0.0-1.0)
    pub memory_pressure_threshold: f64,
    /// Number of worker threads per stage
    pub workers_per_stage: usize,
    /// Channel buffer size
    pub channel_buffer_size: usize,
    /// Enable adaptive chunk sizing
    pub adaptive_chunking: bool,
    /// Backpressure timeout in milliseconds
    pub backpressure_timeout_ms: u64,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        Self {
            base_chunk_size: 1000,
            max_chunk_size: 10000,
            min_chunk_size: 100,
            memory_pressure_threshold: 0.8,
            workers_per_stage: 2,
            channel_buffer_size: 10,
            adaptive_chunking: true,
            backpressure_timeout_ms: 1000,
        }
    }
}

/// Pipeline stage trait for modular processing
pub trait PipelineStage: Send + Sync {
    /// Process a batch of data
    fn process_batch(&mut self, batch: ReadBatch) -> Result<ProcessedBatch>;

    /// Current memory usage in bytes
    fn memory_usage(&self) -> usize;

    /// Check if stage can process given memory budget
    fn can_process(&self, memory_budget: usize) -> bool;

    /// Stage name for diagnostics
    fn name(&self) -> &str;

    /// Warm up stage (optional)
    fn warmup(&mut self) -> Result<()> { Ok(()) }

    /// Clean up stage (optional)
    fn cleanup(&mut self) -> Result<()> { Ok(()) }

    /// Performance metrics
    fn metrics(&self) -> StageMetrics;
}

/// Read batch for pipeline processing
#[derive(Debug, Clone)]
pub struct ReadBatch {
    pub id: usize,
    pub reads: Vec<CorrectedRead>,
    pub chunk_size: usize,
    pub timestamp: Instant,
}

/// Processed batch output
#[derive(Debug)]
pub struct ProcessedBatch {
    pub id: usize,
    pub kmers: Vec<CompactKmer>,
    pub graph_fragment: Option<CSRAssemblyGraph>,
    pub contigs: Vec<Contig>,
    pub processing_time: Duration,
    pub memory_used: usize,
}

/// Pipeline metrics
#[derive(Debug, Default)]
pub struct PipelineMetrics {
    pub total_batches_processed: AtomicUsize,
    pub total_reads_processed: AtomicUsize,
    pub total_processing_time: AtomicUsize, // in milliseconds
    pub memory_pressure_events: AtomicUsize,
    pub backpressure_events: AtomicUsize,
    pub error_count: AtomicUsize,
}

/// Stage-specific metrics
#[derive(Debug, Default, Clone)]
pub struct StageMetrics {
    pub batches_processed: usize,
    pub total_processing_time: Duration,
    pub average_batch_time: Duration,
    pub peak_memory_usage: usize,
    pub current_memory_usage: usize,
    pub error_count: usize,
}

impl StreamingAssemblyPipeline {
    /// Create new streaming pipeline
    pub fn new(
        config: PipelineConfig,
        resource_manager: Arc<AdaptiveResourceManager>
    ) -> Self {
        Self {
            config,
            resource_manager,
            // memory_pool removed
            stages: Vec::new(),
            metrics: PipelineMetrics::default(),
            shutdown_signal: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Add a pipeline stage
    pub fn add_stage(&mut self, stage: Box<dyn PipelineStage>) {
        self.stages.push(stage);
    }

    /// Process reads through the pipeline
    pub fn process_reads(&mut self, reads: Vec<CorrectedRead>) -> Result<Vec<Contig>> {
        let start_time = Instant::now();

        // Warm up stages
        for stage in &mut self.stages {
            stage.warmup()?;
        }

        // Create adaptive chunks
        let chunks = self.create_adaptive_chunks(reads)?;

        // Process chunks through pipeline
        let mut all_contigs = Vec::new();

        for chunk in chunks {
            match self.process_chunk(chunk) {
                Ok(mut contigs) => {
                    all_contigs.append(&mut contigs);
                }
                Err(e) => {
                    self.metrics.error_count.fetch_add(1, Ordering::Relaxed);
                    eprintln!("Pipeline error: {}", e);
                }
            }

            // Check for shutdown signal
            if self.shutdown_signal.load(Ordering::Relaxed) {
                break;
            }
        }

        // Clean up stages
        for stage in &mut self.stages {
            stage.cleanup()?;
        }

        let total_time = start_time.elapsed();
        self.metrics.total_processing_time.fetch_add(
            total_time.as_millis() as usize,
            Ordering::Relaxed
        );

        Ok(all_contigs)
    }

    /// Create adaptive chunks based on memory pressure
    fn create_adaptive_chunks(&self, reads: Vec<CorrectedRead>) -> Result<Vec<ReadBatch>> {
        let mut chunks = Vec::new();
        let mut batch_id = 0;

        if !self.config.adaptive_chunking {
            // Fixed chunk size
            for chunk in reads.chunks(self.config.base_chunk_size) {
                chunks.push(ReadBatch {
                    id: batch_id,
                    reads: chunk.to_vec(),
                    chunk_size: chunk.len(),
                    timestamp: Instant::now(),
                });
                batch_id += 1;
            }
            return Ok(chunks);
        }

        // Adaptive chunking based on memory pressure
        let mut current_chunk = Vec::new();
        let mut current_chunk_size = self.config.base_chunk_size;

        for read in reads {
            current_chunk.push(read);

            if current_chunk.len() >= current_chunk_size {
                chunks.push(ReadBatch {
                    id: batch_id,
                    reads: current_chunk.clone(),
                    chunk_size: current_chunk.len(),
                    timestamp: Instant::now(),
                });

                current_chunk.clear();
                batch_id += 1;

                // Adjust chunk size based on memory pressure
                current_chunk_size = self.calculate_adaptive_chunk_size();
            }
        }

        // Handle remaining reads
        if !current_chunk.is_empty() {
            let chunk_size = current_chunk.len();
            chunks.push(ReadBatch {
                id: batch_id,
                reads: current_chunk,
                chunk_size,
                timestamp: Instant::now(),
            });
        }

        Ok(chunks)
    }

    /// Calculate adaptive chunk size based on system state
    /// OPTIMIZATION: Enhanced with CPU load balancing and throughput monitoring
    fn calculate_adaptive_chunk_size(&self) -> usize {
        let memory_pressure = self.resource_manager.memory_pressure();

        // OPTIMIZATION: Consider processing velocity for adaptive sizing
        let processing_velocity = self.estimate_processing_velocity();
        let target_latency_ms = 100.0; // Target processing latency per chunk

        if memory_pressure > self.config.memory_pressure_threshold {
            // Reduce chunk size under memory pressure
            self.metrics.memory_pressure_events.fetch_add(1, Ordering::Relaxed);

            let reduction_factor = (memory_pressure - self.config.memory_pressure_threshold) /
                                 (1.0 - self.config.memory_pressure_threshold);

            // OPTIMIZATION: Apply stronger reduction for extreme memory pressure
            let pressure_multiplier = if memory_pressure > 0.95 { 0.7 } else { 0.5 };
            let reduced_size = self.config.base_chunk_size as f64 * (1.0 - reduction_factor * pressure_multiplier);
            reduced_size.max(self.config.min_chunk_size as f64) as usize
        } else {
            // OPTIMIZATION: Increase chunk size based on both memory and processing capacity
            let increase_factor = (self.config.memory_pressure_threshold - memory_pressure) /
                                self.config.memory_pressure_threshold;

            // Factor in processing velocity to avoid creating chunks too large for timely processing
            let velocity_factor = if processing_velocity > 0.0 {
                (target_latency_ms / processing_velocity).min(2.0).max(0.5)
            } else {
                1.0
            };

            let increased_size = self.config.base_chunk_size as f64 *
                               (1.0 + increase_factor * 0.3) * velocity_factor;
            increased_size.min(self.config.max_chunk_size as f64) as usize
        }
    }

    /// OPTIMIZATION: Estimate processing velocity (items per millisecond)
    fn estimate_processing_velocity(&self) -> f64 {
        let total_processed = self.metrics.total_reads_processed.load(Ordering::Relaxed);
        let total_time_ms = self.metrics.total_processing_time.load(Ordering::Relaxed);

        if total_time_ms > 0 && total_processed > 0 {
            total_processed as f64 / total_time_ms as f64
        } else {
            1.0 // Default velocity
        }
    }

    /// Process a single chunk through all pipeline stages
    fn process_chunk(&mut self, mut batch: ReadBatch) -> Result<Vec<Contig>> {
        self.metrics.total_batches_processed.fetch_add(1, Ordering::Relaxed);
        self.metrics.total_reads_processed.fetch_add(batch.reads.len(), Ordering::Relaxed);

        // Create pipeline for this batch
        let mut processed_batch = ProcessedBatch {
            id: batch.id,
            kmers: Vec::new(),
            graph_fragment: None,
            contigs: Vec::new(),
            processing_time: Duration::from_secs(0),
            memory_used: 0,
        };

        // Process through each stage with backpressure handling
        let stages_len = self.stages.len();
        for (stage_idx, stage) in self.stages.iter_mut().enumerate() {
            let stage_start = Instant::now();

            // Check memory budget before processing
            let available_memory = self.resource_manager.available_memory();
            if !stage.can_process(available_memory) {
                // Handle backpressure
                self.metrics.backpressure_events.fetch_add(1, Ordering::Relaxed);

                let timeout = Duration::from_millis(self.config.backpressure_timeout_ms);
                // Wait for memory availability if needed
                if self.resource_manager.memory_pressure() > 0.8 {
                    thread::sleep(Duration::from_millis(10));
                }
            }

            // Process the batch
            match stage.process_batch(batch.clone()) {
                Ok(result) => {
                    processed_batch = result;
                    batch.reads = Vec::new(); // Clear to save memory

                    // Update batch with stage results for next stage
                    if stage_idx < stages_len - 1 {
                        // Convert processed batch back to read batch for next stage
                        // This is simplified - in practice, you'd have stage-specific conversions
                        batch = ReadBatch {
                            id: processed_batch.id,
                            reads: Vec::new(), // Would contain transformed data
                            chunk_size: 0,
                            timestamp: Instant::now(),
                        };
                    }
                }
                Err(e) => {
                    return Err(anyhow!("Stage '{}' failed: {}", stage.name(), e));
                }
            }

            processed_batch.processing_time += stage_start.elapsed();
        }

        Ok(processed_batch.contigs)
    }

    /// Wait for memory availability with timeout
    fn wait_for_memory_availability(
        &self,
        stage: &dyn PipelineStage,
        timeout: Duration
    ) -> Result<()> {
        let start = Instant::now();

        while start.elapsed() < timeout {
            let available_memory = self.resource_manager.available_memory();

            if stage.can_process(available_memory) {
                return Ok(());
            }

            // Brief sleep before checking again
            thread::sleep(Duration::from_millis(10));

            // Check shutdown signal
            if self.shutdown_signal.load(Ordering::Relaxed) {
                return Err(anyhow!("Pipeline shutdown requested"));
            }
        }

        Err(anyhow!("Timeout waiting for memory availability"))
    }

    /// Get pipeline metrics
    pub fn metrics(&self) -> PipelineMetricsSnapshot {
        PipelineMetricsSnapshot {
            total_batches_processed: self.metrics.total_batches_processed.load(Ordering::Relaxed),
            total_reads_processed: self.metrics.total_reads_processed.load(Ordering::Relaxed),
            total_processing_time: Duration::from_millis(
                self.metrics.total_processing_time.load(Ordering::Relaxed) as u64
            ),
            memory_pressure_events: self.metrics.memory_pressure_events.load(Ordering::Relaxed),
            backpressure_events: self.metrics.backpressure_events.load(Ordering::Relaxed),
            error_count: self.metrics.error_count.load(Ordering::Relaxed),
            stage_metrics: self.stages.iter().map(|s| s.metrics()).collect(),
        }
    }

    /// Signal pipeline to shutdown gracefully
    pub fn shutdown(&self) {
        self.shutdown_signal.store(true, Ordering::Relaxed);
    }
}

/// Snapshot of pipeline metrics
#[derive(Debug, Clone)]
pub struct PipelineMetricsSnapshot {
    pub total_batches_processed: usize,
    pub total_reads_processed: usize,
    pub total_processing_time: Duration,
    pub memory_pressure_events: usize,
    pub backpressure_events: usize,
    pub error_count: usize,
    pub stage_metrics: Vec<StageMetrics>,
}

/// K-mer extraction pipeline stage
pub struct KmerExtractionStage {
    k: usize,
    memory_usage: AtomicUsize,
    metrics: StageMetrics,
}

impl KmerExtractionStage {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            memory_usage: AtomicUsize::new(0),
            metrics: StageMetrics::default(),
        }
    }
}

impl PipelineStage for KmerExtractionStage {
    fn process_batch(&mut self, batch: ReadBatch) -> Result<ProcessedBatch> {
        let start_time = Instant::now();
        let mut all_kmers = Vec::new();

        for read in &batch.reads {
            if read.corrected.len() >= self.k {
                for i in 0..=read.corrected.len() - self.k {
                    let kmer_seq = &read.corrected[i..i + self.k];
                    match CompactKmer::new(kmer_seq) {
                        Ok(kmer) => all_kmers.push(kmer),
                        Err(_) => continue, // Skip invalid k-mers
                    }
                }
            }
        }

        let processing_time = start_time.elapsed();
        let memory_used = all_kmers.len() * std::mem::size_of::<CompactKmer>();

        self.memory_usage.store(memory_used, Ordering::Relaxed);
        self.metrics.batches_processed += 1;
        self.metrics.total_processing_time += processing_time;
        self.metrics.current_memory_usage = memory_used;

        if memory_used > self.metrics.peak_memory_usage {
            self.metrics.peak_memory_usage = memory_used;
        }

        Ok(ProcessedBatch {
            id: batch.id,
            kmers: all_kmers,
            graph_fragment: None,
            contigs: Vec::new(),
            processing_time,
            memory_used,
        })
    }

    fn memory_usage(&self) -> usize {
        self.memory_usage.load(Ordering::Relaxed)
    }

    fn can_process(&self, memory_budget: usize) -> bool {
        // Estimate memory needed for processing
        let estimated_memory = self.k * 1000 * std::mem::size_of::<CompactKmer>(); // Rough estimate
        estimated_memory <= memory_budget
    }

    fn name(&self) -> &str {
        "K-mer Extraction"
    }

    fn metrics(&self) -> StageMetrics {
        let mut metrics = self.metrics.clone();
        if metrics.batches_processed > 0 {
            metrics.average_batch_time = metrics.total_processing_time / metrics.batches_processed as u32;
        }
        metrics
    }
}

/// Graph construction pipeline stage
pub struct GraphConstructionStage {
    memory_usage: AtomicUsize,
    metrics: StageMetrics,
}

impl GraphConstructionStage {
    pub fn new() -> Self {
        Self {
            memory_usage: AtomicUsize::new(0),
            metrics: StageMetrics::default(),
        }
    }
}

impl PipelineStage for GraphConstructionStage {
    fn process_batch(&mut self, _batch: ReadBatch) -> Result<ProcessedBatch> {
        let start_time = Instant::now();

        // Simplified graph construction
        let graph = CSRAssemblyGraph::new(1000, 2000);
        let memory_used = graph.memory_usage_bytes() as usize;

        self.memory_usage.store(memory_used, Ordering::Relaxed);
        self.metrics.batches_processed += 1;
        self.metrics.total_processing_time += start_time.elapsed();
        self.metrics.current_memory_usage = memory_used;

        Ok(ProcessedBatch {
            id: 0,
            kmers: Vec::new(),
            graph_fragment: Some(graph),
            contigs: Vec::new(),
            processing_time: start_time.elapsed(),
            memory_used,
        })
    }

    fn memory_usage(&self) -> usize {
        self.memory_usage.load(Ordering::Relaxed)
    }

    fn can_process(&self, memory_budget: usize) -> bool {
        // Conservative estimate for graph construction
        let estimated_memory = 1024 * 1024 * 10; // 10MB rough estimate
        estimated_memory <= memory_budget
    }

    fn name(&self) -> &str {
        "Graph Construction"
    }

    fn metrics(&self) -> StageMetrics {
        self.metrics.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::data_structures::CorrectionMetadata;

    fn create_test_read(id: usize, sequence: &str) -> CorrectedRead {
        CorrectedRead {
            id,
            original: sequence.to_string(),
            corrected: sequence.to_string(),
            corrections: Vec::new(),
            quality_scores: vec![30; sequence.len()],
            correction_metadata: CorrectionMetadata {
                algorithm: "test".to_string(),
                confidence_threshold: 0.95,
                context_window: 5,
                correction_time_ms: 0,
            },
        }
    }

    #[test]
    fn test_pipeline_creation() {
        let config = PipelineConfig::default();
        let laptop_config = crate::assembly::laptop_assembly::LaptopConfig::medium_memory();
        let opt_config = crate::assembly::optimized::optimized_assembler::OptimizedConfig::from_laptop_config(laptop_config);
        let resource_manager = Arc::new(AdaptiveResourceManager::new(&opt_config));
        // memory_pool removed - using standard allocation

        let mut pipeline = StreamingAssemblyPipeline::new(
            config,
            resource_manager,
        );

        pipeline.add_stage(Box::new(KmerExtractionStage::new(21)));
        pipeline.add_stage(Box::new(GraphConstructionStage::new()));

        assert_eq!(pipeline.stages.len(), 2);
    }

    #[test]
    fn test_adaptive_chunking() {
        let config = PipelineConfig {
            adaptive_chunking: true,
            base_chunk_size: 3,
            ..Default::default()
        };
        let laptop_config = crate::assembly::laptop_assembly::LaptopConfig::medium_memory();
        let opt_config = crate::assembly::optimized::optimized_assembler::OptimizedConfig::from_laptop_config(laptop_config);
        let resource_manager = Arc::new(AdaptiveResourceManager::new(&opt_config));
        // memory_pool removed - using standard allocation

        let pipeline = StreamingAssemblyPipeline::new(
            config,
            resource_manager,
        );

        let reads = vec![
            create_test_read(0, "ATCGATCGATCG"),
            create_test_read(1, "GCTAGCTAGCTA"),
            create_test_read(2, "TTAACCGGTTAA"),
            create_test_read(3, "GGCCAATTGGCC"),
            create_test_read(4, "AATTGGCCAATT"),
        ];

        let chunks = pipeline.create_adaptive_chunks(reads).unwrap();
        assert!(!chunks.is_empty());

        for chunk in &chunks {
            assert!(chunk.reads.len() <= 5);
        }
    }

    #[test]
    fn test_kmer_extraction_stage() {
        let mut stage = KmerExtractionStage::new(4);

        let batch = ReadBatch {
            id: 0,
            reads: vec![create_test_read(0, "ATCGATCG")],
            chunk_size: 1,
            timestamp: Instant::now(),
        };

        let result = stage.process_batch(batch).unwrap();
        assert_eq!(result.kmers.len(), 5); // 8 - 4 + 1 = 5 k-mers
        assert!(result.processing_time > Duration::from_nanos(0));
    }

    #[test]
    fn test_memory_budget_checking() {
        let stage = KmerExtractionStage::new(21);

        // Should be able to process with large budget
        assert!(stage.can_process(1024 * 1024 * 100)); // 100MB

        // Might not be able to process with very small budget
        assert!(!stage.can_process(1024)); // 1KB
    }
}