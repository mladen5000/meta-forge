// Assembly Performance Benchmarks
// Comprehensive criterion-based benchmarking suite for assembly pipeline phases

use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId, Throughput};
use meta_forge::assembly::laptop_assembly::{LaptopAssembler, LaptopConfig, CompactKmer, RollingKmerHash};
use meta_forge::core::data_structures::CorrectedRead;
use std::time::Duration;

/// Generate synthetic reads for benchmarking
fn generate_synthetic_reads(count: usize, length: usize, seed: u64) -> Vec<CorrectedRead> {
    use rand::{SeedableRng, Rng};
    use rand::rngs::StdRng;

    let mut rng = StdRng::seed_from_u64(seed);
    let bases = b"ACGT";

    (0..count)
        .map(|id| {
            let sequence: String = (0..length)
                .map(|_| bases[rng.gen_range(0..4)] as char)
                .collect();

            CorrectedRead {
                id,
                original: sequence.clone(),
                corrected: sequence,
                corrections: Vec::new(),
                quality_scores: vec![30; length],
                correction_metadata: Default::default(),
            }
        })
        .collect()
}

/// Benchmark k-mer counting with different read counts
fn bench_kmer_counting(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_counting");
    group.sample_size(20);
    group.measurement_time(Duration::from_secs(15));

    for read_count in [1000, 5000, 10000].iter() {
        group.throughput(Throughput::Elements(*read_count as u64));

        let reads = generate_synthetic_reads(*read_count, 150, 42);
        let config = LaptopConfig::auto_detect();
        let assembler = LaptopAssembler::new(config);

        group.bench_with_input(
            BenchmarkId::from_parameter(read_count),
            read_count,
            |b, _| {
                b.iter(|| {
                    // Benchmark full assembly (includes k-mer counting)
                    let _contigs = assembler.assemble(black_box(&reads))
                        .expect("Assembly failed");
                });
            },
        );
    }

    group.finish();
}

/// Benchmark graph construction with different read counts
fn bench_graph_construction(c: &mut Criterion) {
    let mut group = c.benchmark_group("graph_construction");
    group.sample_size(15);
    group.measurement_time(Duration::from_secs(20));

    for read_count in [1000, 5000, 10000].iter() {
        group.throughput(Throughput::Elements(*read_count as u64));

        let reads = generate_synthetic_reads(*read_count, 150, 42);
        let config = LaptopConfig::auto_detect();
        let assembler = LaptopAssembler::new(config);

        group.bench_with_input(
            BenchmarkId::from_parameter(read_count),
            read_count,
            |b, _| {
                b.iter(|| {
                    // Benchmark full assembly (graph construction is main cost)
                    let _contigs = assembler.assemble(black_box(&reads))
                        .expect("Assembly failed");
                });
            },
        );
    }

    group.finish();
}

/// Benchmark full assembly pipeline
fn bench_full_assembly(c: &mut Criterion) {
    let mut group = c.benchmark_group("full_assembly");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(30));

    for read_count in [1000, 5000, 10000].iter() {
        group.throughput(Throughput::Elements(*read_count as u64));

        let reads = generate_synthetic_reads(*read_count, 150, 42);
        let config = LaptopConfig::auto_detect();
        let assembler = LaptopAssembler::new(config);

        group.bench_with_input(
            BenchmarkId::from_parameter(read_count),
            read_count,
            |b, _| {
                b.iter(|| {
                    let _contigs = assembler.assemble(black_box(&reads))
                        .expect("Assembly failed");
                });
            },
        );
    }

    group.finish();
}

/// Benchmark k-mer operations (hash computation, rolling hash)
fn bench_kmer_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("kmer_operations");

    let sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let k = 21;

    // Benchmark hash computation
    group.bench_function("hash_computation", |b| {
        let kmer = CompactKmer::new(&sequence[0..k]).unwrap();
        b.iter(|| {
            black_box(kmer.rolling_hash());
        });
    });

    // Benchmark rolling hash update
    group.bench_function("rolling_hash_update", |b| {
        let mut rolling_hash = RollingKmerHash::new(k);
        rolling_hash.init(sequence[0..k].as_bytes());

        b.iter(|| {
            let _hash = rolling_hash.roll(black_box(b'A'), black_box(b'T'));
        });
    });

    // Benchmark k-mer creation from sequence
    group.bench_function("kmer_creation", |b| {
        b.iter(|| {
            for i in 0..=(sequence.len() - k) {
                let _kmer = CompactKmer::new(black_box(&sequence[i..i+k])).unwrap();
            }
        });
    });

    group.finish();
}

/// Benchmark memory efficiency with different k values
fn bench_memory_efficiency(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_efficiency");
    group.sample_size(10);

    for k in [15, 21, 27, 31].iter() {
        let reads = generate_synthetic_reads(10000, 150, 42);
        let mut config = LaptopConfig::auto_detect();
        config.max_k = *k;
        let assembler = LaptopAssembler::new(config);

        group.bench_with_input(
            BenchmarkId::new("k_value", k),
            k,
            |b, _| {
                b.iter(|| {
                    let _contigs = assembler.assemble(black_box(&reads))
                        .expect("Assembly failed");
                });
            },
        );
    }

    group.finish();
}

/// Benchmark parallel vs sequential processing
fn bench_parallelism(c: &mut Criterion) {
    let mut group = c.benchmark_group("parallelism");
    group.sample_size(10);

    let reads = generate_synthetic_reads(20000, 150, 42);

    for cpu_cores in [1, 2, 4, 8].iter() {
        let mut config = LaptopConfig::auto_detect();
        config.cpu_cores = *cpu_cores;
        let assembler = LaptopAssembler::new(config);

        group.bench_with_input(
            BenchmarkId::new("cpu_cores", cpu_cores),
            cpu_cores,
            |b, _| {
                b.iter(|| {
                    let _contigs = assembler.assemble(black_box(&reads))
                        .expect("Assembly failed");
                });
            },
        );
    }

    group.finish();
}

/// Benchmark contig generation from assembly
fn bench_contig_generation(c: &mut Criterion) {
    let mut group = c.benchmark_group("contig_generation");
    group.sample_size(10);

    let reads = generate_synthetic_reads(10000, 150, 42);
    let config = LaptopConfig::auto_detect();
    let assembler = LaptopAssembler::new(config);

    group.bench_function("generate_contigs", |b| {
        b.iter(|| {
            let _contigs = assembler.assemble(black_box(&reads))
                .expect("Assembly failed");
        });
    });

    group.finish();
}

/// Benchmark read lengths
fn bench_read_lengths(c: &mut Criterion) {
    let mut group = c.benchmark_group("read_lengths");

    for read_length in [100, 150, 250, 300].iter() {
        let reads = generate_synthetic_reads(5000, *read_length, 42);
        let config = LaptopConfig::auto_detect();
        let assembler = LaptopAssembler::new(config);

        group.bench_with_input(
            BenchmarkId::new("read_length", read_length),
            read_length,
            |b, _| {
                b.iter(|| {
                    let _contigs = assembler.assemble(black_box(&reads))
                        .expect("Assembly failed");
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    bench_kmer_counting,
    bench_graph_construction,
    bench_full_assembly,
    bench_kmer_operations,
    bench_memory_efficiency,
    bench_parallelism,
    bench_contig_generation,
    bench_read_lengths
);

criterion_main!(benches);
