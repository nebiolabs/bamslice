use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use std::io::Cursor;

/// Benchmark full processing pipeline with different byte ranges
fn bench_process_blocks(c: &mut Criterion) {
    let test_file = "tests/fixtures/dnbseq-test1.bam";

    if !std::path::Path::new(test_file).exists() {
        eprintln!("Skipping bench_process_blocks: {} not found", test_file);
        return;
    }

    let mut group = c.benchmark_group("process_blocks");

    // Different range sizes to benchmark
    let ranges = [
        ("small_500kb", 10_000u64, 510_000u64),
        ("medium_1mb", 10_000u64, 1_010_000u64),
        ("large_5mb", 10_000u64, 5_010_000u64),
    ];

    for (name, start, end) in ranges {
        let range_size = end - start;
        group.throughput(Throughput::Bytes(range_size));

        group.bench_with_input(
            BenchmarkId::from_parameter(name),
            &(start, end),
            |b, &(start, end)| {
                b.iter(|| {
                    let mut output = Cursor::new(Vec::new());
                    bamslice::process_blocks(
                        test_file,
                        black_box(start),
                        black_box(end),
                        &mut output,
                    )
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_process_blocks);
criterion_main!(benches);
