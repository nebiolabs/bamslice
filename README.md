# bamslice

[![Crates.io](https://img.shields.io/crates/v/bamslice.svg)](https://crates.io/crates/bamslice)
[![Documentation](https://docs.rs/bamslice/badge.svg)](https://docs.rs/bamslice)
[![CI](https://github.com/nebiolabs/bamslice/workflows/CI/badge.svg)](https://github.com/nebiolabs/bamslice/actions/workflows/ci.yml)
[![Coverage](https://github.com/nebiolabs/bamslice/workflows/Coverage/badge.svg)](https://github.com/nebiolabs/bamslice/actions/workflows/coverage.yml)
[![Security Audit](https://github.com/nebiolabs/bamslice/workflows/Security%20audit/badge.svg)](https://github.com/nebiolabs/bamslice/actions/workflows/audit.yml)

Extract specific byte ranges from BAM/CRAM files and convert to interleaved FASTQ format. Designed for parallel processing across compute nodes without requiring pre-indexing.

## Features

- **No pre-indexing required** - accepts approximate byte offsets
- **Auto-aligns to block boundaries** - finds the next valid BGZF block at or after the start offset
- **Byte-range based** - process arbitrary byte ranges for easy parallelization
- **No overlap** - using contiguous byte ranges guarantees no duplicate reads
- **Interleaved FASTQ output** - same format as `samtools fastq`
- **Parallel-ready** - designed for distributed processing

## Installation

```bash
cargo build --release
```

Binary: `target/release/bamslice`

## Usage

```bash
bamslice \
  --input input.bam \
  --start-offset 0 \
  --end-offset 10000000 \
  --output output.fastq
```

### Arguments

- `--input, -i`: Input BAM
- `--start-offset, -s`: Starting byte offset (will find next BGZF block at or after this offset)
- `--end-offset, -e`: Ending byte offset (will stop when reaching a block at or after this offset)
- `--output, -o`: Output FASTQ file (default: stdout)

### Examples

**Extract first half of file:**
```bash
FILE_SIZE=$(stat -f%z input.bam)  # macOS
# FILE_SIZE=$(stat -c%s input.bam)  # Linux
HALF=$((FILE_SIZE / 2))

bamslice -i input.bam -s 0 -e $HALF -o first_half.fastq
```

**Extract second half (no overlap!):**
```bash
bamslice -i input.bam -s $HALF -e $FILE_SIZE -o second_half.fastq
```

**Output to stdout:**
```bash
bamslice -i input.bam -s 0 -e 1000000 | head -n 4
```

## Parallel Processing

The tool uses byte ranges, making it trivial to parallelize without coordination

### Nextflow Example

See `example.nf` for a  pipeline that pipes bamslice output through fastp for QC/filtering.

```bash
nextflow run example.nf --bam input.bam --chunk_size 104857600
```

## How It Works

1. **BGZF Structure**: BAM files use BGZF (Blocked GZIP) - a series of independent compressed blocks
2. **Block Discovery**: Given a start offset, scans forward to find the next valid BGZF block (magic: `0x1f 0x8b 0x08`)
3. **Range Processing**: Processes all reads from blocks starting before `end_offset`
4. **No Overlap**: Each block is processed by exactly one job when using contiguous byte ranges
5. **FASTQ Output**: Converts BAM records to interleaved FASTQ format

### Why Byte Ranges?

- **No indexing overhead**: Don't need to scan the entire file first
- **Trivial parallelization**: Just choose your start/end offsets (see example nextflow)
- **No coordination**: Each process works independently
- **Guaranteed coverage**: Contiguous ranges ensure no reads are skipped
- **No duplication**: Block alignment ensures no reads are processed twice

## Testing

Run the test suite to verify correctness:

```bash
cargo test
```

lint the codebase:

```bash
make lint
```

## Development Commands

Run a coverage analysis:

```bash
make coverage && open target/coverage/html/index.html
```

Build a flamegraph for performance profiling:

```bash
make flamegraph && open flamegraph.svg
```

Run the performance benchmark:
```bash
make bench && open target/criterion/report/index.html
```

Release a new version:
```bash
echo "update Cargo.toml with new version"
git commit -m 'tag release vX.Y.Z'
git tag vX.Y.Z
git push --tags
cargo release
```

## License

AGPLv3 - See LICENSE file for details
