# bamslice

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

- `--input, -i`: Input BAM or CRAM file
- `--start-offset, -s`: Starting byte offset (will find next BGZF block at or after this offset)
- `--end-offset, -e`: Ending byte offset (will stop when reaching a block at or after this offset)
- `--reference, -r`: Reference FASTA for CRAM files (optional)
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

**Process CRAM with reference:**
```bash
bamslice -i input.cram -s 0 -e 1000000 -r reference.fa -o output.fastq
```

**Output to stdout:**
```bash
bamslice -i input.bam -s 0 -e 1000000 | head -n 4
```

## Parallel Processing

The tool uses byte ranges, making it trivial to parallelize without coordination:

```bash
#!/bin/bash
# parallel_bam2fastq.sh

BAM_FILE="input.bam"
CHUNK_SIZE=$((100 * 1024 * 1024))  # 100 MB chunks

# Get file size
FILE_SIZE=$(stat -f%z "$BAM_FILE")  # macOS
# FILE_SIZE=$(stat -c%s "$BAM_FILE")  # Linux

# Submit jobs (example: SLURM)
JOB_ID=0
for START in $(seq 0 $CHUNK_SIZE $FILE_SIZE); do
    END=$((START + CHUNK_SIZE))

    sbatch --job-name=bam2fq_${JOB_ID} \
           --output=logs/job_${JOB_ID}.log \
           --wrap="bamslice \
               -i $BAM_FILE \
               -s $START \
               -e $END \
               -o output/chunk_${JOB_ID}.fastq"
    
    JOB_ID=$((JOB_ID + 1))
done

# After all jobs complete, merge results
cat output/chunk_*.fastq > final.fastq
```

### SGE Example

```bash
JOB_ID=0
for START in $(seq 0 $CHUNK_SIZE $FILE_SIZE); do
    END=$((START + CHUNK_SIZE))

    qsub -N bam2fq_${JOB_ID} -o logs/job_${JOB_ID}.log -b y \
        bamslice -i $BAM_FILE -s $START -e $END \
        -o output/chunk_${JOB_ID}.fastq
    
    JOB_ID=$((JOB_ID + 1))
done
```

### Nextflow Example

See `example.nf` for a complete pipeline that pipes bamslice output through fastp for QC/filtering. The example uses `wc -l` to count reads, but includes a comment showing how to replace this with an aligner (e.g., `bwa mem`) or other downstream tools.

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
- **Trivial parallelization**: Just divide file size by number of jobs
- **No coordination**: Each process works independently
- **Guaranteed coverage**: Contiguous ranges ensure no reads are skipped
- **No duplication**: Block alignment ensures no reads are processed twice

## Testing

Run the test suite to verify correctness:

```bash
./test.sh
```

This compares output against `samtools fastq` using a split-and-merge approach.

## License

AGPLv3 - See LICENSE file for details
