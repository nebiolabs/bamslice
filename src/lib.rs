use anyhow::{Context, Result};
use log::{debug, info};
use rust_htslib::bam::{self, Read, Reader};
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{Read as IoRead, Seek, SeekFrom, Write};

// BGZF block magic bytes: gzip magic (0x1f 0x8b) + deflate method (0x08)
const BGZF_MAGIC: &[u8; 3] = &[0x1f, 0x8b, 0x08];

// FASTQ quality score encoding (Phred+33)
const PHRED_OFFSET: u8 = 33;

/// Find the next BGZF block at or after the given offset
/// # Errors
/// Returns an error if file operations fail or no BGZF block is found within a reasonable range
pub fn find_next_bgzf_block(file_path: &str, approximate_offset: u64) -> Result<u64> {
    const MAX_SCAN: usize = 65536; // Don't scan more than one max BGZF block size

    let mut file =
        File::open(file_path).with_context(|| format!("Failed to open BAM file: {file_path}"))?;

    // Special case: offset 0 is always the start
    if approximate_offset == 0 {
        return Ok(0);
    }

    file.seek(SeekFrom::Start(approximate_offset))
        .with_context(|| format!("Failed to seek to offset {approximate_offset}"))?;

    // Read chunks and scan for BGZF magic bytes
    let mut buffer = vec![0u8; 8192];
    let mut bytes_scanned = 0u64;

    loop {
        let bytes_read = file
            .read(&mut buffer)
            .context("Failed to read from BAM file while scanning for BGZF block")?;
        if bytes_read == 0 {
            anyhow::bail!("Reached EOF without finding BGZF block");
        }

        // Look for BGZF magic bytes (need at least 3 bytes)
        if bytes_read >= 3 {
            for i in 0..=(bytes_read - 3) {
                if &buffer[i..i + 3] == BGZF_MAGIC.as_slice() {
                    let block_offset = approximate_offset + bytes_scanned + i as u64;
                    debug!("Found BGZF block at offset {block_offset}");
                    return Ok(block_offset);
                }
            }
        }

        bytes_scanned += bytes_read as u64;
        if bytes_scanned > MAX_SCAN as u64 {
            anyhow::bail!(
                "Could not find BGZF block within {MAX_SCAN} bytes of offset {approximate_offset}"
            );
        }
    }
}

/// Convert a read to FASTQ format (interleaved - same format as samtools fastq)
#[must_use]
pub fn read_to_fastq(record: &bam::Record) -> String {
    let name = std::str::from_utf8(record.qname()).unwrap_or("non_utf8_read_name");
    let seq = record.seq();
    let qual = record.qual();

    // Pre-allocate with estimated capacity: @name/1\nseq\n+\nqual\n (extra 2 for /1 or /2)
    let mut fastq = String::with_capacity(name.len() + seq.len() * 2 + qual.len() + 6);

    fastq.push('@');
    fastq.push_str(name);

    // Determine read number and add /1 or /2 suffix for paired reads
    let read_num = if record.is_paired() && record.is_last_in_template() {
        fastq.push_str("/2");
        2
    } else {
        fastq.push_str("/1");
        1
    };

    // Add barcode information if present (BC tag in Illumina format: read:filtered:control:barcode)
    if let Ok(bam::record::Aux::String(barcode_str)) = record.aux(b"BC") {
        let filter_flag = if record.is_quality_check_failed() {
            'Y'
        } else {
            'N'
        };
        write!(fastq, " {read_num}:{filter_flag}:0:{barcode_str}").unwrap();
    }

    fastq.push('\n');

    // Sequence - use as_bytes() to get properly decoded sequence
    fastq.extend(seq.as_bytes().iter().map(|&b| b as char));
    fastq.push('\n');
    fastq.push_str("+\n");

    // Quality scores
    fastq.extend(qual.iter().map(|&q| (q + PHRED_OFFSET) as char));
    fastq.push('\n');

    fastq
}

/// Process blocks in the given byte range and output interleaved FASTQ
/// # Errors
/// Returns an error if file operations or BAM reading fails
pub fn process_blocks(
    input_path: &str,
    start_offset: u64,
    end_offset: u64,
    output: &mut dyn Write,
) -> Result<usize> {
    // Find the actual block start from the approximate offset
    let aligned_start = find_next_bgzf_block(input_path, start_offset)?;

    info!(
        "Processing byte range: {start_offset} to {end_offset} (aligned to block at {aligned_start})"
    );

    // Open BAM file
    let mut reader = Reader::from_path(input_path).context("Failed to open BAM file")?;

    // Seek to the aligned starting offset using virtual offset
    // Only seek if not starting from the beginning (offset 0 starts after header automatically)
    if aligned_start > 0 {
        let virtual_offset = aligned_start << 16; // https://samtools.github.io/hts-specs/SAMv1.pdf
        let virtual_offset_i64: i64 = virtual_offset
            .try_into()
            .context("Virtual offset too large to fit in i64")?;
        reader
            .seek(virtual_offset_i64)
            .context("Failed to seek to start offset")?;
    }

    let mut read_count = 0;
    let mut record = bam::Record::new();
    let mut blocks_processed = 0;
    let mut last_block_offset = aligned_start;

    loop {
        // Check current block position before reading
        let tell_result = reader.tell();
        let current_block_offset = u64::try_from(tell_result)
            .context("Invalid negative virtual offset from BAM reader")?
            >> 16;

        // Track when we enter a new block
        if current_block_offset != last_block_offset {
            blocks_processed += 1;
            last_block_offset = current_block_offset;
            debug!("Processing block at offset {current_block_offset}");
        }

        // Stop if we've reached a block at or past the end offset
        if current_block_offset >= end_offset {
            debug!("Reached end offset at block {current_block_offset}");
            break;
        }

        // Read next record
        match reader.read(&mut record) {
            Some(Ok(())) => {}
            None => break, // EOF
            Some(Err(e)) => {
                // If we're at or past the end offset, errors are expected
                // (could be reading beyond the valid region or hitting EOF marker)
                if current_block_offset >= end_offset {
                    debug!("Reached end of processing region at block {current_block_offset}");
                    break;
                }
                // Otherwise this is an unexpected error in the middle of processing
                return Err(e).context("Failed to read BAM record");
            }
        }

        // Write FASTQ record (interleaved format)
        let fastq = read_to_fastq(&record);
        output
            .write_all(fastq.as_bytes())
            .context("Failed to write FASTQ output")?;

        read_count += 1;

        if read_count % 100_000 == 0 {
            info!("Processed {read_count} reads...");
        }
    }

    info!("Completed: {blocks_processed} blocks, {read_count} reads");
    Ok(read_count)
}
