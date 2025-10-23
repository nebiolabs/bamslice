use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read, Reader};
use std::fs::File;
use std::io::{Read as IoRead, Seek, SeekFrom, Write};

// BGZF block magic bytes: gzip magic (0x1f 0x8b) + deflate method (0x08)
const BGZF_MAGIC_BYTE1: u8 = 0x1f;
const BGZF_MAGIC_BYTE2: u8 = 0x8b;
const BGZF_MAGIC_BYTE3: u8 = 0x08;
const BGZF_HEADER_MIN_SIZE: usize = 18;

// FASTQ quality score encoding (Phred+33)
const PHRED_OFFSET: u8 = 33;

/// Find the next BGZF block at or after the given offset
pub fn find_next_bgzf_block(file_path: &str, approximate_offset: u64) -> Result<u64> {
    let mut file = File::open(file_path)?;

    // Special case: offset 0 is always the start
    if approximate_offset == 0 {
        return Ok(0);
    }

    file.seek(SeekFrom::Start(approximate_offset))?;

    // Read chunks and scan for BGZF magic, with overlap to handle headers spanning chunks
    let mut buffer = vec![0u8; 8192];
    let mut overlap = vec![0u8; BGZF_HEADER_MIN_SIZE];
    let mut overlap_len = 0;
    let mut bytes_scanned = 0u64;
    const MAX_SCAN: usize = 65536; // Don't scan more than one max block size

    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            anyhow::bail!("Reached EOF without finding BGZF block");
        }

        // Prepare search buffer combining overlap from previous read and new data
        let search_buf = if overlap_len > 0 {
            let mut combined = Vec::with_capacity(overlap_len + n);
            combined.extend_from_slice(&overlap[..overlap_len]);
            combined.extend_from_slice(&buffer[..n]);
            combined
        } else {
            buffer[..n].to_vec()
        };

        // Calculate the file offset where our search buffer starts
        // If we have overlap, we're searching from (bytes_scanned - overlap_len) into the file
        let search_buffer_file_offset = bytes_scanned.saturating_sub(overlap_len as u64);

        // Look for BGZF magic in search buffer
        let search_end = search_buf.len().saturating_sub(2); // Need at least 3 bytes for magic
        for i in 0..=search_end {
            if search_buf[i] == BGZF_MAGIC_BYTE1
                && i + 1 < search_buf.len()
                && search_buf[i + 1] == BGZF_MAGIC_BYTE2
                && i + 2 < search_buf.len()
                && search_buf[i + 2] == BGZF_MAGIC_BYTE3
            {
                // Found magic at position i in search_buf
                // Absolute file offset = approximate_offset + where search_buf starts + position in buffer
                let block_offset = approximate_offset + search_buffer_file_offset + i as u64;
                eprintln!("Found BGZF block at offset {}", block_offset);
                return Ok(block_offset);
            }
        }

        // Save last BGZF_HEADER_MIN_SIZE bytes as overlap for next iteration
        if n >= BGZF_HEADER_MIN_SIZE {
            overlap[..BGZF_HEADER_MIN_SIZE].copy_from_slice(&buffer[n - BGZF_HEADER_MIN_SIZE..n]);
            overlap_len = BGZF_HEADER_MIN_SIZE;
        } else {
            overlap[..n].copy_from_slice(&buffer[..n]);
            overlap_len = n;
        }

        bytes_scanned += n as u64;
        if bytes_scanned > MAX_SCAN as u64 {
            anyhow::bail!(
                "Could not find BGZF block within {} bytes of offset {}",
                MAX_SCAN,
                approximate_offset
            );
        }
    }
}

/// Convert a read to FASTQ format (interleaved - same format as samtools fastq)
pub fn read_to_fastq(record: &bam::Record) -> String {
    let name = std::str::from_utf8(record.qname()).unwrap_or("unknown");
    let seq = record.seq();
    let qual = record.qual();

    // Pre-allocate with estimated capacity: @name\nseq\n+\nqual\n
    let mut fastq = String::with_capacity(name.len() + seq.len() * 2 + qual.len() + 4);

    fastq.push('@');
    fastq.push_str(name);
    fastq.push('\n');

    // Sequence
    for i in 0..seq.len() {
        fastq.push(match seq[i] {
            1 => 'A',
            2 => 'C',
            4 => 'G',
            8 => 'T',
            15 => 'N',
            _ => 'N',
        });
    }
    fastq.push('\n');
    fastq.push_str("+\n");

    // Quality scores
    for &q in qual {
        fastq.push((q + PHRED_OFFSET) as char);
    }
    fastq.push('\n');

    fastq
}

/// Process blocks in the given byte range and output interleaved FASTQ
pub fn process_blocks(
    input_path: &str,
    start_offset: u64,
    end_offset: u64,
    reference: Option<&str>,
    output: &mut dyn Write,
) -> Result<usize> {
    // Find the actual block start from the approximate offset
    let aligned_start = find_next_bgzf_block(input_path, start_offset)?;

    eprintln!(
        "Processing byte range: {} to {} (aligned to block at {})",
        start_offset, end_offset, aligned_start
    );

    // Open BAM/CRAM file
    let mut reader = Reader::from_path(input_path).context("Failed to open BAM/CRAM file")?;

    // Set reference if provided (for CRAM files)
    if let Some(ref_path) = reference {
        reader
            .set_reference(ref_path)
            .context("Failed to set reference for CRAM")?;
    }

    // Seek to the aligned starting offset using virtual offset
    // Only seek if not starting from the beginning (offset 0 starts after header automatically)
    if aligned_start > 0 {
        let virtual_offset = aligned_start << 16; // https://samtools.github.io/hts-specs/SAMv1.pdf
        reader
            .seek(virtual_offset as i64)
            .context("Failed to seek to start offset")?;
    }

    let mut read_count = 0;
    let mut record = bam::Record::new();
    let mut blocks_processed = 0;
    let mut last_block_offset = aligned_start;

    loop {
        // Check current block position before reading
        let current_block_offset = (reader.tell() as u64) >> 16;

        // Track when we enter a new block
        if current_block_offset != last_block_offset {
            blocks_processed += 1;
            last_block_offset = current_block_offset;
            eprintln!("  Processing block at offset {}...", current_block_offset);
        }

        // Stop if we've reached a block at or past the end offset
        if current_block_offset >= end_offset {
            eprintln!("Reached end offset at block {}", current_block_offset);
            break;
        }

        // Read next record
        match reader.read(&mut record) {
            Some(Ok(())) => {}
            None => break, // EOF
            Some(Err(e)) => {
                // Check if we've hit the end of valid data (EOF marker or beyond)
                // This can happen when end_offset extends to the file size
                let err_msg = format!("{}", e);
                if err_msg.contains("invalid record") || err_msg.contains("EOF") {
                    eprintln!("Reached end of valid BAM records");
                    break;
                }
                return Err(e.into());
            }
        }

        // Write FASTQ record (interleaved format)
        let fastq = read_to_fastq(&record);
        output
            .write_all(fastq.as_bytes())
            .context("Failed to write FASTQ output")?;

        read_count += 1;

        if read_count % 100000 == 0 {
            eprintln!("  Processed {} reads...", read_count);
        }
    }

    eprintln!(
        "Processed {} blocks, {} reads",
        blocks_processed, read_count
    );

    eprintln!("Completed: {} reads processed", read_count);
    Ok(read_count)
}
