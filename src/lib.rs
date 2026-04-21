//! Extract byte ranges from BAM files and output as interleaved FASTQ or BAM.
//!
//! This library provides efficient extraction of specific byte ranges from BAM (Binary Alignment/Map)
//! files and outputs them as interleaved FASTQ or BAM. This enables parallel processing of large
//! BAM files by splitting them into chunks that can be processed independently.
//!
//! # Key Features
//!
//! - **Block-aligned extraction**: Automatically aligns to BGZF block boundaries for valid data
//! - **Paired-read aware**: Ensures read pairs are kept together across chunk boundaries
//! - **Interleaved FASTQ output**: Compatible with `samtools fastq` interleaved format
//! - **BAM output**: Raw-copies middle BGZF blocks, only recompressing at slice boundaries
//! - **Barcode support**: Preserves BC tags in FASTQ headers
//!
//! # BGZF Block Format
//!
//! BAM files use BGZF (Blocked GNU Zip Format) compression. Each block is independently
//! compressed, allowing random access. This library scans for valid BGZF block headers
//! using an 8-byte signature and validates blocks by decompressing and checking for valid BAM records.
//!
//! # Implementation Notes
//!
//! - Read pairs are kept together by reading one extra record past the end boundary if needed

use anyhow::{Context, Result};
use log::{debug, info};
use noodles::bam;
use noodles::bgzf;
use noodles::sam;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};

/// Output format for extracted reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    /// Interleaved FASTQ format (compatible with `samtools fastq`)
    Fastq,
    /// BAM format
    Bam,
}

// BGZF block magic bytes: gzip magic (0x1f 0x8b) + deflate method (0x08) + extra field flag (0x04)
const BGZF_MAGIC: &[u8; 4] = &[0x1f, 0x8b, 0x08, 0x04];

// Maximum compressed BGZF block size per spec (64KB)
const BGZF_MAX_BLOCK_SIZE: usize = 65536;

// FASTQ quality score encoding (Phred+33)
const PHRED_OFFSET: u8 = 33;

/// Check if buffer at position i contains a valid BGZF block header.
///
/// BGZF header structure per [SAM spec section 4.1](https://samtools.github.io/hts-specs/SAMv1.pdf):
/// - Bytes 0-3: gzip magic (1f 8b) + compression method (08) + flags (04)
/// - Bytes 12-13: subfield identifier SI1='B' (0x42), SI2='C' (0x43)
/// - Bytes 14-15: subfield length SLEN=2 (0x02, 0x00)
///
/// This checks 8 specific bytes (64 bits), giving a false positive rate of ~5e-11 per GB.
/// A 3-byte gzip magic check (24 bits) would have ~60 false positives per GB.
fn is_valid_bgzf_header(buffer: &[u8], i: usize) -> bool {
    if i + 16 > buffer.len() {
        return false;
    }

    // Check gzip magic + extra field flag
    if &buffer[i..i + 4] != BGZF_MAGIC.as_slice() {
        return false;
    }

    // Check BGZF subfield identifier 'BC' at bytes 12-13
    if buffer[i + 12] != 0x42 || buffer[i + 13] != 0x43 {
        return false;
    }

    // Check SLEN = 2 at bytes 14-15
    if buffer[i + 14] != 0x02 || buffer[i + 15] != 0x00 {
        return false;
    }

    true
}

/// Find the next BGZF block at or after the given offset.
/// Returns `None` if EOF is reached. Offset 0 is always valid.
fn find_next_bgzf_block<R>(
    reader: &mut R,
    approximate_offset: u64,
    buffer: &mut Vec<u8>,
) -> Result<Option<u64>>
where
    R: Read + Seek,
{
    if approximate_offset == 0 {
        return Ok(Some(0));
    }

    reader
        .seek(SeekFrom::Start(approximate_offset))
        .with_context(|| format!("Failed to seek to offset {approximate_offset}"))?;

    // Ensure buffer is large enough
    let required_size = BGZF_MAX_BLOCK_SIZE + 15;
    // Resize to exact required size (truncates if larger, extends with 0 if smaller)
    buffer.resize(required_size, 0);

    let mut bytes_read = 0;

    while bytes_read < buffer.len() {
        let n = reader.read(&mut buffer[bytes_read..])?;
        if n == 0 {
            if bytes_read == 0 {
                // EOF
                info!("Reached EOF at offset {approximate_offset} without finding BGZF block");
                return Ok(None);
            }
            break; // EOF
        }
        bytes_read += n;
    }

    // Scan for BGZF header
    for i in 0..bytes_read {
        if is_valid_bgzf_header(buffer, i) {
            let block_offset = approximate_offset + i as u64;
            debug!("Found BGZF block signature at offset {block_offset}");
            return Ok(Some(block_offset));
        }
    }

    // If we hit EOF (read less than requested), then no block exists.
    if bytes_read < buffer.len() {
        return Ok(None);
    }

    anyhow::bail!(
        "Scanned {bytes_read} bytes from offset {approximate_offset} without finding a BGZF block."
    );
}

/// Write a BAM record as interleaved FASTQ (compatible with `samtools fastq`).
/// Header: `@NAME/[12] [12]:FILTER:0:BARCODE` where FILTER is Y/N for QC fail.
fn write_fastq_to(record: &bam::Record, output: &mut dyn Write) -> std::io::Result<()> {
    // Read name
    let name = record
        .name()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::InvalidData, "Missing read name"))?;

    output.write_all(b"@")?;
    output.write_all(name.as_ref())?;

    // Determine read number and add /1 or /2 suffix only for paired reads
    // bam::Record::flags() returns Flags directly (inherent method)
    let flags = record.flags();
    let is_paired = flags.is_segmented();

    if is_paired {
        if flags.is_last_segment() {
            output.write_all(b"/2")?;
        } else {
            output.write_all(b"/1")?;
        }
    }

    // Add barcode information if present (BC tag)
    let data = record.data();
    for result in data.iter() {
        let (tag, value) =
            result.map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        if tag.as_ref() == b"BC" {
            // Value is an enum.
            use noodles::sam::alignment::record::data::field::Value;
            if let Value::String(s) = value {
                let filter_flag = if flags.is_qc_fail() { 'Y' } else { 'N' };
                // s is likely &BStr or similar
                let bc_str = std::str::from_utf8(s.as_ref()).unwrap_or("INVALID");
                if is_paired {
                    let read_num = if flags.is_last_segment() { 2 } else { 1 };
                    write!(output, " {read_num}:{filter_flag}:0:{bc_str}")?;
                } else {
                    write!(output, " 0:{filter_flag}:0:{bc_str}")?;
                }
            }
            break;
        }
    }

    output.write_all(b"\n")?;

    // Sequence
    let sequence = record.sequence();
    for base in sequence.iter() {
        output.write_all(&[base])?;
    }

    output.write_all(b"\n+\n")?;

    // Quality scores
    let quality_scores = record.quality_scores();
    for score in quality_scores.iter() {
        output.write_all(&[score + PHRED_OFFSET])?;
    }
    output.write_all(b"\n")?;

    Ok(())
}

/// Check if buffer contains a valid BAM record header per SAM spec Section 4.2.
fn is_valid_record_start(buf: &[u8], len: usize) -> bool {
    if len < 36 {
        return false;
    }

    // Sanity check block size: must be at least 32 (fixed fields) and not absurdly large
    let block_size = i32::from_le_bytes(buf[0..4].try_into().unwrap());
    if !(32..=100_000_000).contains(&block_size) {
        return false;
    }

    let ref_id = i32::from_le_bytes(buf[4..8].try_into().unwrap());
    if ref_id < -1 {
        return false;
    }

    let pos = i32::from_le_bytes(buf[8..12].try_into().unwrap());
    if pos < -1 {
        return false;
    }

    let l_read_name = buf[12];
    if l_read_name < 2 {
        return false;
    }

    let n_cigar_op = u16::from_le_bytes(buf[16..18].try_into().unwrap());
    let l_seq = i32::from_le_bytes(buf[20..24].try_into().unwrap());

    if l_seq < 0 {
        return false;
    }

    // Spec: read_name is NUL-terminated char[l_read_name]
    // The last byte of the read name must be 0.
    // Fixed fields size is 36 bytes (up to start of read_name).
    // read_name starts at offset 36.
    if len >= 36 + l_read_name as usize && buf[36 + l_read_name as usize - 1] != 0 {
        return false;
    }

    // Consistency check: are the numbers we read from the header patterns self-consistent?
    // block_size = 32 + l_read_name + n_cigar_op * 4 + (l_seq + 1) / 2 + l_seq + aux_len
    // Use i64 to prevent overflow
    let min_size = 32
        + i64::from(l_read_name)
        + i64::from(n_cigar_op) * 4
        + (i64::from(l_seq) + 1) / 2
        + i64::from(l_seq);
    if i64::from(block_size) < min_size {
        return false;
    }

    true
}

/// Decompress a BGZF block and scan for the first valid BAM record header.
/// Returns the byte offset within the decompressed block, or `None`.
fn find_first_record_in_block<R>(
    reader: &mut R,
    block_start: u64,
    buffer: &mut Vec<u8>,
) -> Result<Option<u64>>
where
    R: Read + Seek,
{
    reader.seek(SeekFrom::Start(block_start))?;
    let mut reader = flate2::read::GzDecoder::new(reader);

    // Read enough data to find a record. 64KB is typical max block size.
    buffer.resize(BGZF_MAX_BLOCK_SIZE, 0);

    let bytes_read = std::io::Read::read(&mut reader, buffer)?;

    if bytes_read < 36 {
        return Ok(None);
    }

    // Try every byte offset
    // We scan byte-by-byte because a BAM record can start at any offset within the decompressed block.
    // There is no alignment guarantee for records within the block.
    for i in 0..bytes_read - 36 {
        if is_valid_record_start(&buffer[i..], bytes_read - i) {
            return Ok(Some(i as u64));
        }
    }
    Ok(None)
}

/// Validate a BGZF block by decompressing it and finding the first BAM record.
/// Returns the virtual position of that record, or `None` if the block is invalid.
fn validate_block<R>(
    reader: &mut R,
    aligned_start: u64,
    buffer: &mut Vec<u8>,
) -> Result<Option<u64>>
where
    R: Read + Seek,
{
    reader.seek(SeekFrom::Start(aligned_start))?;

    // Try to find a record start within the decompressed block
    match find_first_record_in_block(reader, aligned_start, buffer) {
        Ok(Some(offset)) => {
            let offset_u16 =
                u16::try_from(offset).map_err(|e| anyhow::anyhow!("Invalid offset: {e}"))?;
            let new_virtual_pos = bgzf::VirtualPosition::try_from((aligned_start, offset_u16))
                .map_err(|e| anyhow::anyhow!("Invalid virtual position: {e}"))?;
            Ok(Some(u64::from(new_virtual_pos)))
        }
        Ok(None) | Err(_) => Ok(None), // Treats errors as invalid blocks
    }
}

struct ProcessResult {
    parsed_reads: usize,
    /// block offset where raw copying can start
    clean_boundary_offset: Option<u64>,
}

/// Read records and write them using the provided write function.
///
/// For BAM output, returns early at the first clean BGZF block boundary
/// (record at decompressed offset 0, at a pair boundary) to enable raw block copy.
fn process_records<R, F>(
    reader: &mut bam::io::Reader<bgzf::io::Reader<R>>,
    start_aligned_offset: u64,
    end_offset: u64,
    write_record: &mut F,
    format: OutputFormat,
) -> Result<ProcessResult>
where
    R: Read + Seek,
    F: FnMut(&bam::Record) -> std::io::Result<()>,
{
    let mut parsed_reads = 0;
    let mut record = bam::Record::default();
    let mut blocks_processed = 0;
    let mut prev_block_offset = start_aligned_offset;
    let mut prev_was_read1: Option<bool> = None;

    loop {
        let virtual_pos = reader.get_ref().virtual_position();
        let current_block_offset = virtual_pos.compressed();

        if current_block_offset != prev_block_offset {
            blocks_processed += 1;

            // For BAM output: check for clean boundary before consuming data from the new block.
            // A clean boundary means the next record starts at byte 0 of this block
            // (no record spans across the block boundary) and we're at a pair boundary.
            if matches!(format, OutputFormat::Bam)
                && virtual_pos.uncompressed() == 0
                && parsed_reads > 0
                && prev_was_read1 != Some(true)
                && current_block_offset < end_offset
            {
                info!(
                    "Clean block boundary at offset {current_block_offset} \
                     after {blocks_processed} blocks, {parsed_reads} records"
                );
                return Ok(ProcessResult {
                    parsed_reads,
                    clean_boundary_offset: Some(current_block_offset),
                });
            }

            prev_block_offset = current_block_offset;
            // If this block starts at/after end_offset
            if current_block_offset >= end_offset {
                // Check if we just wrote a read1 that needs its mate
                // If so, we must read one more record (the mate) before stopping
                if prev_was_read1 == Some(true) {
                    debug!(
                        "Reached end_offset at block {current_block_offset}, but need to read mate for previous read1"
                    );
                    // Don't break yet - continue to read the mate
                } else {
                    debug!("Reached end of processing region at block {current_block_offset}");
                    break;
                }
            } else {
                debug!("Processing block at offset {current_block_offset}");
            }
        }

        match reader.read_record(&mut record) {
            Ok(0) => break, // EOF
            Ok(_) => {}
            Err(e) => {
                if current_block_offset >= end_offset {
                    debug!("Reached end of processing region at block {current_block_offset}");
                    break;
                }
                return Err(anyhow::Error::new(e)).context("Failed to read BAM record");
            }
        }

        let flags = record.flags();
        let current_is_read1 = flags.is_first_segment();

        // Make sure that both mates stay together
        if flags.is_segmented() {
            // Skip orphaned read2 at the start of processing region
            if parsed_reads == 0 && !current_is_read1 {
                debug!(
                    "Skipping read2 at start of processing region, it was handled in the previous region"
                );
                continue;
            }

            // Strict validation for collated reads
            if prev_was_read1 == Some(true) && current_is_read1 {
                anyhow::bail!(
                    "Read '{}' has READ_PAIRED set but follows another R1 with no \
                     intervening R2. Either the BAM is not collated, or it contains \
                     only R1 reads from a paired-end run (e.g. extracted with \
                     `samtools view -f 64`).",
                    String::from_utf8_lossy(record.name().unwrap_or_default().as_ref())
                );
            }
            if prev_was_read1 == Some(false) && !current_is_read1 {
                anyhow::bail!(
                    "Found R2 read ({}) without preceding R1. Input BAM must be collated.",
                    String::from_utf8_lossy(record.name().unwrap_or_default().as_ref())
                );
            }
            prev_was_read1 = Some(current_is_read1);
        } else {
            prev_was_read1 = Some(false);
        }
        write_record(&record).context("Failed to write record")?;

        parsed_reads += 1;

        if parsed_reads % 100_000 == 0 {
            info!("Processed {parsed_reads} reads...");
        }

        // After writing a record, check if we're beyond end_offset and:
        // - For paired reads: just completed a pair (wrote read2)
        // - For single-end reads: just wrote a read (no pairing needed)
        if current_block_offset >= end_offset {
            let can_stop = if flags.is_segmented() {
                !current_is_read1 // Stop after read2 (completed pair)
            } else {
                true // Stop after any single-end read
            };
            if can_stop {
                debug!("Completed beyond end_offset, stopping");
                break;
            }
        }
    }

    info!("Completed: {blocks_processed} blocks, {parsed_reads} reads");
    Ok(ProcessResult {
        parsed_reads,
        clean_boundary_offset: None,
    })
}

/// Read the BAM header from the beginning of a BAM file.
fn read_bam_header(input_path: &str) -> Result<sam::Header> {
    let file =
        File::open(input_path).with_context(|| format!("Failed to open {input_path} for header"))?;
    let mut reader = bam::io::Reader::new(std::io::BufReader::new(file));
    reader.read_header().context("Failed to read BAM header")
}

/// Raw-copy middle BGZF blocks and handle the end boundary.
///
/// Copies raw BGZF blocks from `raw_start` to the first block at/after `end_offset`,
/// then checks whether an R2 mate record from the end boundary block needs to be
/// appended to complete a pair. Finishes by writing the BGZF EOF block.
fn raw_copy_middle_and_end(
    output: &mut dyn Write,
    header: &sam::Header,
    raw_start: u64,
    end_offset: u64,
    input_path: &str,
    buffer: &mut Vec<u8>,
) -> Result<(usize, u64)> {
    // Find end of raw copy range: the first BGZF block at/after end_offset
    let mut end_scan = std::io::BufReader::new(
        File::open(input_path).with_context(|| format!("Failed to open {input_path}"))?,
    );
    let end_block = find_next_bgzf_block(&mut end_scan, end_offset, buffer)?;
    let raw_end = end_block.unwrap_or(end_offset);

    let mut extra_reads = 0;
    let mut bytes_copied = 0;

    if raw_start < raw_end {
        // Copy raw BGZF blocks directly to output — no decompression or recompression needed
        let mut copy_file =
            File::open(input_path).with_context(|| format!("Failed to open {input_path}"))?;
        copy_file.seek(SeekFrom::Start(raw_start))?;
        bytes_copied = raw_end - raw_start;
        std::io::copy(&mut (&mut copy_file).take(bytes_copied), &mut *output)?;

        // Count records in the raw-copied range (decompression only, no recompression)
        let mut count_file =
            File::open(input_path).with_context(|| format!("Failed to open {input_path}"))?;
        count_file.seek(SeekFrom::Start(raw_start))?;
        let mut count_reader =
            bam::io::Reader::from(bgzf::io::Reader::new(std::io::BufReader::new(count_file)));
        let mut count_record = bam::Record::default();
        loop {
            let vpos = count_reader.get_ref().virtual_position();
            if vpos.compressed() >= raw_end {
                break;
            }
            match count_reader.read_record(&mut count_record) {
                Ok(0) | Err(_) => break,
                Ok(_) => extra_reads += 1,
            }
        }

        info!("Raw-copied {bytes_copied} bytes ({raw_start}..{raw_end}), {extra_reads} records");
    }

    // Handle end boundary: only include the first record from the end block
    // if it's R2 (completing a pair whose R1 was in the last raw-copied block).
    // The next chunk will process this block from the start, skipping the orphaned R2.
    if let Some(end_block_start) = end_block {
        let mut end_file = std::io::BufReader::new(
            File::open(input_path).with_context(|| format!("Failed to open {input_path}"))?,
        );
        if let Some(record_offset) =
            find_first_record_in_block(&mut end_file, end_block_start, buffer)?
        {
            end_file.seek(SeekFrom::Start(end_block_start))?;
            let mut bgzf_end = bgzf::io::Reader::new(end_file);
            let offset_u16 = u16::try_from(record_offset)
                .map_err(|e| anyhow::anyhow!("Invalid offset: {e}"))?;
            let vpos = bgzf::VirtualPosition::try_from((end_block_start, offset_u16))
                .map_err(|e| anyhow::anyhow!("Invalid virtual position: {e}"))?;
            bgzf_end.seek(vpos)?;
            let mut bam_end = bam::io::Reader::from(bgzf_end);

            let mut end_record = bam::Record::default();
            if bam_end.read_record(&mut end_record)? > 0 {
                let flags = end_record.flags();
                if flags.is_segmented() && !flags.is_first_segment() {
                    let mut end_writer = bam::io::Writer::new(&mut *output);
                    end_writer
                        .write_record(header, &end_record)
                        .context("Failed to write end-boundary R2")?;
                    end_writer
                        .try_finish()
                        .context("Failed to finalize end-boundary BGZF")?;
                    extra_reads += 1;
                    info!("Included R2 mate from end boundary block at {end_block_start}");
                    return Ok((extra_reads, bytes_copied));
                }
            }
        }
    }

    // Write BGZF EOF block
    bgzf::io::Writer::new(&mut *output)
        .try_finish()
        .context("Failed to write BGZF EOF")?;

    Ok((extra_reads, bytes_copied))
}

/// Process a byte range from a BAM file and output in the specified format.
///
/// For BAM output, middle BGZF blocks are raw-copied from the input when a clean
/// block boundary is found, avoiding recompression for the bulk of the data.
///
/// # Errors
///
/// Returns an error if the input file cannot be opened or read, the BAM data is
/// malformed, or writing to the output stream fails.
pub fn process_blocks(
    input_path: &str,
    start_offset: u64,
    end_offset: u64,
    output: &mut dyn Write,
    format: OutputFormat,
) -> Result<usize> {
    let file = File::open(input_path).with_context(|| format!("Failed to open {input_path}"))?;
    let mut reader = std::io::BufReader::new(file);

    let mut current_offset = start_offset;
    let mut buffer = vec![0u8; BGZF_MAX_BLOCK_SIZE + 15];

    loop {
        let Some(block_start) = find_next_bgzf_block(&mut reader, current_offset, &mut buffer)?
        else {
            break;
        };

        if block_start >= end_offset {
            break;
        }

        if let Some(virtual_pos) = validate_block(&mut reader, block_start, &mut buffer)? {
            let aligned_start = virtual_pos >> 16;
            info!(
                "Processing byte range: {start_offset} to {end_offset} (aligned to block at {aligned_start})"
            );

            reader.seek(SeekFrom::Start(aligned_start))?;
            let mut bgzf_reader = bgzf::io::Reader::new(&mut reader);
            bgzf_reader.seek(bgzf::VirtualPosition::from(virtual_pos))?;
            let mut bam_reader = bam::io::Reader::from(bgzf_reader);

            return match format {
                OutputFormat::Fastq => {
                    let mut write_fn = |record: &bam::Record| write_fastq_to(record, output);
                    let result = process_records(
                        &mut bam_reader,
                        aligned_start,
                        end_offset,
                        &mut write_fn,
                        OutputFormat::Fastq,
                    )?;
                    Ok(result.parsed_reads)
                }
                OutputFormat::Bam => {
                    let header = read_bam_header(input_path)?;
                    let mut bam_writer = bam::io::Writer::new(output);
                    bam_writer.write_header(&header)?;

                    let result;
                    {
                        let mut write_fn = |record: &bam::Record| -> std::io::Result<()> {
                            bam_writer.write_record(&header, record)
                        };
                        result = process_records(
                            &mut bam_reader,
                            aligned_start,
                            end_offset,
                            &mut write_fn,
                            OutputFormat::Bam,
                        )?;
                    }

                    let mut total = result.parsed_reads;

                    if let Some(raw_start) = result.clean_boundary_offset {
                        // Decompose the BAM writer to access the raw output stream.
                        // flush() compresses any buffered data without writing an EOF block.
                        // into_inner() then yields the underlying writer.
                        let mut bgzf_writer = bam_writer.into_inner();
                        bgzf_writer.flush()?;
                        let raw_output = bgzf_writer.into_inner();

                        let (extra, _bytes) = raw_copy_middle_and_end(
                            raw_output,
                            &header,
                            raw_start,
                            end_offset,
                            input_path,
                            &mut buffer,
                        )?;
                        total += extra;
                    } else {
                        bam_writer
                            .try_finish()
                            .context("Failed to finalize BGZF output")?;
                    }

                    Ok(total)
                }
            };
        }
        debug!("Invalid block at {block_start}, continuing search");
        current_offset = block_start + 1;
    }
    Ok(0)
}
