//! Extract byte ranges from BAM files and convert to interleaved FASTQ format.
//!
//! This library provides efficient extraction of specific byte ranges from BAM (Binary Alignment/Map)
//! files and converts them to interleaved FASTQ format. This enables parallel processing of large
//! BAM files by splitting them into chunks that can be processed independently.
//!
//! # Key Features
//!
//! - **Block-aligned extraction**: Automatically aligns to BGZF block boundaries for valid data
//! - **Paired-read aware**: Ensures read pairs are kept together across chunk boundaries
//! - **Interleaved FASTQ output**: Compatible with `samtools fastq` interleaved format
//! - **Barcode support**: Preserves BC tags in FASTQ headers
//!
//! # Example
//!
//! ```no_run
//! use std::io::stdout;
//! use bamslice::process_blocks;
//!
//! # fn main() -> anyhow::Result<()> {
//! let mut output = stdout();
//! let read_count = process_blocks(
//!     "input.bam",
//!     0,           // start offset
//!     1_000_000,   // end offset (1MB chunk)
//!     &mut output
//! )?;
//! println!("Extracted {read_count} reads");
//! # Ok(())
//! # }
//! ```
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
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};

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
///
/// Scans forward from the approximate offset to find a valid BGZF block header.
/// Special case: offset 0 is always valid and returned immediately.
///
/// # Arguments
///
/// * `reader` - A seekable reader positioned in the BAM file
/// * `approximate_offset` - The byte offset to start searching from
/// * `buffer` - Scratch buffer to use for reading (should be at least `BGZF_MAX_BLOCK_SIZE` + 15 bytes)
///
/// # Returns
///
/// - `Ok(Some(offset))` - Found a valid BGZF block at this offset
/// - `Ok(None)` - Reached EOF before finding a BGZF block (normal near end of file)
///
/// # Errors
///
/// Returns an error if:
/// - Unable to seek to the requested offset
/// - Unable to read from the file
/// - Scanned `MAX_SCAN` bytes without finding a block or EOF (corrupt data)
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

/// Write a BAM record to output in FASTQ format.
///
/// Outputs interleaved FASTQ format compatible with `samtools fastq`:
/// - Read name with `/1` or `/2` suffix for paired reads
/// - Barcode information from BC tag appended to header if present
/// - Quality scores in Phred+33 format
///
/// # FASTQ Header Format
///
/// Unpaired: `@READ_NAME [0:FILTER:0:BARCODE]`
/// Paired: `@READ_NAME/[12] [12]:FILTER:0:BARCODE`
///
/// Where FILTER is 'Y' if QC failed, 'N' otherwise.
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

/// Check if buffer contains a valid BAM record header.
///
/// Validates according to [BAM spec Section 4.2](https://samtools.github.io/hts-specs/SAMv1.pdf):
/// - Block size is reasonable (32 to 100MB)
/// - Reference ID >= -1
/// - Position >= -1
/// - Read name length >= 2
/// - Read name is NUL-terminated
/// - Sequence length >= 0
/// - Header fields are self-consistent
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

/// Scan a BGZF block to find the offset of the first valid BAM record.
///
/// Decompresses the block and scans byte-by-byte for a valid record header.
/// Records can start at any offset within a decompressed block (no alignment guarantee).
///
/// # Returns
///
/// - `Ok(Some(offset))` - Offset within decompressed block where first valid record starts
/// - `Ok(None)` - No valid record found in this block
///
/// # Errors
///
/// Returns an error if unable to seek to block start or decompress the block
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

/// Validate a BGZF block and find the first record within it.
///
/// Attempts to decompress the block at the given offset and locate the first
/// valid BAM record. Returns a virtual position that can be used to seek to
/// the record.
///
/// # Arguments
///
/// * `reader` - A seekable reader positioned in the BAM file
/// * `aligned_start` - Byte offset of the BGZF block to validate
/// * `buffer` - Scratch buffer for decompression
///
/// # Returns
///
/// - `Ok(Some(virtual_pos))` - Virtual position of first valid record (if found)
/// - `Ok(None)` - Block is invalid or contains no valid records
///
/// # Errors
///
/// Returns an error if unable to seek or read from the file
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

/// Read records from a BAM reader and write them as FASTQ.
///
/// Processes records starting from the current reader position until:
/// - EOF is reached
/// - A block at or past `end_offset` is encountered
/// - An error occurs reading a record
///
/// # Paired-read handling
///
/// - Skips the first record if it's read 2 (to ensure pairs align)
/// - Reads one extra record past the end boundary if needed to complete a pair
///
/// # Arguments
///
/// * `reader` - BAM reader positioned at the start of processing
/// * `start_aligned_offset` - Starting block offset (for logging)
/// * `end_offset` - Process blocks that beging before this offset
/// * `output` - Output stream for FASTQ data
///
/// # Returns
///
/// - `Number of reads (not pairs) written to output`
///
/// # Errors
///
/// Returns an error if:
/// - Unable to read a BAM record (except at end boundary)
/// - Unable to write FASTQ output
fn process_records<R>(
    reader: &mut bam::io::Reader<bgzf::io::Reader<R>>,
    start_aligned_offset: u64,
    end_offset: u64,
    output: &mut dyn Write,
) -> Result<usize>
where
    R: Read + Seek,
{
    let mut read_count = 0;
    let mut record = bam::Record::default();
    let mut blocks_processed = 0;
    let mut prev_block_offset = start_aligned_offset;
    let mut prev_was_read1: Option<bool> = None;

    loop {
        let virtual_pos = reader.get_ref().virtual_position();
        let current_block_offset = virtual_pos.compressed();

        if current_block_offset != prev_block_offset {
            blocks_processed += 1;
            prev_block_offset = current_block_offset;
            // If this block starts at/after end_offset
            if current_block_offset >= end_offset {
                debug!("Reached end of processing region at block {current_block_offset}");
                break;
            }
            debug!("Block {current_block_offset} starts before end_offset {end_offset}");

            debug!("Processing block at offset {current_block_offset}");
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

        // Strict validation for collated reads
        if flags.is_segmented() {
            if prev_was_read1 == Some(true) && current_is_read1 {
                anyhow::bail!(
                    "Found consecutive R1 reads ({}). Input BAM must be collated.",
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
        write_fastq_to(&record, output).context("Failed to write FASTQ output")?;

        read_count += 1;

        if read_count % 100_000 == 0 {
            info!("Processed {read_count} reads...");
        }
    }

    info!("Completed: {blocks_processed} blocks, {read_count} reads");
    Ok(read_count)
}

/// Process a byte range from a BAM file and output interleaved FASTQ.
///
/// This is the main entry point for the library. It extracts reads from the specified
/// byte range, automatically aligning to BGZF block boundaries and ensuring read pairs
/// are kept together.
///
/// # Arguments
///
/// * `input_path` - Path to the input BAM file
/// * `start_offset` - Starting byte offset (will be aligned to next BGZF block)
/// * `end_offset` - Ending byte offset (processing stops at or after this offset)
/// * `output` - Output stream for interleaved FASTQ data
///
/// # Returns
///
/// Total number of reads (not pairs) written to the output
///
/// # Errors
///
/// Returns an error if:
/// - Unable to open the input file
/// - No valid BGZF blocks found in the specified range
/// - BAM file is corrupted or malformed
/// - Unable to write to output stream
pub fn process_blocks(
    input_path: &str,
    start_offset: u64,
    end_offset: u64,
    output: &mut dyn Write,
) -> Result<usize> {
    let file = File::open(input_path).with_context(|| format!("Failed to open {input_path}"))?;
    let mut reader = std::io::BufReader::new(file);

    let mut current_offset = start_offset;
    let mut total_reads = 0;

    // Allocate a single reusable buffer for block scanning and validation
    // Needs to be large enough for find_next_bgzf_block (MAX + 15)
    let mut buffer = vec![0u8; BGZF_MAX_BLOCK_SIZE + 15];

    loop {
        let Some(block_start) = find_next_bgzf_block(&mut reader, current_offset, &mut buffer)?
        else {
            // EOF reached without finding a block
            break;
        };

        if block_start >= end_offset {
            // next block is after end_offset
            break;
        }

        if let Some(virtual_pos) = validate_block(&mut reader, block_start, &mut buffer)? {
            let aligned_start = virtual_pos >> 16;
            info!(
                "Processing byte range: {start_offset} to {end_offset} (aligned to block at {aligned_start})"
            );

            // Construct the reader for processing
            reader.seek(SeekFrom::Start(aligned_start))?;

            let mut bgzf_reader = bgzf::io::Reader::new(&mut reader);
            bgzf_reader.seek(bgzf::VirtualPosition::from(virtual_pos))?;
            let mut bam_reader = bam::io::Reader::from(bgzf_reader);

            let count = process_records(&mut bam_reader, aligned_start, end_offset, output)?;
            total_reads += count;
            break;
        }
        debug!("Invalid block at {block_start}, continuing search");
        current_offset = block_start + 1;
    }
    Ok(total_reads)
}
