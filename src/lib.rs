use anyhow::{Context, Result};
use log::{debug, info};
use noodles::bam;
use noodles::bgzf;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};

// BGZF block magic bytes: gzip magic (0x1f 0x8b) + deflate method (0x08) + extra field flag (0x04)
const BGZF_MAGIC: &[u8; 4] = &[0x1f, 0x8b, 0x08, 0x04];

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

/// Find the next BGZF block at or after the given offset
/// # Errors
/// Returns an error if file operations fail or no BGZF block is found within a reasonable range
pub fn find_next_bgzf_block<R>(reader: &mut R, approximate_offset: u64) -> Result<u64>
where
    R: Read + Seek,
{
    const MAX_SCAN: usize = 65536; // Don't scan more than one max BGZF block size

    // Special case: offset 0 is always the start
    if approximate_offset == 0 {
        return Ok(0);
    }

    reader
        .seek(SeekFrom::Start(approximate_offset))
        .with_context(|| format!("Failed to seek to offset {approximate_offset}"))?;

    let mut buffer = vec![0u8; 8192];
    let mut bytes_scanned = 0u64;

    loop {
        let bytes_read = reader
            .read(&mut buffer)
            .context("Failed to read from BAM file while scanning for BGZF block")?;
        if bytes_read == 0 {
            anyhow::bail!("Reached EOF without finding BGZF block");
        }

        // Look for valid BGZF block headers
        if bytes_read >= 16 {
            for i in 0..=(bytes_read - 16) {
                let block_offset = approximate_offset + bytes_scanned + i as u64;
                if is_valid_bgzf_header(&buffer, i) {
                    debug!("Found BGZF block signature at offset {block_offset}");
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

/// Check if a record is read 1 in a paired-end read
/// # Errors
/// Returns an error if the record flags cannot be parsed
pub fn is_read1<R>(record: &R) -> std::io::Result<bool>
where
    R: noodles::sam::alignment::Record,
{
    let flags = record.flags()?;
    if flags.is_segmented() {
        Ok(flags.is_first_segment())
    } else {
        Ok(true) // unpaired reads are treated as read 1
    }
}

/// Write a read directly to output in FASTQ format (interleaved - same format as samtools fastq)
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

/// Scan a BGZF block to find the offset of the first valid BAM record
fn find_first_record_in_block<R>(reader: &mut R, block_start: u64) -> Result<Option<u64>>
where
    R: Read + Seek,
{
    reader.seek(SeekFrom::Start(block_start))?;
    let mut reader = flate2::read::GzDecoder::new(reader);

    // Read enough data to find a record. 64KB is typical max block size.
    let mut buffer = vec![0u8; 65536];
    let bytes_read = std::io::Read::read(&mut reader, &mut buffer)?;

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

/// Checks if a record is valid according to the (BAM spec)[<https://samtools.github.io/hts-specs/SAMv1.pdf>] (Section 4.2)
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

/// Attempt to open a valid BGZF block at the given offset.
///
/// Returns `Ok(Some((reader, aligned_start)))` if successful.
/// Returns `Ok(None)` if the block is invalid or no record found.
/// # Errors
/// Returns an error if unable to open the file or seek to the desired position
/// Attempt to validate a BGZF block at the given offset and find the first record.
///
/// Returns `Ok(Some(virtual_position))` if successful.
/// Returns `Ok(None)` if the block is invalid or no record found.
/// # Errors
/// Returns an error if unable to seek or read
pub fn validate_block<R>(reader: &mut R, aligned_start: u64) -> Result<Option<u64>>
where
    R: Read + Seek,
{
    reader.seek(SeekFrom::Start(aligned_start))?;

    // Try to find a record start within the decompressed block
    match find_first_record_in_block(reader, aligned_start) {
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

/// Reads records from a BAM file starting at a specific offset and writes tem to the output stream
///
/// # Arguments
/// * `reader` - A mutable reference to a BAM reader
/// * `start_aligned_offset` - The start offset of the aligned region
/// * `end_offset` - The end offset of the aligned region
/// * `output` - A mutable reference to a writeable stream
///
/// # Returns
/// * `Ok(usize)` - The number of records processed
/// * `Err(anyhow::Error)` - An error if the reader fails to read a record
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
    let mut first_read_in_range = true;
    let mut prev_was_read1 = false;

    //
    loop {
        let virtual_pos = reader.get_ref().virtual_position();
        let current_block_offset = virtual_pos.compressed();

        if current_block_offset != prev_block_offset {
            blocks_processed += 1;
            prev_block_offset = current_block_offset;
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

        // bam::Record::flags() returns Flags directly
        let is_paired = record.flags().is_segmented();
        // is_read1 returns Result<bool>
        let current_is_read1 = is_read1(&record)?;

        if first_read_in_range {
            first_read_in_range = false;
            if is_paired {
                let current_is_read2 = !current_is_read1;
                if current_is_read2 {
                    debug!("Skipping first read (read 2) to align on read pairs");
                    continue;
                }
            }
        }

        let past_end_offset = current_block_offset >= end_offset;
        let need_mate_for_read1 = is_paired && prev_was_read1;

        if past_end_offset && !need_mate_for_read1 {
            debug!("Reached end offset at block {current_block_offset}");
            break;
        }

        write_fastq_to(&record, output).context("Failed to write FASTQ output")?;

        read_count += 1;
        prev_was_read1 = current_is_read1;

        let current_is_read2 = !current_is_read1;
        if past_end_offset && is_paired && current_is_read2 {
            debug!("Wrote mate read 2 from block {current_block_offset}, stopping");
            break;
        }

        if read_count % 100_000 == 0 {
            info!("Processed {read_count} reads...");
        }
    }

    info!("Completed: {blocks_processed} blocks, {read_count} reads");
    Ok(read_count)
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
    let file = File::open(input_path).with_context(|| format!("Failed to open {input_path}"))?;
    let mut reader = std::io::BufReader::new(file);

    let mut current_offset = start_offset;
    let mut total_reads = 0;

    loop {
        let block_start = find_next_bgzf_block(&mut reader, current_offset)?;
        if block_start >= end_offset {
            break;
        }

        if let Some(virtual_pos) = validate_block(&mut reader, block_start)? {
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

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam::alignment::RecordBuf;
    use noodles::sam::alignment::record::Flags;

    #[test]
    fn test_is_read1_unpaired() {
        let record = RecordBuf::default();
        // Unpaired reads should be treated as read 1
        assert!(is_read1(&record).unwrap(), "Unpaired read should be read 1");
    }

    #[test]
    fn test_is_read1_paired_first() {
        let flags = Flags::SEGMENTED | Flags::FIRST_SEGMENT;
        let record = RecordBuf::builder().set_flags(flags).build();

        assert!(
            is_read1(&record).unwrap(),
            "Paired read marked as first in template should be read 1"
        );
    }

    #[test]
    fn test_is_read1_paired_last() {
        let flags = Flags::SEGMENTED | Flags::LAST_SEGMENT;
        let record = RecordBuf::builder().set_flags(flags).build();

        assert!(
            !is_read1(&record).unwrap(),
            "Paired read marked as last in template should NOT be read 1"
        );
    }
}
