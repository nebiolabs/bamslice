use std::collections::HashSet;
use std::fs;

const TEST_BAM: &str = "tests/fixtures/illumina-test1.bam";
const TEST_ONT_BAM: &str = "tests/fixtures/ont-test1.bam";
const EXPECTED_READS: usize = 9534; // Pre-calculated: samtools fastq emseq-test1.bam | wc -l / 4

fn count_fastq_reads(content: &str) -> usize {
    content.lines().count() / 4
}

fn get_file_size(path: &str) -> u64 {
    fs::metadata(path).unwrap().len()
}

/// Extract read names from FASTQ content
fn extract_read_names(fastq_content: &str) -> HashSet<String> {
    fastq_content
        .lines()
        .enumerate()
        .filter_map(|(i, line)| {
            // Every 4th line starting at 0 is a read name (@ line)
            if i % 4 == 0 && line.starts_with('@') {
                Some(line[1..].to_string()) // Remove '@' prefix
            } else {
                None
            }
        })
        .collect()
}

#[test]
fn test_split_processing_produces_all_reads_no_dupes() {
    // Check that test file exists
    assert!(
        std::path::Path::new(TEST_BAM).exists(),
        "Test BAM file not found: {TEST_BAM}"
    );

    let file_size = get_file_size(TEST_BAM);
    let half = file_size / 2;

    eprintln!("File size: {file_size} bytes");
    eprintln!("Split point: {half} bytes");

    // Process first half in memory
    let mut chunk1_buffer = Vec::new();
    let chunk1_reads = bamslice::process_blocks(TEST_BAM, 0, half, &mut chunk1_buffer).unwrap();

    // Process second half in memory
    let mut chunk2_buffer = Vec::new();
    let chunk2_reads =
        bamslice::process_blocks(TEST_BAM, half, file_size, &mut chunk2_buffer).unwrap();

    // Convert to strings
    let chunk1_content = String::from_utf8(chunk1_buffer).unwrap();
    let chunk2_content = String::from_utf8(chunk2_buffer).unwrap();

    // Extract read names from both chunks
    let chunk1_names = extract_read_names(&chunk1_content);
    let chunk2_names = extract_read_names(&chunk2_content);

    // Find any duplicate read names (intersection of the two sets)
    let duplicates: Vec<_> = chunk1_names.intersection(&chunk2_names).collect();

    assert!(
        duplicates.is_empty(),
        "Found {} duplicate read(s) appearing in both chunks: {:?}",
        duplicates.len(),
        duplicates.iter().take(10).collect::<Vec<_>>() // Show first 10 duplicates
    );

    let chunk1_fastq_reads = count_fastq_reads(&chunk1_content);
    let chunk2_fastq_reads = count_fastq_reads(&chunk2_content);
    let total_reads = chunk1_fastq_reads + chunk2_fastq_reads;

    // Verify that return values match FASTQ content
    assert_eq!(
        chunk1_reads, chunk1_fastq_reads,
        "Chunk 1: returned count doesn't match FASTQ content"
    );
    assert_eq!(
        chunk2_reads, chunk2_fastq_reads,
        "Chunk 2: returned count doesn't match FASTQ content"
    );

    // Verify total matches expected from samtools
    assert_eq!(
        total_reads, EXPECTED_READS,
        "Total read count mismatch: {total_reads} (ours) vs {EXPECTED_READS} (expected from samtools)"
    );

    // Verify both chunks have reads (no empty split)
    assert!(chunk1_reads > 0, "Chunk 1 should have reads");
    assert!(chunk2_reads > 0, "Chunk 2 should have reads");

    println!("  Chunk 1: {chunk1_reads} reads");
    println!("  Chunk 2: {chunk2_reads} reads");
}

#[test]
fn test_whole_file_processing() {
    let file_size = get_file_size(TEST_BAM);

    // Test processing the entire file in one go, in memory
    let mut buffer = Vec::new();
    let read_count = bamslice::process_blocks(TEST_BAM, 0, file_size, &mut buffer).unwrap();

    let content = String::from_utf8(buffer).unwrap();
    let fastq_reads = count_fastq_reads(&content);

    assert_eq!(
        read_count, fastq_reads,
        "Returned read count doesn't match FASTQ content {read_count} vs {fastq_reads}"
    );

    assert_eq!(
        fastq_reads, EXPECTED_READS,
        "FASTQ record count mismatch: {fastq_reads} vs {EXPECTED_READS}"
    );
}

#[test]
fn test_first_read_content() {
    let file_size = get_file_size(TEST_BAM);

    // Process the entire file
    let mut buffer = Vec::new();
    bamslice::process_blocks(TEST_BAM, 0, file_size, &mut buffer).unwrap();

    let content = String::from_utf8(buffer).unwrap();
    let mut lines = content.lines();

    // Expected first read (with barcode in Illumina format: read:filtered:control:barcode)
    let expected_name = "@NS500355:NS500355:HHVN5AFX7:1:11101:11181:6634/1 1:N:0:CGTCAAGA-GGGTTGTT";
    let expected_seq =
        "AAATTTTAGAAAATTGTTATATTATTTGGGTTATTAGTGGAGATATTTGTTATAATTTTTTTTTAGGCGTAATTTG";
    let expected_qual =
        "AAAAAEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";

    // Read first FASTQ record (4 lines)
    let name = lines.next().expect("Missing read name line");
    let seq = lines.next().expect("Missing sequence line");
    let plus = lines.next().expect("Missing + line");
    let qual = lines.next().expect("Missing quality line");

    assert_eq!(name, expected_name, "First read name mismatch");
    assert_eq!(seq, expected_seq, "First read sequence mismatch");
    assert_eq!(plus, "+", "Expected + separator line");
    assert_eq!(qual, expected_qual, "First read quality scores mismatch");
}

#[test]
fn test_blocks_start_with_read1() {
    let file_size = get_file_size(TEST_BAM);

    // Test multiple split points to ensure blocks always start with read 1
    let split_points = [
        file_size / 4,
        file_size / 3,
        file_size / 2,
        file_size * 2 / 3,
        file_size * 3 / 4,
    ];

    for &split_point in &split_points {
        let mut buffer = Vec::new();
        bamslice::process_blocks(TEST_BAM, split_point, file_size, &mut buffer).unwrap();

        let content = String::from_utf8(buffer).unwrap();
        let mut lines = content.lines();

        if let Some(first_line) = lines.next() {
            // First line should be a read name starting with @
            assert!(
                first_line.starts_with('@'),
                "First line should be read name at split {split_point}"
            );

            // Check if it's read 1 by looking for /1 suffix
            assert!(
                first_line.contains("/1 "),
                "First read should be read 1 (have /1 suffix) at split {split_point}, got: {first_line}"
            );
        }
    }
}

#[test]
fn test_interleaved_pairs_in_output() {
    let file_size = get_file_size(TEST_BAM);
    let half = file_size / 2;

    // Process second half
    let mut buffer = Vec::new();
    bamslice::process_blocks(TEST_BAM, half, file_size, &mut buffer).unwrap();

    let content = String::from_utf8(buffer).unwrap();
    let lines: Vec<&str> = content.lines().collect();

    // Check that reads alternate between /1 and /2
    let mut expect_read1 = true;
    for i in (0..lines.len()).step_by(4) {
        if i >= lines.len() {
            break;
        }

        let read_name = lines[i];
        if expect_read1 {
            assert!(
                read_name.contains("/1 "),
                "Expected read 1 at line {i}, got: {read_name}"
            );
        } else {
            assert!(
                read_name.contains("/2 "),
                "Expected read 2 at line {i}, got: {read_name}"
            );
        }
        expect_read1 = !expect_read1;
    }
}

#[test]
fn test_complete_pairs_at_boundaries() {
    let file_size = get_file_size(TEST_BAM);
    let split = file_size / 2;

    // Process first chunk
    let mut chunk1 = Vec::new();
    let chunk1_count = bamslice::process_blocks(TEST_BAM, 0, split, &mut chunk1).unwrap();

    // Process second chunk
    let mut chunk2 = Vec::new();
    let chunk2_count = bamslice::process_blocks(TEST_BAM, split, file_size, &mut chunk2).unwrap();

    // Both chunks should have even number of reads (complete pairs)
    assert_eq!(
        chunk1_count % 2,
        0,
        "First chunk should have complete pairs (even count), got {chunk1_count}"
    );
    assert_eq!(
        chunk2_count % 2,
        0,
        "Second chunk should have complete pairs (even count), got {chunk2_count}"
    );
}

#[test]
fn test_ont_reads() {
    let split = 180;
    // Process first chunk
    let mut chunk1 = Vec::new();
    let chunk1_count = bamslice::process_blocks(TEST_ONT_BAM, 0, split, &mut chunk1).unwrap();

    // Process second chunk
    let mut chunk2 = Vec::new();
    let chunk2_count =
        bamslice::process_blocks(TEST_ONT_BAM, split, 1_000_000, &mut chunk2).unwrap();

    assert_eq!(chunk1_count, 84);
    assert_eq!(chunk2_count, 125 - 84);
}

#[test]
fn test_skips_false_positive_gzip_magic() {
    // Offset 132321 contains bytes "1f 8b 08" which match gzip magic but is not a valid
    // BGZF block (4th byte is 0x40 not 0x04). This test verifies we correctly skip it
    // and find the next valid BGZF block.
    let false_positive_offset = 132_321;
    let end_offset = 200_000;

    let mut buffer = Vec::new();
    let result = bamslice::process_blocks(TEST_BAM, false_positive_offset, end_offset, &mut buffer);
    assert!(
        result.is_ok(),
        "Should skip false positive gzip magic at {} and find valid block, got error: {:?}",
        false_positive_offset,
        result.err()
    );
}

#[test]
fn test_dnbseq_seek_issue() {
    // This test reproduces the issue where seeking to offset 1000 in dnbseq-test1.bam
    // lands in the middle of a record or finds a false positive BGZF header,
    // requiring robust retry logic to find the next valid block.
    const TEST_FILE: &str = "tests/fixtures/dnbseq-test1.bam";

    // Parameters from the failing command
    let start_offset = 1000;
    let end_offset = 100_000;

    let mut buffer = Vec::new();
    let result = bamslice::process_blocks(TEST_FILE, start_offset, end_offset, &mut buffer);

    assert!(
        result.is_ok(),
        "Should successfully process blocks even with tricky seek, got error: {:?}",
        result.err()
    );

    let read_count = result.unwrap();
    assert!(
        read_count > 0,
        "Should have extracted reads, got {read_count}"
    );
}

#[test]
fn test_process_near_eof() {
    const TEST_FILE: &str = "tests/fixtures/dnbseq-test1.bam";

    let mut output = Vec::new();

    // Start scanning near the end of the file (file is ~13.6MB)
    // This should still find the last block and extract remaining reads
    let result = bamslice::process_blocks(TEST_FILE, 13_590_000, 13_700_000, &mut output);

    assert!(result.is_ok(), "Should handle near-EOF offsets gracefully");
    let read_count = result.unwrap();

    // Should get the last 147 reads
    assert_eq!(read_count, 147, "Should extract remaining reads near EOF");

    // Verify we got the last read in the file
    let output_str = String::from_utf8(output).unwrap();
    assert!(
        output_str.contains("V350145865L1C001R01500464560"),
        "Should include the last read from the file"
    );
}

#[test]
fn test_window_after_last_block() {
    const TEST_FILE: &str = "tests/fixtures/dnbseq-test1.bam";

    let mut output = Vec::new();

    // Start beyond all BGZF blocks. 13_598_802 finds 146 reads,so we start just after that
    let result = bamslice::process_blocks(TEST_FILE, 13_598_803, u64::MAX, &mut output);

    assert!(
        result.is_ok(),
        "Should not error when starting beyond last block"
    );
    let read_count = result.unwrap();

    assert_eq!(read_count, 0, "Should return 0 reads when no blocks found");
}

#[test]
fn test_window_after_end_of_file() {
    const TEST_FILE: &str = "tests/fixtures/dnbseq-test1.bam";

    let mut output = Vec::new();

    // Start beyond all BGZF blocks. 13_598_802 finds 146 reads,so we start just after that
    let result = bamslice::process_blocks(TEST_FILE, 18_000_000, u64::MAX, &mut output);

    assert!(
        result.is_ok(),
        "Should not error when starting beyond last block"
    );
    let read_count = result.unwrap();

    assert_eq!(read_count, 0, "Should return 0 reads when no blocks found");
}

#[test]
fn test_dnbseq_100_chunks_exactly_200k_reads() {
    // This test confirms that processing dnbseq-test1.bam in 100000-byte chunks
    // produces exactly 200,000 reads with the same read names as samtools
    // Truth hash: samtools fastq dnbseq-test1.bam | awk 'NR % 4 == 1' | cut -d' ' -f1 | sed 's/^@//' | md5
    const TEST_FILE: &str = "tests/fixtures/dnbseq-test1.bam";
    const EXPECTED_READS: usize = 200_000;
    const TRUTH_HASH: &str = "25fe990b5612285e8b0fa04e2bf0e612"; // MD5 hash from samtools
    const CHUNK_SIZE: u64 = 100_000;

    let file_size = get_file_size(TEST_FILE);

    let mut total_reads = 0;
    let mut hasher = md5::Context::new();
    let mut chunk_num = 0;

    let mut start = 0;
    while start < file_size {
        let end = (start + CHUNK_SIZE).min(file_size);

        let mut buffer = Vec::new();
        let chunk_reads = bamslice::process_blocks(TEST_FILE, start, end, &mut buffer)
            .unwrap_or_else(|e| panic!("Failed to process chunk {chunk_num} ({start}-{end}): {e}"));

        total_reads += chunk_reads;

        // Hash read names from this chunk
        let content = String::from_utf8(buffer).unwrap();
        for (idx, line) in content.lines().enumerate() {
            if idx % 4 == 0 && line.starts_with('@') {
                let read_name = &line[1..];
                let core_name = read_name.split_whitespace().next().unwrap_or(read_name);
                hasher.consume(core_name.as_bytes());
                hasher.consume(b"\n"); // Match the newline from awk output
            }
        }

        chunk_num += 1;
        if chunk_num % 10 == 0 {
            println!("  Processed chunk {chunk_num}: {total_reads} reads so far");
        }

        start = end;
    }

    println!("  Processed {chunk_num} total chunks");

    let result_hash = format!("{:x}", hasher.compute());

    assert_eq!(
        total_reads, EXPECTED_READS,
        "Total reads from all chunks should be exactly 200,000, got {total_reads}"
    );

    assert_eq!(
        result_hash, TRUTH_HASH,
        "MD5 hash of read names doesn't match samtools output.\nGot:      {result_hash}\nExpected: {TRUTH_HASH}"
    );

    println!(
        "  ✓ {chunk_num} chunks (100KB each) produced exactly {total_reads} reads with correct hash: {result_hash}"
    );
}
