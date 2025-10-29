use std::collections::HashSet;
use std::fs;

const TEST_BAM: &str = "tests/fixtures/emseq-test1.bam";
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
