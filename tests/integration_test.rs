use std::fs;

const TEST_BAM: &str = "tests/fixtures/emseq-test1.bam";
const EXPECTED_READS: usize = 9534; // Pre-calculated: samtools fastq emseq-test1.bam | wc -l / 4

fn count_fastq_reads(content: &str) -> usize {
    content.lines().count() / 4
}

fn get_file_size(path: &str) -> u64 {
    fs::metadata(path).unwrap().len()
}

#[test]
fn test_split_processing_produces_all_reads() {
    // Check that test file exists
    assert!(
        std::path::Path::new(TEST_BAM).exists(),
        "Test BAM file not found: {}",
        TEST_BAM
    );

    let file_size = get_file_size(TEST_BAM);
    let half = file_size / 2;

    eprintln!("File size: {} bytes", file_size);
    eprintln!("Split point: {} bytes", half);

    // Process first half in memory
    let mut chunk1_buffer = Vec::new();
    let chunk1_reads =
        bam_block_extractor::process_blocks(TEST_BAM, 0, half, None, &mut chunk1_buffer).unwrap();

    // Process second half in memory
    let mut chunk2_buffer = Vec::new();
    let chunk2_reads =
        bam_block_extractor::process_blocks(TEST_BAM, half, file_size, None, &mut chunk2_buffer)
            .unwrap();

    // Convert to strings
    let chunk1_content = String::from_utf8(chunk1_buffer).unwrap();
    let chunk2_content = String::from_utf8(chunk2_buffer).unwrap();

    let chunk1_fastq_reads = count_fastq_reads(&chunk1_content);
    let chunk2_fastq_reads = count_fastq_reads(&chunk2_content);
    let total_reads = chunk1_fastq_reads + chunk2_fastq_reads;

    eprintln!("Chunk 1 reads (returned): {}", chunk1_reads);
    eprintln!("Chunk 1 reads (FASTQ): {}", chunk1_fastq_reads);
    eprintln!("Chunk 2 reads (returned): {}", chunk2_reads);
    eprintln!("Chunk 2 reads (FASTQ): {}", chunk2_fastq_reads);
    eprintln!("Total reads: {}", total_reads);
    eprintln!("Expected reads: {}", EXPECTED_READS);

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
        "Total read count mismatch: {} (ours) vs {} (expected from samtools)",
        total_reads, EXPECTED_READS
    );

    // Verify both chunks have reads (no empty split)
    assert!(chunk1_reads > 0, "Chunk 1 should have reads");
    assert!(chunk2_reads > 0, "Chunk 2 should have reads");

    println!(
        "✓ Test passed: split processing produces all {} reads",
        EXPECTED_READS
    );
    println!("  Chunk 1: {} reads", chunk1_reads);
    println!("  Chunk 2: {} reads", chunk2_reads);
}

#[test]
fn test_whole_file_processing() {
    let file_size = get_file_size(TEST_BAM);

    // Test processing the entire file in one go, in memory
    let mut buffer = Vec::new();
    let read_count =
        bam_block_extractor::process_blocks(TEST_BAM, 0, file_size, None, &mut buffer).unwrap();

    let content = String::from_utf8(buffer).unwrap();
    let fastq_reads = count_fastq_reads(&content);

    assert_eq!(
        read_count, EXPECTED_READS,
        "Whole file read count mismatch: {} vs {}",
        read_count, EXPECTED_READS
    );

    assert_eq!(
        fastq_reads, EXPECTED_READS,
        "FASTQ record count mismatch: {} vs {}",
        fastq_reads, EXPECTED_READS
    );

    println!(
        "✓ Test passed: whole file processing produces {} reads",
        EXPECTED_READS
    );
}
