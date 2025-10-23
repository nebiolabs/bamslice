#!/usr/bin/env nextflow

/*
 * bamslice - Parallel BAM to FASTQ conversion using byte ranges
 *
 * This example pipeline demonstrates how to use bamslice to convert
 * BAM files to FASTQ in parallel by splitting on byte ranges, then
 * piping through fastp for QC/filtering.
 */

params.bam = "tests/fixtures/emseq-test1.bam"  // Can be a glob pattern like "*.bam"
params.chunk_size = "100 MB"  // Supports: 100 KB, 100 MB, 1 GB, etc.
params.outdir = "results"
params.bamslice = "${projectDir}/target/release/bamslice"  // Path to bamslice binary

workflow {
    def bamslice_bin = file(params.bamslice)
    if (!bamslice_bin.exists()) {
        error "ERROR: bamslice binary not found at: ${params.bamslice}\n" +
              "Please build it with: cargo build --release\n" +
              "Or specify the correct path with: --bamslice /path/to/bamslice"
    }

    if ("which fastp".execute().waitFor() != 0) {
        log.error "ERROR: fastp not found on PATH"
    }

    // Create channel of BAM files
    bam_ch = Channel.fromPath(params.bam, checkIfExists: true)

    chunk_size_bytes = new nextflow.util.MemoryUnit(params.chunk_size).toBytes()
    println "Chunk size: ${params.chunk_size} (${chunk_size_bytes} bytes)"

    // Create chunks for each BAM file
    chunks_ch = bam_ch
        .map { bam ->
            def file_size = bam.size()
            def num_chunks = Math.ceil(file_size / chunk_size_bytes).intValue()
            println "Input BAM: ${bam} (size: ${file_size} bytes, ${num_chunks} chunks)"

            // Return list of [bam, start, end] tuples for each chunk
            (0..<num_chunks).collect { i ->
                def start = i * chunk_size_bytes
                def end = Math.min(start + chunk_size_bytes, file_size)
                [bam, start, end]
            }
        }
        .flatMap()  // Flatten the list of chunks into individual items

    // Process each chunk: bamslice -> fastp -> count (or align)
    results = processChunk(chunks_ch).view()


}

process processChunk {
    publishDir params.outdir, mode: 'copy', pattern: "*.{json,html}"

    input:
    tuple path(bam), val(start), val(end)

    output:
    stdout emit: counts

    script:
    """
    # Extract chunk and pipe through fastp for QC/filtering
    # Then count reads with wc -l
    # NOTE: Replace 'wc -l' with an aligner (e.g., 'bwa mem ref.fa - | samtools ...')
    #       or other downstream tool that accepts interleaved FASTQ on stdin
    ${params.bamslice} -i ${bam} -s ${start} -e ${end} \
        | fastp --stdin --interleaved_in \
                --stdout \
                --disable_quality_filtering \
                --disable_length_filtering \
                --json ${bam.simpleName}_chunk_${start}.json \
                --html ${bam.simpleName}_chunk_${start}.html \
        | wc -l \
        | awk '{printf "${bam.simpleName} chunk ${start}: %d reads\\n", \$1/4}'
    """
}
