use anyhow::{Context, Result};
use bam_block_extractor::process_blocks;
use clap::Parser;

/// Extract specific BGZF blocks from BAM/CRAM files and convert to interleaved FASTQ
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input BAM or CRAM file
    #[arg(short, long)]
    input: String,

    /// Starting byte offset (will find next BGZF block at or after this offset)
    #[arg(short = 's', long)]
    start_offset: u64,

    /// Ending byte offset (will process until reaching a block at or after this offset)
    #[arg(short = 'e', long)]
    end_offset: u64,

    /// Reference file for CRAM (optional)
    #[arg(short, long)]
    reference: Option<String>,

    /// Output file (default: stdout)
    #[arg(short, long)]
    output: Option<String>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Validate arguments
    if args.start_offset >= args.end_offset {
        anyhow::bail!(
            "start_offset ({}) must be less than end_offset ({})",
            args.start_offset,
            args.end_offset
        );
    }

    // Process blocks and output FASTQ
    let read_count = if let Some(output_path) = &args.output {
        let mut file = std::fs::File::create(output_path)
            .context(format!("Failed to create output file {}", output_path))?;
        process_blocks(
            &args.input,
            args.start_offset,
            args.end_offset,
            args.reference.as_deref(),
            &mut file,
        )?
    } else {
        let stdout = std::io::stdout();
        let mut handle = stdout.lock();
        process_blocks(
            &args.input,
            args.start_offset,
            args.end_offset,
            args.reference.as_deref(),
            &mut handle,
        )?
    };

    eprintln!("Total reads extracted: {}", read_count);
    Ok(())
}
