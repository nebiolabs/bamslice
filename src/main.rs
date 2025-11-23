use anyhow::{Context, Result};
use clap::Parser;
use log::info;

/// Extract specific BGZF blocks from BAM files and convert to interleaved FASTQ
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input BAM file
    #[arg(short, long)]
    input: String,

    /// Starting byte offset (will find next BGZF block at or after this offset)
    #[arg(short = 's', long)]
    start_offset: u64,

    /// Ending byte offset (will process until reaching a block at or after this offset)
    #[arg(short = 'e', long)]
    end_offset: u64,

    /// Output file (default: stdout)
    #[arg(short, long)]
    output: Option<String>,

    /// Log level (off, error, warn, info, debug, trace)
    #[arg(short = 'l', long, default_value = "info")]
    log_level: String,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Initialize logger with specified level (can be overridden by RUST_LOG env var)
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or(&args.log_level))
        .init();

    // Validate arguments
    if args.start_offset >= args.end_offset {
        anyhow::bail!(
            "start_offset ({}) must be less than end_offset ({})",
            args.start_offset,
            args.end_offset
        );
    }

    // Uses Box<dyn Write> to handle both File and Stdout via dynamic dispatch.
    // This simplifies the code compared to generics, at the cost of slight runtime overhead.
    let mut writer: Box<dyn std::io::Write> = if let Some(output_path) = &args.output {
        let file = std::fs::File::create(output_path)
            .with_context(|| format!("Failed to create output file {output_path}"))?;
        Box::new(std::io::BufWriter::new(file))
    } else {
        let stdout = std::io::stdout();
        Box::new(std::io::BufWriter::new(stdout.lock()))
    };

    let read_count =
        bamslice::process_blocks(&args.input, args.start_offset, args.end_offset, &mut writer)?;

    info!("Total reads extracted: {read_count}");
    Ok(())
}
