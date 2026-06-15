use anyhow::{Context, Result};
use clap::Parser;
use std::{
    fs,
    io::{self, BufWriter, Write},
    path::PathBuf,
};

#[derive(Parser)]
#[command(
    name = "merge_fastp_json",
    about = "Merge multiple fastp JSON files from chunked BAM processing into a single result"
)]
struct Args {
    /// Input fastp JSON files to merge
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Output file path (default: stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let data_list = args
        .inputs
        .iter()
        .map(|path| {
            let content = fs::read_to_string(path)
                .with_context(|| format!("Failed to read {}", path.display()))?;
            serde_json::from_str(&content)
                .with_context(|| format!("Failed to parse JSON from {}", path.display()))
        })
        .collect::<Result<Vec<_>>>()?;

    let merged = bamslice::fastp_merge::merge_fastp_jsons(&data_list)
        .context("Failed to merge fastp JSON files")?;
    let json_str = serde_json::to_string_pretty(&merged)?;

    if let Some(ref path) = args.output {
        let file = fs::File::create(path)
            .with_context(|| format!("Failed to create {}", path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "{json_str}")?;
    } else {
        let mut writer = BufWriter::new(io::stdout().lock());
        writeln!(writer, "{json_str}")?;
    }

    Ok(())
}
