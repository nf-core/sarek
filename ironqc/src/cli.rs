//! Clap argument definitions for `ironqc` subcommands.

use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

/// Top-level CLI parser.
#[derive(Debug, Parser)]
#[command(name = "ironqc")]
#[command(about = "Unified QC scaffold for samtools stats, mosdepth, and indexcov")]
pub struct Cli {
    /// Subcommand to execute.
    #[command(subcommand)]
    pub command: Commands,
}

/// All supported `ironqc` subcommands.
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Placeholder for `samtools stats`-compatible output plumbing.
    Stats(StatsArgs),
    /// Placeholder for `mosdepth`-compatible output plumbing.
    Mosdepth(MosdepthArgs),
    /// Placeholder for `goleft indexcov`-compatible output plumbing.
    Indexcov(IndexcovArgs),
    /// Run all three placeholder tools using one command.
    Bundle(BundleArgs),
}

/// Arguments for the `stats` subcommand.
#[derive(Debug, Args)]
pub struct StatsArgs {
    /// Input BAM/CRAM file.
    pub bam: PathBuf,
    /// Reference FASTA used for CRAM decoding and stats parity.
    #[arg(long)]
    pub reference: PathBuf,
    /// Number of worker threads.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
    /// Output prefix used to create `<prefix>.stats`.
    #[arg(long)]
    pub prefix: String,
}

/// Arguments for the `mosdepth` subcommand.
#[derive(Debug, Args)]
pub struct MosdepthArgs {
    /// Input BAM/CRAM file.
    pub bam: PathBuf,
    /// Reference FASTA used by mosdepth.
    #[arg(long)]
    pub fasta: PathBuf,
    /// Number of worker threads.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
    /// Bin size or BED path.
    #[arg(long)]
    pub by: Option<String>,
    /// Disable per-base output (`-n` in upstream mosdepth).
    #[arg(short = 'n')]
    pub no_per_base: bool,
    /// Enable fast mode.
    #[arg(long)]
    pub fast_mode: bool,
    /// Output prefix used to create `<prefix>.mosdepth.*` outputs.
    pub prefix: String,
}

/// Arguments for the `indexcov` subcommand.
#[derive(Debug, Args)]
pub struct IndexcovArgs {
    /// One or more BAM files.
    #[arg(required = true)]
    pub bams: Vec<PathBuf>,
    /// Reference FAI path.
    #[arg(long)]
    pub fai: PathBuf,
    /// Output directory for indexcov artifacts.
    #[arg(long)]
    pub directory: PathBuf,
    /// Output prefix for `<prefix>-indexcov.*` files.
    #[arg(long)]
    pub prefix: Option<String>,
}

/// Arguments for the `bundle` subcommand.
#[derive(Debug, Args)]
pub struct BundleArgs {
    /// Input BAM/CRAM file.
    pub bam: PathBuf,
    /// Reference FASTA for stats and mosdepth paths.
    #[arg(long)]
    pub reference: PathBuf,
    /// Reference FAI for indexcov path.
    #[arg(long)]
    pub fai: PathBuf,
    /// Number of worker threads.
    #[arg(long, default_value_t = 1)]
    pub threads: usize,
    /// Output prefix for stats and mosdepth outputs.
    #[arg(long)]
    pub prefix: String,
    /// Output directory for indexcov artifacts.
    #[arg(long)]
    pub indexcov_dir: PathBuf,
}
