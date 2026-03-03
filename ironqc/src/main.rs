//! Entry point for the `ironqc` command-line interface.

mod bundle;
mod cli;
mod indexcov;
mod io;
mod mosdepth;
mod stats;

use anyhow::Result;
use clap::Parser;
use cli::{Cli, Commands};

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Stats(args) => stats::run(args),
        Commands::Mosdepth(args) => mosdepth::run(args),
        Commands::Indexcov(args) => indexcov::run(args),
        Commands::Bundle(args) => bundle::run(args),
    }
}
