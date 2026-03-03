//! Placeholder implementation of the `stats` subcommand.

use crate::cli::StatsArgs;
use crate::io;
use anyhow::Result;
use std::path::PathBuf;

/// Create deterministic placeholder output for `stats`.
pub fn run(args: StatsArgs) -> Result<()> {
    let output = PathBuf::from(format!("{}.stats", args.prefix));
    io::touch(output)
}
