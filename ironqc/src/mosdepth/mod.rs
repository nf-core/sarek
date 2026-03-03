//! Placeholder implementation of the `mosdepth` subcommand.

use crate::cli::MosdepthArgs;
use crate::io;
use anyhow::Result;
use std::path::PathBuf;

/// Create deterministic placeholder outputs for `mosdepth`.
pub fn run(args: MosdepthArgs) -> Result<()> {
    let prefix = args.prefix;
    let outputs = [
        format!("{prefix}.mosdepth.global.dist.txt"),
        format!("{prefix}.mosdepth.region.dist.txt"),
        format!("{prefix}.mosdepth.summary.txt"),
    ];

    for output in outputs {
        io::touch(PathBuf::from(output))?;
    }

    Ok(())
}
