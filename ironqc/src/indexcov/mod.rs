//! Placeholder implementation of the `indexcov` subcommand.

use crate::cli::IndexcovArgs;
use crate::io;
use anyhow::Result;
use std::path::PathBuf;

/// Create deterministic placeholder outputs for `indexcov`.
pub fn run(args: IndexcovArgs) -> Result<()> {
    io::ensure_dir(&args.directory)?;

    let prefix = args.prefix.unwrap_or_else(|| "indexcov".to_string());
    let base = args.directory.join(format!("{prefix}-indexcov"));
    let outputs = [
        PathBuf::from(format!("{}.ped", base.display())),
        PathBuf::from(format!("{}.roc", base.display())),
        PathBuf::from(format!("{}.bed.gz", base.display())),
        PathBuf::from(format!("{}.html", base.display())),
        PathBuf::from(format!("{}.sex.png", base.display())),
        PathBuf::from(format!("{}.roc.png", base.display())),
    ];

    for output in outputs {
        io::touch(output)?;
    }

    Ok(())
}
