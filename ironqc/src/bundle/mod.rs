//! Placeholder implementation of the `bundle` subcommand.

use crate::cli::{BundleArgs, IndexcovArgs, MosdepthArgs, StatsArgs};
use crate::{indexcov, mosdepth, stats};
use anyhow::Result;

/// Run all placeholder subcommands and create their expected outputs.
pub fn run(args: BundleArgs) -> Result<()> {
    let stats_args = StatsArgs {
        bam: args.bam.clone(),
        reference: args.reference.clone(),
        threads: args.threads,
        prefix: args.prefix.clone(),
    };
    stats::run(stats_args)?;

    let mosdepth_args = MosdepthArgs {
        bam: args.bam.clone(),
        fasta: args.reference.clone(),
        threads: args.threads,
        by: None,
        no_per_base: false,
        fast_mode: false,
        prefix: args.prefix.clone(),
    };
    mosdepth::run(mosdepth_args)?;

    let indexcov_args = IndexcovArgs {
        bams: vec![args.bam],
        fai: args.fai,
        directory: args.indexcov_dir,
        prefix: Some(args.prefix),
    };
    indexcov::run(indexcov_args)
}
