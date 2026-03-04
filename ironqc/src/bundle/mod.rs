use crate::cli::BundleArgs;
use crate::indexcov;
use crate::mosdepth::MosdepthAccumulator;
use crate::stats::StatsAccumulator;
use anyhow::{Context, Result};
use noodles::bam;
use std::path::PathBuf;

pub fn run(args: BundleArgs) -> Result<()> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("failed to open BAM {}", args.bam.display()))?;

    let header = reader.read_header()?;

    let mut stats_acc = StatsAccumulator::new();
    let mut mosdepth_acc = MosdepthAccumulator::new(&header, 500, false, false);

    for result in reader.records() {
        let record = result.context("failed to read BAM record")?;
        stats_acc.process_record(&record, &header);
        mosdepth_acc.process_record(&record, &header);
    }

    let finalized = stats_acc.finalize();
    let stats_output = PathBuf::from(format!("{}.stats", args.prefix));
    finalized.write_output(&stats_output)?;

    mosdepth_acc.write_outputs(&args.prefix)?;

    // Extract just the filename stem from the full prefix path for indexcov,
    // since indexcov joins directory + prefix internally.
    let indexcov_prefix = std::path::Path::new(&args.prefix)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(&args.prefix)
        .to_string();

    let indexcov_args = crate::cli::IndexcovArgs {
        bams: vec![args.bam],
        fai: args.fai,
        directory: args.indexcov_dir,
        prefix: Some(indexcov_prefix),
    };
    indexcov::run(indexcov_args)?;

    Ok(())
}
