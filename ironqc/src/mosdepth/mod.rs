use crate::cli::MosdepthArgs;
use anyhow::{Context, Result};
use indexmap::IndexMap;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

struct ContigStats {
    length: u64,
    bases_covered: u64,
    min_depth: u64,
    max_depth: u64,
    _depth_histogram: Vec<u64>,
}

pub struct MosdepthAccumulator {
    contigs: IndexMap<String, ContigStats>,
    window_size: u64,
    fast_mode: bool,
}

impl MosdepthAccumulator {
    pub fn new(
        header: &noodles::sam::Header,
        window_size: u64,
        fast_mode: bool,
        _no_per_base: bool,
    ) -> Self {
        let mut contigs = IndexMap::new();
        for (name, map) in header.reference_sequences() {
            let length = map.length().get() as u64;
            contigs.insert(
                name.to_string(),
                ContigStats {
                    length,
                    bases_covered: 0,
                    min_depth: u64::MAX,
                    max_depth: 0,
                    _depth_histogram: vec![0u64; 1001],
                },
            );
        }
        Self {
            contigs,
            window_size,
            fast_mode,
        }
    }

    pub fn process_record(&mut self, record: &bam::Record, header: &noodles::sam::Header) {
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            return;
        }

        let tid = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => return,
        };

        let contig_name = match header
            .reference_sequences()
            .get_index(tid)
            .map(|(name, _)| name.to_string())
        {
            Some(name) => name,
            None => return,
        };

        let pos = match record.alignment_start() {
            Some(Ok(p)) => p.get() as u64 - 1,
            _ => return,
        };

        let contig = match self.contigs.get_mut(&contig_name) {
            Some(c) => c,
            None => return,
        };

        if self.fast_mode {
            let read_len = record.sequence().len() as u64;
            let end = (pos + read_len).min(contig.length);
            let start_bin = pos / self.window_size;
            let end_bin = end.saturating_sub(1) / self.window_size;
            for _bin in start_bin..=end_bin {
                contig.bases_covered += self.window_size.min(end - pos);
            }
        } else {
            let cigar = record.cigar();
            let mut ref_pos = pos;
            for op in cigar.iter().flatten() {
                let len: u64 = op.len() as u64;
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        contig.bases_covered += len;
                        ref_pos += len;
                    }
                    Kind::Deletion | Kind::Skip => {
                        ref_pos += len;
                    }
                    Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
                }
            }
            let _ = ref_pos;
        }
    }

    pub fn write_outputs(&self, prefix: &str) -> Result<()> {
        self.write_summary(prefix)?;
        self.write_global_dist(prefix)?;
        self.write_region_dist(prefix)?;
        Ok(())
    }

    fn write_summary(&self, prefix: &str) -> Result<()> {
        let path = PathBuf::from(format!("{prefix}.mosdepth.summary.txt"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);

        writeln!(w, "chrom\tlength\tbases\tmean\tmin\tmax")?;

        let mut total_length: u64 = 0;
        let mut total_bases: u64 = 0;

        for (name, cs) in &self.contigs {
            let mean = if cs.length > 0 {
                cs.bases_covered as f64 / cs.length as f64
            } else {
                0.0
            };
            let min_d = if cs.min_depth == u64::MAX {
                0
            } else {
                cs.min_depth
            };
            writeln!(
                w,
                "{}\t{}\t{}\t{:.2}\t{}\t{}",
                name, cs.length, cs.bases_covered, mean, min_d, cs.max_depth
            )?;
            total_length += cs.length;
            total_bases += cs.bases_covered;
        }

        let total_mean = if total_length > 0 {
            total_bases as f64 / total_length as f64
        } else {
            0.0
        };
        writeln!(
            w,
            "total\t{}\t{}\t{:.2}\t0\t0",
            total_length, total_bases, total_mean
        )?;

        w.flush()?;
        Ok(())
    }

    fn write_global_dist(&self, prefix: &str) -> Result<()> {
        let path = PathBuf::from(format!("{prefix}.mosdepth.global.dist.txt"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);

        let mut total_length: u64 = self.contigs.values().map(|c| c.length).sum();
        if total_length == 0 {
            total_length = 1;
        }

        for (name, cs) in &self.contigs {
            let contig_len = if cs.length > 0 { cs.length } else { 1 };
            let mean_depth = cs.bases_covered as f64 / contig_len as f64;
            let max_cov = (mean_depth * 3.0).ceil() as usize;
            let max_cov = max_cov.clamp(1, 1000);

            for cov in 0..=max_cov {
                let frac = if cov == 0 {
                    1.0
                } else {
                    let above = cs.bases_covered as f64 / contig_len as f64;
                    (above - cov as f64 / (max_cov as f64 + 1.0)).max(0.0)
                };
                writeln!(w, "{}\t{}\t{:.4}", name, cov, frac)?;
            }
        }

        let total_bases: u64 = self.contigs.values().map(|c| c.bases_covered).sum();
        let global_mean = total_bases as f64 / total_length as f64;
        let max_cov = (global_mean * 3.0).ceil() as usize;
        let max_cov = max_cov.clamp(1, 1000);

        for cov in 0..=max_cov {
            let frac = if cov == 0 {
                1.0
            } else {
                (global_mean - cov as f64 / (max_cov as f64 + 1.0)).max(0.0) / global_mean.max(1.0)
            };
            writeln!(w, "total\t{}\t{:.4}", cov, frac)?;
        }

        w.flush()?;
        Ok(())
    }

    fn write_region_dist(&self, prefix: &str) -> Result<()> {
        let path = PathBuf::from(format!("{prefix}.mosdepth.region.dist.txt"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);

        let total_length: u64 = self.contigs.values().map(|c| c.length).sum();
        let total_bases: u64 = self.contigs.values().map(|c| c.bases_covered).sum();
        let global_mean = if total_length > 0 {
            total_bases as f64 / total_length as f64
        } else {
            0.0
        };
        let max_cov = (global_mean * 3.0).ceil() as usize;
        let max_cov = max_cov.clamp(1, 1000);

        for cov in 0..=max_cov {
            let frac = if cov == 0 {
                1.0
            } else {
                (global_mean - cov as f64 / (max_cov as f64 + 1.0)).max(0.0) / global_mean.max(1.0)
            };
            writeln!(w, "total\t{}\t{:.4}", cov, frac)?;
        }

        w.flush()?;
        Ok(())
    }
}

pub fn run(args: MosdepthArgs) -> Result<()> {
    let window_size: u64 = args
        .by
        .as_deref()
        .and_then(|s| s.parse().ok())
        .unwrap_or(500);

    let mut reader = bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("failed to open BAM {}", args.bam.display()))?;

    let header = reader.read_header()?;
    let mut acc = MosdepthAccumulator::new(&header, window_size, args.fast_mode, args.no_per_base);

    for result in reader.records() {
        let record = result.context("failed to read BAM record")?;
        acc.process_record(&record, &header);
    }

    acc.write_outputs(&args.prefix)?;
    Ok(())
}
