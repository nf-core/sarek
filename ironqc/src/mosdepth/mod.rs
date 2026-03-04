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
    summary_deltas: Option<Vec<i32>>,
    dist_deltas: Option<Vec<i32>>,
}

struct ComputedContigStats {
    length: u64,
    bases_covered: u64,
    min_depth: u32,
    max_depth: u32,
    global_histogram: Vec<u64>,
    region_histogram: Vec<u64>,
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
                    summary_deltas: None,
                    dist_deltas: None,
                },
            );
        }
        Self {
            contigs,
            window_size: window_size.max(1),
            fast_mode,
        }
    }

    pub fn process_record(&mut self, record: &bam::Record, header: &noodles::sam::Header) {
        let flags = record.flags();
        if flags.is_unmapped()
            || flags.is_secondary()
            || flags.is_supplementary()
            || flags.is_duplicate()
            || flags.is_qc_fail()
        {
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

        let len = match usize::try_from(contig.length) {
            Ok(v) => v,
            Err(_) => return,
        };

        let summary_deltas = contig
            .summary_deltas
            .get_or_insert_with(|| vec![0; len.saturating_add(1)]);
        let dist_deltas = contig
            .dist_deltas
            .get_or_insert_with(|| vec![0; len.saturating_add(1)]);

        let dist_start = pos.min(contig.length);
        let dist_end = pos
            .saturating_add(record.sequence().len() as u64)
            .min(contig.length);
        Self::add_delta_interval(dist_deltas, dist_start, dist_end);

        if self.fast_mode {
            Self::add_delta_interval(summary_deltas, dist_start, dist_end);
            return;
        }

        let cigar = record.cigar();
        let mut ref_pos = pos;

        for op in cigar.iter().flatten() {
            let op_len = op.len() as u64;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                    let start = ref_pos.min(contig.length);
                    let end = ref_pos.saturating_add(op_len).min(contig.length);
                    Self::add_delta_interval(summary_deltas, start, end);

                    ref_pos = ref_pos.saturating_add(op_len);
                }
                Kind::Deletion | Kind::Skip => {
                    let start = ref_pos.min(contig.length);
                    let end = ref_pos.saturating_add(op_len).min(contig.length);
                    Self::add_delta_interval(summary_deltas, start, end);

                    ref_pos = ref_pos.saturating_add(op_len);
                }
                Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
            }
        }
    }

    pub fn write_outputs(&self, prefix: &str) -> Result<()> {
        let computed = self.computed_contigs();
        self.write_summary(prefix, &computed)?;
        self.write_global_dist(prefix, &computed)?;
        self.write_region_dist(prefix, &computed)?;
        Ok(())
    }

    fn write_summary(&self, prefix: &str, computed: &[(&str, ComputedContigStats)]) -> Result<()> {
        let path = PathBuf::from(format!("{prefix}.mosdepth.summary.txt"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);

        writeln!(w, "chrom\tlength\tbases\tmean\tmin\tmax")?;

        let mut total_length: u64 = 0;
        let mut total_bases: u64 = 0;
        let mut total_min: Option<u32> = None;
        let mut total_max: u32 = 0;

        for (name, cs) in computed {
            let mean = if cs.length > 0 {
                cs.bases_covered as f64 / cs.length as f64
            } else {
                0.0
            };

            writeln!(
                w,
                "{}\t{}\t{}\t{:.2}\t{}\t{}",
                name, cs.length, cs.bases_covered, mean, cs.min_depth, cs.max_depth
            )?;

            writeln!(
                w,
                "{}_region\t{}\t{}\t{:.2}\t{}\t{}",
                name, cs.length, cs.bases_covered, mean, cs.min_depth, cs.max_depth
            )?;

            total_length += cs.length;
            total_bases += cs.bases_covered;
            total_min = Some(match total_min {
                Some(current) => current.min(cs.min_depth),
                None => cs.min_depth,
            });
            total_max = total_max.max(cs.max_depth);
        }

        let total_mean = if total_length > 0 {
            total_bases as f64 / total_length as f64
        } else {
            0.0
        };

        let total_min = total_min.unwrap_or(0);
        writeln!(
            w,
            "total\t{}\t{}\t{:.2}\t{}\t{}",
            total_length, total_bases, total_mean, total_min, total_max
        )?;

        writeln!(
            w,
            "total_region\t{}\t{}\t{:.2}\t{}\t{}",
            total_length, total_bases, total_mean, total_min, total_max
        )?;

        w.flush()?;
        Ok(())
    }

    fn write_global_dist(
        &self,
        prefix: &str,
        computed: &[(&str, ComputedContigStats)],
    ) -> Result<()> {
        let path = PathBuf::from(format!("{prefix}.mosdepth.global.dist.txt"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);

        let mut total_length: u64 = 0;
        let mut total_histogram: Vec<u64> = Vec::new();

        for (name, cs) in computed {
            total_length += cs.length;
            Self::merge_histogram(&mut total_histogram, &cs.global_histogram);
            Self::write_distribution_rows(&mut w, name, cs.length, &cs.global_histogram)?;
        }

        if total_length > 0 {
            Self::write_distribution_rows(&mut w, "total", total_length, &total_histogram)?;
        }

        w.flush()?;
        Ok(())
    }

    fn write_region_dist(
        &self,
        prefix: &str,
        computed: &[(&str, ComputedContigStats)],
    ) -> Result<()> {
        let path = PathBuf::from(format!("{prefix}.mosdepth.region.dist.txt"));
        let file =
            File::create(&path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);

        let mut total_length: u64 = 0;
        let mut total_histogram: Vec<u64> = Vec::new();

        for (name, cs) in computed {
            total_length += cs.length;
            Self::merge_histogram(&mut total_histogram, &cs.region_histogram);
            Self::write_distribution_rows(&mut w, name, cs.length, &cs.region_histogram)?;
        }

        if total_length > 0 {
            Self::write_distribution_rows(&mut w, "total", total_length, &total_histogram)?;
        }

        w.flush()?;
        Ok(())
    }

    fn computed_contigs(&self) -> Vec<(&str, ComputedContigStats)> {
        let mut computed = Vec::new();

        for (name, contig) in &self.contigs {
            let Some(summary_deltas) = &contig.summary_deltas else {
                continue;
            };

            let Some(dist_deltas) = &contig.dist_deltas else {
                continue;
            };

            let length_usize = match usize::try_from(contig.length) {
                Ok(v) => v,
                Err(_) => continue,
            };

            if length_usize == 0 {
                continue;
            }

            let mut summary_depth: i64 = 0;
            let mut bases_covered: u64 = 0;
            let mut min_depth: u32 = u32::MAX;
            let mut max_depth: u32 = 0;

            for delta in summary_deltas.iter().take(length_usize) {
                summary_depth += i64::from(*delta);
                let d = summary_depth.max(0) as u32;

                bases_covered += u64::from(d);
                min_depth = min_depth.min(d);
                max_depth = max_depth.max(d);
            }

            let mut global_histogram: Vec<u64> = vec![0];
            let mut region_histogram: Vec<u64> = vec![0];
            let mut dist_depth: i64 = 0;

            let mut window_sum: u64 = 0;
            let mut window_bases: u64 = 0;

            for delta in dist_deltas.iter().take(length_usize) {
                dist_depth += i64::from(*delta);
                let d = dist_depth.max(0) as u32;

                let global_depth = d as usize;
                if global_depth >= global_histogram.len() {
                    global_histogram.resize(global_depth + 1, 0);
                }
                global_histogram[global_depth] += 1;

                window_sum += u64::from(d);
                window_bases += 1;

                if window_bases == self.window_size {
                    let region_depth =
                        (window_sum as f64 / self.window_size as f64).floor() as usize;
                    if region_depth >= region_histogram.len() {
                        region_histogram.resize(region_depth + 1, 0);
                    }
                    region_histogram[region_depth] += window_bases;
                    window_sum = 0;
                    window_bases = 0;
                }
            }

            if window_bases > 0 {
                let region_depth = (window_sum as f64 / self.window_size as f64).floor() as usize;
                if region_depth >= region_histogram.len() {
                    region_histogram.resize(region_depth + 1, 0);
                }
                region_histogram[region_depth] += window_bases;
            }

            computed.push((
                name.as_str(),
                ComputedContigStats {
                    length: contig.length,
                    bases_covered,
                    min_depth,
                    max_depth,
                    global_histogram,
                    region_histogram,
                },
            ));
        }

        computed
    }

    fn add_delta_interval(deltas: &mut [i32], start: u64, end: u64) {
        if start >= end {
            return;
        }

        let start_idx = start as usize;
        let end_idx = end as usize;
        deltas[start_idx] += 1;
        deltas[end_idx] -= 1;
    }

    fn merge_histogram(total: &mut Vec<u64>, histogram: &[u64]) {
        if histogram.len() > total.len() {
            total.resize(histogram.len(), 0);
        }

        for (idx, value) in histogram.iter().enumerate() {
            total[idx] += value;
        }
    }

    fn write_distribution_rows(
        w: &mut BufWriter<File>,
        label: &str,
        length: u64,
        histogram: &[u64],
    ) -> Result<()> {
        if length == 0 {
            return Ok(());
        }

        let mut max_depth = histogram.len().saturating_sub(1);
        while max_depth > 0 && histogram[max_depth] == 0 {
            max_depth -= 1;
        }

        let mut cumulative: u64 = 0;
        for depth in (0..=max_depth).rev() {
            cumulative += histogram[depth];
            let frac = cumulative as f64 / length as f64;
            writeln!(w, "{label}\t{depth}\t{frac:.2}")?;
        }

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
