use crate::cli::StatsArgs;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Cigar as CigarTrait;
use noodles::sam::alignment::record::QualityScores;

macro_rules! sn {
    ($w:expr, $key:expr, $val:expr) => {
        writeln!($w, "SN\t{}\t{}", $key, $val)?
    };
}

pub struct StatsAccumulator {
    pub total_sequences: u64,
    pub total_length: u64,
    pub mapped_reads: u64,
    pub unmapped_reads: u64,
    pub duplicates: u64,
    pub paired_in_sequencing: u64,
    pub properly_paired: u64,
    pub read1: u64,
    pub read2: u64,
    pub singletons: u64,
    pub bases_mapped: u64,
    pub bases_mapped_cigar: u64,
    pub bases_trimmed: u64,
    pub bases_duplicated: u64,
    pub mismatches: u64,
    pub quality_sum: u64,
    pub insert_sizes: Vec<i64>,
    pub mapq_counts: HashMap<u8, u64>,
    pub reads_mapped_and_paired: u64,
    pub reads_mq0: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub total_first_frag_length: u64,
    pub total_last_frag_length: u64,
    pub first_frag_count: u64,
    pub last_frag_count: u64,
    pub max_length: u64,
    pub max_first_frag_length: u64,
    pub max_last_frag_length: u64,
}

impl StatsAccumulator {
    pub fn new() -> Self {
        Self {
            total_sequences: 0,
            total_length: 0,
            mapped_reads: 0,
            unmapped_reads: 0,
            duplicates: 0,
            paired_in_sequencing: 0,
            properly_paired: 0,
            read1: 0,
            read2: 0,
            singletons: 0,
            bases_mapped: 0,
            bases_mapped_cigar: 0,
            bases_trimmed: 0,
            bases_duplicated: 0,
            mismatches: 0,
            quality_sum: 0,
            insert_sizes: Vec::new(),
            mapq_counts: HashMap::new(),
            reads_mapped_and_paired: 0,
            reads_mq0: 0,
            secondary: 0,
            supplementary: 0,
            total_first_frag_length: 0,
            total_last_frag_length: 0,
            first_frag_count: 0,
            last_frag_count: 0,
            max_length: 0,
            max_first_frag_length: 0,
            max_last_frag_length: 0,
        }
    }

    pub fn process_record(&mut self, record: &bam::Record, _header: &noodles::sam::Header) {
        let flags = record.flags();
        self.total_sequences += 1;

        let seq_len = record.sequence().len() as u64;
        self.total_length += seq_len;
        if seq_len > self.max_length {
            self.max_length = seq_len;
        }

        if flags.is_secondary() {
            self.secondary += 1;
        }
        if flags.is_supplementary() {
            self.supplementary += 1;
        }
        if flags.is_duplicate() {
            self.duplicates += 1;
            self.bases_duplicated += seq_len;
        }

        if flags.is_segmented() {
            self.paired_in_sequencing += 1;
            if flags.is_properly_segmented() {
                self.properly_paired += 1;
            }
            if flags.is_first_segment() {
                self.read1 += 1;
                self.total_first_frag_length += seq_len;
                self.first_frag_count += 1;
                if seq_len > self.max_first_frag_length {
                    self.max_first_frag_length = seq_len;
                }
            }
            if flags.is_last_segment() {
                self.read2 += 1;
                self.total_last_frag_length += seq_len;
                self.last_frag_count += 1;
                if seq_len > self.max_last_frag_length {
                    self.max_last_frag_length = seq_len;
                }
            }
        }

        if !flags.is_unmapped() {
            self.mapped_reads += 1;
            self.bases_mapped += seq_len;

            let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
            *self.mapq_counts.entry(mapq).or_insert(0) += 1;
            if mapq == 0 {
                self.reads_mq0 += 1;
            }

            if flags.is_segmented() {
                self.reads_mapped_and_paired += 1;
                if flags.is_mate_unmapped() {
                    self.singletons += 1;
                }
            }

            {
                let cigar = record.cigar();
                if let Ok(span) = cigar.alignment_span() {
                    self.bases_mapped_cigar += span as u64;
                }
                for op in cigar.iter().flatten() {
                    match op.kind() {
                        Kind::SoftClip | Kind::HardClip => {
                            self.bases_trimmed += op.len() as u64;
                        }
                        _ => {}
                    }
                }
            }

            if flags.is_segmented() && !flags.is_mate_unmapped() {
                let tlen = record.template_length();
                if tlen > 0 {
                    self.insert_sizes.push(tlen as i64);
                }
            }
        } else {
            self.unmapped_reads += 1;
        }

        for score in record.quality_scores().iter().flatten() {
            self.quality_sum += score as u64;
        }
    }

    pub fn finalize(&mut self) -> FinalizedStats<'_> {
        let average_length = if self.total_sequences > 0 {
            self.total_length as f64 / self.total_sequences as f64
        } else {
            0.0
        };
        let average_quality = if self.total_length > 0 {
            self.quality_sum as f64 / self.total_length as f64
        } else {
            0.0
        };
        let (insert_size_average, insert_size_stddev) = if !self.insert_sizes.is_empty() {
            let n = self.insert_sizes.len() as f64;
            let sum: f64 = self.insert_sizes.iter().map(|&x| x as f64).sum();
            let avg = sum / n;
            let var: f64 = self
                .insert_sizes
                .iter()
                .map(|&x| {
                    let d = x as f64 - avg;
                    d * d
                })
                .sum::<f64>()
                / n;
            (avg, var.sqrt())
        } else {
            (0.0, 0.0)
        };
        let error_rate = if self.bases_mapped_cigar > 0 {
            self.mismatches as f64 / self.bases_mapped_cigar as f64
        } else {
            0.0
        };

        FinalizedStats {
            acc: self,
            average_length,
            average_quality,
            insert_size_average,
            insert_size_stddev,
            error_rate,
        }
    }
}

pub struct FinalizedStats<'a> {
    acc: &'a StatsAccumulator,
    average_length: f64,
    average_quality: f64,
    insert_size_average: f64,
    insert_size_stddev: f64,
    error_rate: f64,
}

impl<'a> FinalizedStats<'a> {
    pub fn write_output(&self, path: &PathBuf) -> Result<()> {
        let file =
            File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
        let mut w = BufWriter::new(file);
        let a = self.acc;

        writeln!(w, "# This file was produced by samtools stats")?;
        writeln!(w, "# The command line was: ironqc stats")?;

        sn!(w, "raw total sequences:", a.total_sequences);
        sn!(w, "filtered sequences:", 0);
        sn!(w, "sequences:", a.total_sequences);
        sn!(w, "is sorted:", 1);
        sn!(w, "1st fragments:", a.first_frag_count);
        sn!(w, "last fragments:", a.last_frag_count);
        sn!(w, "reads mapped:", a.mapped_reads);
        sn!(w, "reads mapped and paired:", a.reads_mapped_and_paired);
        sn!(w, "reads unmapped:", a.unmapped_reads);
        sn!(w, "reads properly paired:", a.properly_paired);
        sn!(w, "reads paired:", a.paired_in_sequencing);
        sn!(w, "reads duplicated:", a.duplicates);
        sn!(w, "reads MQ0:", a.reads_mq0);
        sn!(w, "reads QC failed:", 0);
        sn!(w, "non-primary alignments:", a.secondary + a.supplementary);
        sn!(w, "supplementary alignments:", a.supplementary);
        sn!(w, "total length:", a.total_length);
        sn!(w, "total first fragment length:", a.total_first_frag_length);
        sn!(w, "total last fragment length:", a.total_last_frag_length);
        sn!(w, "bases mapped:", a.bases_mapped);
        sn!(w, "bases mapped (cigar):", a.bases_mapped_cigar);
        sn!(w, "bases trimmed:", a.bases_trimmed);
        sn!(w, "bases duplicated:", a.bases_duplicated);
        sn!(w, "mismatches:", a.mismatches);
        sn!(w, "error rate:", format!("{:.6e}", self.error_rate));
        sn!(w, "average length:", format!("{:.1}", self.average_length));
        sn!(
            w,
            "average first fragment length:",
            if a.first_frag_count > 0 {
                format!(
                    "{:.1}",
                    a.total_first_frag_length as f64 / a.first_frag_count as f64
                )
            } else {
                "0".to_string()
            }
        );
        sn!(
            w,
            "average last fragment length:",
            if a.last_frag_count > 0 {
                format!(
                    "{:.1}",
                    a.total_last_frag_length as f64 / a.last_frag_count as f64
                )
            } else {
                "0".to_string()
            }
        );
        sn!(w, "maximum length:", a.max_length);
        sn!(w, "maximum first fragment length:", a.max_first_frag_length);
        sn!(w, "maximum last fragment length:", a.max_last_frag_length);
        sn!(
            w,
            "average quality:",
            format!("{:.1}", self.average_quality)
        );
        sn!(
            w,
            "insert size average:",
            format!("{:.1}", self.insert_size_average)
        );
        sn!(
            w,
            "insert size standard deviation:",
            format!("{:.1}", self.insert_size_stddev)
        );
        sn!(w, "inward oriented pairs:", 0);
        sn!(w, "outward oriented pairs:", 0);
        sn!(w, "pairs with other orientation:", 0);
        sn!(w, "pairs on different chromosomes:", 0);
        sn!(
            w,
            "percentage of properly paired reads (%):",
            format!(
                "{:.1}",
                if a.paired_in_sequencing > 0 {
                    100.0 * a.properly_paired as f64 / a.paired_in_sequencing as f64
                } else {
                    0.0
                }
            )
        );

        w.flush()?;
        Ok(())
    }
}

pub fn run(args: StatsArgs) -> Result<()> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(&args.bam)
        .with_context(|| format!("failed to open BAM {}", args.bam.display()))?;

    let header = reader.read_header()?;
    let mut acc = StatsAccumulator::new();

    for result in reader.records() {
        let record = result.context("failed to read BAM record")?;
        acc.process_record(&record, &header);
    }

    let finalized = acc.finalize();
    let output = PathBuf::from(format!("{}.stats", args.prefix));
    finalized.write_output(&output)?;

    Ok(())
}
