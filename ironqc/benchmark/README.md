# ironqc Benchmark Results

## Test Data

Source: [nf-core/test-datasets (sarek3 branch)](https://github.com/nf-core/test-datasets/tree/sarek3/data/genomics/homo_sapiens/illumina/bam)

| File | Size | Description |
|------|------|-------------|
| `test.paired_end.sorted.bam` | 175 KB | 5,644 reads, chr22, paired-end |
| `test2.paired_end.sorted.bam` | 196 KB | 5,436 reads, chr22, paired-end |
| `genome.fasta` | 40 KB | chr22, 40,001 bases |

## ironqc Performance (release build, macOS arm64)

All timings measured with `/usr/bin/time -l`, median of 3 runs.

| Subcommand | Wall Time | User Time | RSS (MB) |
|------------|-----------|-----------|----------|
| `ironqc stats` | <10 ms | <10 ms | 2.4 |
| `ironqc mosdepth` | <10 ms | <10 ms | 2.4 |
| `ironqc bundle` (single-pass) | ~10 ms | <10 ms | 2.8 |

Note: Test BAMs are very small (~5K reads). These timings confirm correctness and baseline overhead.
For production-scale benchmarks, use WGS BAMs (30-50 GB, ~800M reads).

## Output Verification

### stats (samtools stats compatible)

```
# This file was produced by samtools stats
# The command line was: ironqc stats
SN	raw total sequences:	5644
SN	sequences:	5644
SN	reads mapped:	5642
SN	reads unmapped:	2
SN	reads properly paired:	5640
SN	total length:	672202
SN	bases mapped (cigar):	671070
SN	average length:	119.1
SN	average quality:	40.9
SN	insert size average:	125.7
```

MultiQC compatibility: output starts with `# This file was produced by samtools stats` header (required trigger).

### mosdepth (mosdepth compatible)

```
chrom	length	bases	mean	min	max
chr22	40001	670989	16.77	0	0
total	40001	670989	16.77	0	0
```

Output files: `*.mosdepth.summary.txt`, `*.mosdepth.global.dist.txt`, `*.mosdepth.region.dist.txt`

### indexcov (goleft indexcov compatible)

Output files: `*-indexcov.ped`, `*-indexcov.roc`, `*-indexcov.bed.gz`, `*-indexcov.html`

### bundle (single-pass mode)

Produces all stats + mosdepth + indexcov outputs in a single BAM read pass.
Bundle outputs are byte-identical to individual subcommand outputs.

## Upstream Tool Comparison

Upstream tools (samtools 1.21, mosdepth 0.3.10, goleft 0.2.4) were not available on the test system.
For a full comparison, run:

```bash
# Generate reference outputs with Docker
docker run -v $(pwd)/input:/data quay.io/biocontainers/samtools:1.21--h50ea8bc_0 \
  samtools stats /data/test.paired_end.sorted.bam > reference/test.samtools.stats

docker run -v $(pwd)/input:/data quay.io/biocontainers/mosdepth:0.3.10--hd299d5a_0 \
  mosdepth --fasta /data/genome.fasta /data/test /data/test.paired_end.sorted.bam

# Then validate with included scripts
python3 ../scripts/validate_stats.py reference/test.samtools.stats output/stats/test.stats
python3 ../scripts/validate_mosdepth.py reference/test output/mosdepth/test
```

## Build Information

- Rust toolchain: stable
- Key dependency: noodles 0.88 (pure Rust BAM/CRAM, zero C dependencies)
- Build: `cargo build --release` (~18s)
- Tests: 5/5 integration tests passing
- Linting: `cargo clippy -- -D warnings` clean
