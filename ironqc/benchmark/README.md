# ironqc Benchmark Results

## Test Data

### nf-core/test-datasets raredisease (primary benchmark)

Source: [nf-core/test-datasets (raredisease branch)](https://github.com/nf-core/test-datasets/tree/raredisease/testdata)

| File | Size | Reads | Contigs |
|------|------|-------|---------|
| `earlycasualcaiman_sorted_md.bam` | 69 MB | 1,351,181 | chr21 + MT + 22 stub contigs |
| `hugelymodelbat_sorted_md.bam` | 56 MB | 1,048,330 | chr21 + MT + 22 stub contigs |
| `slowlycivilbuck_sorted_md.bam` | 52 MB | 992,459 | chr21 + MT + 22 stub contigs |
| `reference.fasta` | 47 MB | — | 24 contigs (chr21: 48M bases, MT: 16.5K bases) |

### nf-core/test-datasets sarek3 (micro benchmark)

Source: [nf-core/test-datasets (sarek3 branch)](https://github.com/nf-core/test-datasets/tree/sarek3/data/genomics/homo_sapiens/illumina/bam)

| File | Size | Reads |
|------|------|-------|
| `test.paired_end.sorted.bam` | 175 KB | 5,644 |
| `test2.paired_end.sorted.bam` | 196 KB | 5,436 |

## Performance Results (release build, macOS arm64 — Apple Silicon)

All timings measured with `/usr/bin/time -l`.

### Raredisease BAMs (~1M reads each)

| Sample | Subcommand | Wall Time | User Time | RSS (MB) |
|--------|------------|-----------|-----------|----------|
| earlycasualcaiman (69 MB, 1.35M reads) | `stats` | 1.02s | 0.98s | 10.6 |
| earlycasualcaiman | `mosdepth` | 0.93s | 0.89s | 2.6 |
| earlycasualcaiman | `bundle` | 1.06s | 1.02s | 9.9 |
| hugelymodelbat (56 MB, 1.05M reads) | `stats` | 0.80s | 0.76s | 8.5 |
| hugelymodelbat | `mosdepth` | 0.74s | 0.71s | 2.6 |
| hugelymodelbat | `bundle` | 0.84s | 0.81s | 8.8 |
| slowlycivilbuck (52 MB, 992K reads) | `stats` | 0.75s | 0.72s | 8.3 |
| slowlycivilbuck | `mosdepth` | 0.69s | 0.66s | 2.6 |
| slowlycivilbuck | `bundle` | 0.79s | 0.76s | 8.8 |

**Key observations:**
- **Bundle is faster than stats + mosdepth separately** (single BAM read pass vs two passes)
- **Throughput**: ~1.3M reads/sec (stats), ~1.5M reads/sec (mosdepth), ~1.3M reads/sec (bundle)
- **Memory**: stats ~9 MB (quality/cycle histograms), mosdepth ~2.6 MB (coverage arrays only), bundle ~9 MB (dominated by stats histograms)

### Sarek3 BAMs (~5K reads, baseline overhead)

| Subcommand | Wall Time | RSS (MB) |
|------------|-----------|----------|
| `stats` | <10 ms | 2.4 |
| `mosdepth` | <10 ms | 2.4 |
| `bundle` | ~10 ms | 2.8 |

## Output Verification

### stats (samtools stats compatible)

```
# This file was produced by samtools stats
SN	raw total sequences:	1351181
SN	reads mapped:	1347040
SN	reads properly paired:	1342930
SN	total length:	202104875
SN	bases mapped (cigar):	198254694
SN	average length:	149.6
SN	average quality:	35.5
SN	insert size average:	51802.9
```

MultiQC compatibility: output starts with `# This file was produced by samtools stats` (required trigger).

### mosdepth (mosdepth compatible)

```
chrom	length	bases	mean	min	max
21	48129895	52191506	1.08	0	0
MT	16569	146179119	8822.45	0	0
total	48168464	198370625	4.12	0	0
```

Output files: `*.mosdepth.summary.txt`, `*.mosdepth.global.dist.txt`, `*.mosdepth.region.dist.txt`

### indexcov (goleft indexcov compatible)

Output files: `*-indexcov.ped`, `*-indexcov.roc`, `*-indexcov.bed.gz`, `*-indexcov.html`

### bundle (single-pass mode)

Produces all stats + mosdepth + indexcov outputs in one BAM read pass.
Bundle outputs are byte-identical to individual subcommand outputs (verified on sarek3 test data).

## Upstream Tool Comparison

Upstream tools (samtools 1.21, mosdepth 0.3.10, goleft 0.2.4) were not available on the test system.
To generate reference outputs for numeric comparison:

```bash
docker run -v $(pwd)/input/raredisease:/data quay.io/biocontainers/samtools:1.21--h50ea8bc_0 \
  samtools stats -r /data/reference.fasta /data/hugelymodelbat_sorted_md.bam \
  > reference/hugelymodelbat.samtools.stats

docker run -v $(pwd)/input/raredisease:/data quay.io/biocontainers/mosdepth:0.3.10--hd299d5a_0 \
  mosdepth --fasta /data/reference.fasta /data/hugelymodelbat \
  /data/hugelymodelbat_sorted_md.bam

python3 ../scripts/validate_stats.py \
  reference/hugelymodelbat.samtools.stats \
  output/raredisease/stats/hugelymodelbat_sorted_md.stats

python3 ../scripts/validate_mosdepth.py \
  reference/hugelymodelbat \
  output/raredisease/mosdepth/hugelymodelbat_sorted_md
```

## Build Information

- **Rust toolchain**: stable
- **Key dependency**: noodles 0.88 (pure Rust BAM/CRAM, zero C dependencies)
- **Build time**: `cargo build --release` ~18s
- **Tests**: 5/5 integration tests passing
- **Linting**: `cargo clippy -- -D warnings` clean
