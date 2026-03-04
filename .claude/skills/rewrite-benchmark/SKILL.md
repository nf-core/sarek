---
name: rewrite-benchmark
description: Performance benchmarking and comparison against upstream tools. Use after validation passes to measure speedup, memory usage, and thread scaling. Produces reproducible benchmark scripts and structured results documentation.
model: sonnet
user-invocable: true
argument-hint: "[--small <small-input>] [--large <large-input>]"
---

# /rewrite-benchmark — Performance Benchmarking

You are benchmarking a rewritten bioinformatics tool against the upstream tool(s) it
replaces. Your job is to produce rigorous, reproducible performance comparisons and
document them clearly.

**Arguments:** `$ARGUMENTS` may specify input files. If not provided, look for test
data paths in `TODO-rewrite-plan.md`.

## Pre-Benchmark Checklist

1. **Validation passes** — Run `/rewrite-validate` first. Never benchmark incorrect code.
2. **Release build** — Benchmark the release binary, not debug:
   ```bash
   cargo build --release
   ```
3. **Upstream tool available** — Ensure the upstream tool is installed and runnable.
4. **Reference outputs from containers** — Reference outputs generated from Docker containers with pinned version tags — never use locally-installed tools for reference generation.
5. **Test data ready** — Both small (quick) and large (realistic) inputs.
6. **System is idle** — Close unnecessary applications. Note background load.

## Test Conditions

Document ALL of these in the benchmark report:

```markdown
## Test Conditions
- **Hardware**: {CPU model}, {cores}C/{threads}T, {RAM}
- **OS**: {distro} {version}, kernel {version}
- **Storage**: {SSD/HDD/NFS}, {filesystem}
- **Rewritten tool**: {name} v{version}, commit {hash}
- **Upstream tool(s)**: {name} v{version}
- **Input data (small)**: {description}, {size}, {read count if BAM}
- **Input data (large)**: {description}, {size}, {read count if BAM}
```

## Benchmarking Method

### Wall Clock Time

Use `hyperfine` for robust timing (handles warmup, multiple runs, statistical analysis):

```bash
# Install if needed
cargo install hyperfine

# Benchmark with warmup and multiple runs
hyperfine \
  --warmup 1 \
  --min-runs 3 \
  --export-json benchmark/results.json \
  --export-markdown benchmark/results.md \
  'upstream-tool input.bam annotation.gtf' \
  './target/release/rewritten-tool input.bam --gtf annotation.gtf'
```

If `hyperfine` is unavailable, use `/usr/bin/time -v` (Linux) or `gtime -v` (macOS):

```bash
/usr/bin/time -v ./target/release/tool args 2> benchmark/timing.txt
```

### Peak Memory

Extract from `/usr/bin/time -v` output:
- `Maximum resident set size (kbytes)` — peak RSS

Or use `/usr/bin/time -l` on macOS:
- `maximum resident set size` — in bytes

#### macOS vs Linux Notes

- **Linux**: `/usr/bin/time -v` reports kbytes. Parse with `grep "Maximum resident"`.
- **macOS**: `/usr/bin/time -l` (not `-v`) reports bytes. Parse with `grep "maximum resident"` (lowercase). Divide by 1024 to get kbytes for consistent reporting.
- If `gtime` (GNU time from Homebrew) is available on macOS, prefer it for Linux-compatible output: `gtime -v`.

### CPU Utilization

From `/usr/bin/time -v`:
- `Percent of CPU this job got` — should be close to N*100% for N threads

### Thread Scaling

Test with different thread counts to measure parallelism efficiency:

```bash
for threads in 1 2 4 8; do
  hyperfine \
    --warmup 1 \
    --min-runs 3 \
    "./target/release/tool --threads $threads input.bam --gtf annotation.gtf" \
    --export-json "benchmark/threads_${threads}.json"
done
```

## Comparison Rules

- **Define per-tool comparison rules upfront.** Different output types need different comparison strategies — exact match for tab-delimited counts, numeric tolerance (e.g., ±0.01) for floating-point statistics, line filtering for outputs with non-deterministic headers or timestamps. Never apply a single comparison strategy across all output types.
- **Separate correctness validation from performance benchmarking.** These are independent concerns. Automate correctness checks first (via `/rewrite-validate`); only benchmark after you trust the outputs are correct.
- **Reference data must match your tool's exact granularity and configuration.** Transcript-level vs gene-level quantification, the same annotation files, the same flags. A granularity or configuration mismatch makes the comparison meaningless — you'll chase "errors" that are actually expected differences.

## Benchmark Scripts

Write reproducible scripts in `benchmark/`:

```bash
#!/bin/bash
# benchmark/run_benchmark.sh — Reproducible benchmark script
set -euo pipefail

INPUT_SMALL="benchmark/input/small/test.bam"
INPUT_LARGE="benchmark/input/large/sample.bam"
ANNOTATION="benchmark/input/annotation.gtf"
REWRITTEN="./target/release/tool"

echo "=== Small input ==="
hyperfine --warmup 1 --min-runs 5 \
  "$REWRITTEN $INPUT_SMALL --gtf $ANNOTATION --outdir /tmp/bench_rewritten" \
  "upstream-tool $INPUT_SMALL $ANNOTATION /tmp/bench_upstream"

echo "=== Large input ==="
hyperfine --warmup 1 --min-runs 3 \
  "$REWRITTEN $INPUT_LARGE --gtf $ANNOTATION --outdir /tmp/bench_rewritten" \
  "upstream-tool $INPUT_LARGE $ANNOTATION /tmp/bench_upstream"
```

## Documenting Results

### benchmark/README.md

Write results using `benchmark-template.md` in this skill directory.

### Top-level README.md

Add a "Performance" section to the project README with headline numbers:

```markdown
## Performance

| Input | Upstream | Rewritten | Speedup |
|-------|----------|-----------|---------|
| Small (50K reads) | 12.3s | 0.8s | 15x |
| Large (100M reads) | 3h 47m | 4m 12s | 54x |

Benchmarked on {CPU}, {RAM}. See `benchmark/README.md` for methodology.
```

### Cross-Reference Check

**CRITICAL**: Ensure numbers in the top-level README match `benchmark/README.md` exactly.
These documents must always be updated together. When re-running benchmarks, update BOTH.

### Additional Guidelines

- **Generate benchmark visualizations from machine-measured data** (e.g., `results.json` from hyperfine). Never hand-craft benchmark charts or tables with manually typed numbers — they fall out of sync immediately.
- **Keep a SINGLE source of truth for benchmark numbers.** If numbers appear in both `benchmark/README.md` and the top-level `README.md`, they MUST match. Always update both together as the very last step after benchmarking is complete.
- **Mark unmeasured tools as "not measured"** rather than using approximate (`~`) estimates. If a tool hasn't been profiled, say so explicitly — guesses erode trust in the entire benchmark.

## What Constitutes 'Good Enough'

- **I/O-bound tools (BAM processing):** 2–5x speedup is realistic from single-pass architecture alone. 10x+ requires eliminating redundant file passes (e.g., bundling multiple tools into one pass).
- **CPU-bound tools (statistical computation):** 10–100x speedup over Python/R is typical.
- **The primary value of a Rust rewrite often isn't raw speed** — it's eliminating redundant BAM passes, reducing dependencies (single static binary vs. a Python/R ecosystem), and combining multiple tools into one.
- **Always measure total pipeline impact, not just individual tool speedup.** Eliminating 4 task submissions (container pulls, file staging, scheduler overhead) can matter more than 2x per-task speedup. Report both isolated tool speedup and estimated pipeline-level improvement.

## Output

- Write `TODO-rewrite-benchmark.md` with results and methodology
- Write/update `benchmark/README.md` using the template
- Update top-level `README.md` with headline performance numbers
- Save benchmark scripts in `benchmark/`
- Save raw benchmark data (JSON from hyperfine) in `benchmark/`
