---
name: rewrite-optimize
description: Iteratively optimize rewritten tool for maximum performance without output regressions. Use after initial benchmarking when further speedup is desired. Runs a profile-optimize-validate-benchmark loop with common bioinformatics Rust optimization patterns (single-pass BAM, rayon parallelism, allocation reduction, I/O buffering).
model: sonnet
user-invocable: true
argument-hint: "[--target <module>]"
---

# Performance Optimization Loop

Iteratively optimize the rewritten tool for speed, memory efficiency, and thread scaling.
$ARGUMENTS may specify a target module to focus on.

**CARDINAL RULE: No output regressions. Run full validation after every change.**

## How This Works

Optimization is an **iterative loop**, not a single pass. You will:

1. Profile to find the current bottleneck
2. Apply ONE targeted optimization
3. Validate outputs (zero regressions)
4. Benchmark to confirm improvement
5. Commit if improved, revert if not
6. Repeat from step 1

**Stop when:** The top hotspot in the flamegraph accounts for < 5% of total time AND
no single optimization opportunity would yield > 5% speedup. Document this stopping
decision explicitly in the optimization log.

Also compare against the **upstream tool** timing. Even if self-optimization has
plateaued, the overall speedup vs upstream is what matters for the rewrite justification.

## Pre-Optimization Baseline (Do This Once)

1. Run full validation first (`/rewrite-validate`). Every test must pass.
2. Record current timing on both small and large test data:
   ```bash
   hyperfine --warmup 1 --min-runs 3 './target/release/<binary> <args>'
   ```
3. Record peak memory:
   ```bash
   # Linux:
   /usr/bin/time -v ./target/release/<binary> <args> 2>&1 | grep "Maximum resident"
   # macOS:
   /usr/bin/time -l ./target/release/<binary> <args> 2>&1 | grep "maximum resident"
   # Note: macOS reports bytes, Linux reports kbytes
   ```
4. Generate initial flamegraph:
   ```bash
   cargo install flamegraph  # if not already installed
   cargo flamegraph --release -- <args>
   ```
5. Save baseline numbers to the optimization log (see Output section).

## The Optimization Loop

### Step 1: Profile — Measure, Don't Guess

**Do not optimize based on intuition. Profile first every iteration.**

```bash
# CPU flamegraph (Linux)
cargo flamegraph --release -- <args>

# macOS
cargo instruments --release -t "CPU Profiler" -- <args>

# Memory profiling (when memory is the concern)
valgrind --tool=dhat ./target/release/<binary> <args>
heaptrack ./target/release/<binary> <args>
```

> **macOS note:** `cargo flamegraph` requires `dtrace` permissions and may fail.
> Alternatives:
> - `cargo instruments --release -t "CPU Profiler"` (requires Xcode Instruments)
> - `samply record ./target/release/<binary> <args>` (samply crate, works without root)
>
> For memory profiling on macOS, `valgrind` and `heaptrack` are unavailable. Use:
> - `cargo instruments -t Allocations`
> - Manual peak-memory tracking with `/usr/bin/time -l`

Identify the SINGLE biggest hotspot:
- Functions consuming > 5% of total time
- Allocation-heavy paths (many small allocations in tight loops)
- I/O wait time vs compute time
- Lock contention in parallel code

### Step 2: Apply ONE Targeted Optimization

Pick the highest-impact optimization for the identified bottleneck. Make ONE change
at a time so you can measure its effect cleanly and revert if it doesn't help.

Consult the optimization patterns below for common wins.

### Step 3: Validate

```bash
cargo fmt && cargo clippy -- -D warnings && cargo test --release
```

Then run full output validation against reference. **Zero regressions allowed.**
If validation fails, fix the regression or revert the optimization entirely.

### Step 4: Benchmark

```bash
hyperfine --warmup 1 --min-runs 3 './target/release/<binary> <args>'
```

Compare against the previous iteration's timing (not just the original baseline).

### Step 5: Commit or Revert

- If timing improved AND validation passes: `git commit` with descriptive message
- If timing didn't improve: `git checkout -- .` to revert
- If timing improved but validation fails: revert and investigate

Log the result in the optimization log (see Output section).

### Step 6: Decide — Continue or Stop

Re-profile (step 1). Look at the new flamegraph. Ask:
- Is there still a hotspot consuming > 5% of total time?
- Is there a concrete optimization that would yield > 5% speedup?
- Have we reached diminishing returns?

**If yes to the first two:** continue the loop.
**If no:** stop and document the stopping decision.

## Common Optimization Patterns for Bioinformatics Rust

Apply these in roughly this order — earlier items tend to have larger impact.

### 1. Single-Pass Architecture (Biggest Win)
If processing BAM files with multiple analyses, iterate the BAM **once** and dispatch
each record to all analyses in the same pass. This eliminates redundant I/O and
BAM decompression. For multi-tool rewrites, this alone can give 5-10x speedup.

### 2. Parallelism with Rayon
```rust
use rayon::prelude::*;

// Per-chromosome parallelism
chromosomes.par_iter().for_each(|chrom| {
    process_chromosome(chrom);
});

// Per-record parallelism (careful with ordering requirements)
records.par_chunks(1000).for_each(|chunk| {
    process_batch(chunk);
});
```
- Test thread scaling: 1, 2, 4, 8 threads
- Use `--threads` CLI flag, wire to rayon's thread pool
- Some analyses require ordered output — collect results and sort after parallel phase

### 3. BAM I/O Optimization
- `rust-htslib::bam::Reader::set_threads()` for parallel decompression
- Avoid repeated header lookups — cache chromosome name-to-tid mapping
- Filter early: check flags before extracting expensive fields (sequence, quality)
- Use `record.pos()` (i64) directly, avoid converting to other types

### 4. Reduce Allocations in Hot Loops
```rust
// BAD: allocates String per record
for record in bam_reader {
    let name = String::from_utf8_lossy(record.qname()).to_string();
}

// GOOD: work with &[u8] directly
for record in bam_reader {
    let name = record.qname(); // &[u8], no allocation
}
```

```rust
// BAD: clone() when you're done with the original
let data = source.clone();

// GOOD: take() for zero-cost ownership transfer
let data = std::mem::take(&mut source);

// BAD: String formatting for hash keys in hot loops
let key = format!("{}:{}", chrom, pos);

// GOOD: hash-based keys (FNV-1a or similar)
let key = hash_position(tid, pos);

// BAD: per-read chromosome name lookup via string
let chrom = header.tid2name(record.tid());

// GOOD: pre-compute per-TID lookup table
let lookup: Vec<_> = (0..n_refs).map(|tid| resolve_name(tid)).collect();
let result = lookup[record.tid() as usize];
```

- Pre-allocate `HashMap` with `HashMap::with_capacity(expected_size)`
- Reuse buffers across loop iterations
- Use `&str` / `&[u8]` instead of `String` / `Vec<u8>` where possible
- Avoid `format!()` in hot paths

> In per-read BAM processing, every allocation in the hot loop multiplies by millions.
> A 10GB BAM has ~200M reads. Even a single `Vec<u8>` allocation per read creates
> 200M+ heap allocations.

### 5. Efficient Data Structures
- `HashMap` for unordered lookups (fastest)
- `IndexMap` ONLY when insertion order matters (gene ordering for output)
- `coitrees::COITree` for interval queries (cache-oblivious, very fast)
- `Vec` over `LinkedList` always
- Consider `ahash::AHashMap` for faster hashing (no cryptographic guarantees needed)

### 6. General I/O
- Use `BufWriter` for all output (default buffer or larger)
- Write output files in large batches, not per-record
- For gzipped annotation files, decompress once and cache in memory

### 7. Split Threads Between Application and Library I/O
```rust
// When using rust-htslib, budget threads for both decompression and application logic
let hts_threads = std::cmp::max(1, total_threads / 4);
let app_threads = total_threads - hts_threads;
reader.set_threads(hts_threads)?;
rayon::ThreadPoolBuilder::new().num_threads(app_threads).build_global()?;
```

> HTSlib BAM decompression is its own form of parallelism. Giving all threads to rayon
> while starving HTSlib (or vice versa) hurts overall throughput.

## Thread Scaling Analysis

Run once the parallelism optimization is in place:

```bash
for threads in 1 2 4 8; do
  echo "Threads: $threads"
  hyperfine --warmup 1 --min-runs 3 \
    "./target/release/<binary> --threads $threads <args>"
done
```

Document the scaling curve. Linear speedup to physical core count is ideal.

- Common bottleneck: if scaling plateaus at 2-4 threads, suspect I/O contention or lock contention.
- For chromosome-level parallelism with indexed BAMs: each thread opens its own `IndexedReader`, distributing chromosomes round-robin. This avoids file-level contention.

## Output

Update `TODO-rewrite-benchmark.md` with final timing numbers.

Maintain an optimization log in the project (e.g. in `TODO-rewrite-benchmark.md`
or a dedicated section). Every loop iteration gets a row:

```markdown
## Optimization Log

| # | Change | Before | After | Delta | Status | Notes |
|---|--------|--------|-------|-------|--------|-------|
| 1 | Add rayon parallelism | 45.0s | 12.0s | -73% | Kept | Linear to 4 threads |
| 2 | Switch to AHashMap | 12.0s | 11.5s | -4% | Kept | |
| 3 | Pre-allocate buffers | 11.5s | 11.2s | -3% | Kept | |
| 4 | Inline hot function | 11.2s | 11.3s | +1% | Reverted | No improvement |
| 5 | (stopped) | 11.2s | — | — | — | Top hotspot < 5%, diminishing returns |
```

**Final entry must document the stopping decision** — why further optimization
is not worth pursuing. Include the final flamegraph's top hotspots.
