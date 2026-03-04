---
name: rewrite-implement
description: Write the Rust implementation guided by the rewrite plan. Use when ready to write code after /rewrite-plan has produced TODO-rewrite-plan.md, or when resuming implementation work on an existing rewrite project. Each invocation handles one focused session (one phase or sub-phase).
model: sonnet
user-invocable: true
argument-hint: "<plan-file-or-description>"
---

# /rewrite-implement — Build the Rewritten Tool

You are implementing a bioinformatics tool rewrite in Rust. You have a plan document
(`TODO-rewrite-plan.md` or as provided in `$ARGUMENTS`) that specifies what to build,
the upstream source code URLs, test data locations, and phased acceptance criteria.

Your job is to write correct, efficient Rust code that produces **byte-identical output**
(or numerically equivalent within tolerance) to the upstream tool.

## Pre-Implementation Checklist

Before writing ANY Rust code, complete these steps:

1. **Read the plan** — Load `TODO-rewrite-plan.md` (or the provided plan). Understand the
   phases, acceptance criteria, and source code URLs.
2. **Verify test data** — Confirm small and large test inputs exist at the documented paths.
3. **Generate reference outputs** — Run the upstream tool on the small test input. Capture
   ALL output files. Store in a `tests/expected/` or `benchmark/reference/` directory.
   These are your ground truth. NEVER modify them by hand.
   - **Check the pipeline's actual invocation first.** Read the module's `main.nf` and any
     `ext.args` to see exactly which flags, modes (e.g., paired-end, strandedness), and
     annotation parameters the pipeline passes. These can change output significantly.
   - **Ensure test data exercises relevant code paths.** A BAM with 0 mapped reads is
     useless for testing a duplicate-counting algorithm. Verify that small test inputs
     contain enough variety (mapped/unmapped, paired/single, supplementary, etc.) to
     exercise the features the pipeline actually uses.
4. **Write validation scripts** — Before writing Rust, write Python or shell scripts that
   compare your output against the reference. Save as `scripts/validate_*.py` or
   `scripts/validate_*.sh`. These scripts persist for CI and future validation passes.
5. **Verify upstream source code** — Fetch the upstream source from the URL in the plan.
   Read the relevant algorithms. Always prefer finding the original source code.
   - **Decompilation is a last resort.** If source is truly unavailable (e.g., only
     distributed as compiled `.class` files), decompiling (e.g., CFR for Java) is
     acceptable. Decompiled code often reveals undocumented behaviors and edge cases.
     Budget extra time — decompiled output lacks comments and uses mangled names.

## Session Discipline

This skill is designed for focused, bounded sessions. Each invocation should accomplish
ONE phase or sub-phase from the plan.

### Starting a Session

1. Read the plan (`TODO-rewrite-plan.md`)
2. Read the previous session log (`TODO-rewrite-session-N.md`) if it exists
3. Identify the next task to work on
4. Create a new session log file: `TODO-rewrite-session-{N+1}.md`

### During a Session

- Focus on ONE task from the plan
- Build incrementally: parse input -> core algorithm -> format output -> edge cases
- Compile early and often — use `cargo check` after every few edits to catch errors before they compound
- After every significant change, run:
  ```bash
  cargo fmt && cargo clippy -- -D warnings && cargo test
  ```
- Test against reference outputs frequently using your validation scripts
- Use the TaskCreate / TaskUpdate tools to track sub-tasks within the session

### Ending a Session

1. **Commit work** — Always commit before ending. Use descriptive messages.
2. **Update the plan** — Mark completed phases, note any blockers or discoveries.
3. **Write session log** — Document what was done, what works, what's next, any blockers.
4. **Estimate remaining work** — Update session/cost estimates in the plan.

## Git Commit Discipline

**Commit early, commit often.** Long-running agentic sessions risk crashes and disconnects.
Never accumulate more than ~2 hours of uncommitted work.

- Commit after every successful phase or milestone
- Commit after generating reference outputs
- Commit after writing validation scripts
- Commit after each test passes for the first time
- Use descriptive commit messages: `"Phase 1: core counting algorithm matches reference"`
- If a session ends unexpectedly, the git history should preserve all meaningful progress

## Coding Practices

### Project Structure
```
src/
  main.rs           # Entry point, CLI parsing
  <module>.rs       # Per-tool or per-feature modules
tests/
  integration_test.rs
  data/             # Small test inputs
  expected/         # Reference outputs from upstream tool
scripts/
  validate_*.py     # Validation scripts
benchmark/
  input/            # Test data
  reference/        # Upstream tool outputs for comparison
```

### Rust Conventions
- `anyhow::Result<T>` for all fallible functions, propagate with `?`
- `#[derive(Debug)]` on all structs
- `///` doc comments on all public items
- `//!` module doc comment at top of every file
- No `unwrap()` or `expect()` in production code (OK in tests)
- Use `IndexMap` when insertion order matters, `HashMap` otherwise
- `u64` for counts/positions, `f64` for metrics

### Recommended Crates
- **`clap`** v4 (derive API) — CLI argument parsing
- **`rust-htslib`** — BAM/CRAM I/O (use static linking for portable binaries)
- **`plotters`** — Chart generation (PNG + SVG backends)
- **`rayon`** — Data parallelism (parallel iterators)
- **`coitrees`** — Interval queries (cache-oblivious interval trees, fast for genomic overlaps)
- **`anyhow`** — Error handling with context
- **`indexmap`** — Insertion-order-preserving maps and sets
- **`rand`** + **`rand_chacha`** — Reproducible random sampling (seed with ChaCha for determinism)

### Build & Check Cycle
```bash
cargo fmt                        # Format
cargo clippy -- -D warnings      # Lint (zero warnings)
cargo test                       # Unit + integration tests
cargo build --release            # Release build
./target/release/<tool> <args>   # Run on test data
python scripts/validate_*.py     # Compare against reference
```

## Plot Generation

See **[plot-generation.md](plot-generation.md)** for comprehensive guidance on generating
plots with `plotters`, including workarounds for known issues (blocky lines, broken rotated
labels, missing chart types), resolution-aware scaling patterns, R parameter matching, and
density estimation.

## Working with Upstream Source Code

### Finding Source Code
The plan should have URLs. If not, find them:
- **Python tools**: `pip show <package>` -> Homepage/Source URL
- **R packages**: Check CRAN (`https://cran.r-project.org/package=<name>`) or
  Bioconductor (`https://bioconductor.org/packages/<name>`)
- **Command-line tools**: Check GitHub/GitLab, or the tool's documentation site

### Reading Source Code
- Read the upstream source DIRECTLY. Fetch from GitHub/Bitbucket raw URLs.
- Use the `explore` subagent type for searching large codebases
- When behavior is unclear, write a small test case and run the upstream tool
- Document any upstream bugs you discover — you must replicate them to match output

### CRITICAL: Match Upstream Exactly
- NEVER assume the upstream tool is wrong
- NEVER "improve" on the upstream behavior
- If the upstream tool has an off-by-one error, replicate it
- If the upstream tool formats floats with trailing zeros, do the same
- If the upstream tool sorts unstably, match its sort order
- If the upstream tool has actual bugs (e.g., CIGAR parsing errors), you must REPRODUCE
  them to match output. Document the bug with a `// BUG COMPAT:` comment referencing the
  upstream source file and line number.
- **Check how the pipeline invokes the tool** — flags like paired-end mode, strandedness,
  and annotation parameters change output significantly. Match the pipeline's invocation.
- **For paired-end handling:** verify what counting strategy is used (union vs intersection
  of gene hits), whether both mates are counted, and which BAM flags control PE/SE mode.
- The ONLY acceptable differences are documented in the plan's acceptance criteria

## Cost Efficiency Rules

These rules minimize token spend and maximize progress per session:

1. **Use `explore` subagents** for codebase searches instead of sequential file reads
2. **Write comparison scripts** (Python/shell) instead of doing interactive diffs.
   "Write a script that compares these two TSV files" = 1 prompt.
   Iterating on the comparison interactively = 10 prompts.
3. **When stuck on a discrepancy:**
   - Add diagnostic logging or counters
   - Dump data to files
   - Write an analysis script to compare the dumps
   - Make ONE informed fix based on the analysis
   - Do NOT iterate build-run-check blindly
4. **Batch related changes** — Don't make one edit, build, test, repeat 10 times.
   Plan the change, make all edits, then build and test once.
5. **Know when to escalate** — If a discrepancy requires deep algorithmic understanding
   of the upstream tool, note it in the session log and recommend an Opus session.
6. **Reuse validation scripts** — Once written, run them after every change. Don't
   manually inspect output files.

## Context and Memory Management

- **In long agentic sessions, context limits are a real constraint.** Prune and distill
  aggressively after completing each phase. Summarize findings to session logs rather than
  keeping large code blocks in working memory.
- **Re-read TODO/plan files from disk** rather than relying on earlier context. Files on disk
  are the source of truth — your memory of them may be stale or truncated.
- **Focus on ONE phase at a time.** Complete it, validate it, commit it, then move on.
  Do not interleave work across phases.
- **Memory budgets matter at runtime too.** Calculate expected `HashMap` sizes for large
  inputs before implementing. A `HashMap` keyed by read name on a 10 GB BAM can easily
  consume 1–2 GB of RAM. Consider streaming approaches, `BTreeMap` for sorted access,
  or sharding strategies when memory is a concern.

## Handling Discrepancies

When your output doesn't match the reference:

0. **Enumerate ALL differences first** — Before fixing anything, systematically audit ALL
   differences between your implementation and the reference. Write them to a tracking file
   (e.g., `TODO-fix-<tool>.md`). Group by category (formatting, numeric, missing data, etc.).
   Fix in groups by category. This prevents whack-a-mole debugging where fixing one thing
   breaks another.
1. **Quantify the difference** — Run validation script, get exact numbers
2. **Classify** — Is it formatting? Rounding? Wrong value? Missing data?
3. **Locate** — Which specific output, which line, which column?
4. **Trace upstream** — Read the upstream source for that specific computation
5. **Hypothesize** — Form ONE hypothesis for the cause
6. **Test** — Write a targeted diagnostic (counter, log, dump) to test the hypothesis
7. **Fix** — Make ONE informed change
8. **Verify** — Run validation. If fixed, commit. If not, go to step 5.

Do NOT skip steps 3-5 and jump straight to trying random fixes.

## Reference Material

Read these supporting files **when the trigger condition applies**:

- **`bioinformatics-tips.md`** — Read at the START of any session working with BAM files,
  GTF/GFF annotations, or interval operations. Contains coordinate system gotchas, CIGAR
  parsing bugs to replicate, PE read edge cases, and float formatting lessons.
- **`session-management.md`** — Read when starting your FIRST session to understand the
  session brief/log format. Skim when starting subsequent sessions.
- **`wave-integration.md`** — Read only when you need to cross-compile the binary for
  Linux x86_64 (e.g., testing on a cloud executor). Not needed for local development.

## Output

- Update `TODO-rewrite-plan.md` with progress after each session
- Write `TODO-rewrite-session-{N}.md` for each session
- Commit all code changes with descriptive messages
