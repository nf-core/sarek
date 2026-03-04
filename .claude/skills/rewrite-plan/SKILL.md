---
name: rewrite-plan
description: Build a detailed implementation plan for rewriting pipeline tools in Rust. Use after /rewrite-identify has selected targets, or when the user specifies which tools to rewrite. Covers source code acquisition, test data generation, Rust project scaffolding, phased implementation breakdown, and cost estimation.
model: opus
user-invocable: true
argument-hint: "<target-description> [--pipeline <path>] [--tools <tool1,tool2,...>]"
---

# Plan the Rewrite

Produce a complete, actionable implementation plan that allows an autonomous agent to execute
the rewrite with minimal human intervention. The plan covers source acquisition, test data,
project scaffolding, phased implementation, and cost estimation.

## Arguments

- `$ARGUMENTS` — Description of the rewrite target (e.g., "Rewrite dupRadar and featureCounts")
- `--pipeline <path>` — Path to the pipeline repository
- `--tools <list>` — Comma-separated list of specific tools to rewrite

Read `TODO-rewrite-identify.md` if it exists for context from the identification phase.

## Procedure

### Step 1: Acquire upstream source code

**This is the most critical step. ALWAYS find the source code. NEVER reverse-engineer.**

For each tool being rewritten:

1. Check GitHub/GitLab/Bitbucket for the tool's repository
2. Check package managers: `pip show <pkg>` (PyPI), CRAN package page, Bioconductor
3. Check the tool's documentation/website for a "Source" or "Repository" link
4. Check the pipeline's container Dockerfile for version pins
5. Document the EXACT version being used in the pipeline (from container or conda env)
6. Check out the EXACT version tag/commit used by the pipeline — do not use `main` if the
   pipeline pins a specific release
7. Document the EXACT source code URL and commit hash

If source code truly cannot be found (rare), flag as HIGH RISK and document why.

**Decompiling as last resort:** If source is unavailable (e.g., Java tools distributed as JARs
with no public repo), download the package from conda/bioconda, extract the archive, and use
a decompiler (e.g., CFR for Java `.class` files). This is expensive but sometimes necessary —
it reveals undocumented behaviors that no amount of documentation reading will uncover.
Decompiled code is harder to read, so budget extra time.

### Step 2: Generate reference outputs

Run each upstream tool on test data and capture ALL output files. This is ground truth.

1. **Check pipeline invocation first**: Read the pipeline's module `main.nf` file to see how
   the tool is actually called. Check `ext.args` in config. Check for strandedness/PE handling
   in the module script or templates. The actual invocation often differs significantly from
   the tool's defaults (e.g., special flags, different modes).
2. **Find test data**: Check the pipeline's test profile, CI configuration, or test fixtures.
   If a Seqera Platform run ID is available, use the Seqera MCP to identify input file paths.
3. **Small test data** (required): Must run in < 1 minute. Used for rapid iteration during
   development. If the pipeline has a `test` profile, use those inputs. Ensure the test data
   actually exercises the code paths being validated (e.g., a BAM with 0 mapped reads is
   useless for testing a duplicate-counting algorithm).
4. **Large test data** (required): Realistic-sized input. Used for final validation and
   benchmarking. If a previous pipeline run exists, use those inputs.
5. **Run upstream tools**: Execute each tool with the EXACT same parameters the pipeline uses
   (from step 1). Capture every output file, not just the "main" one. Use Docker containers
   with pinned versions for reproducibility.
6. **Document output format**: For each output file, record:
   - File name pattern (MultiQC looks for specific names — check MultiQC source for search patterns)
   - Format (TSV, CSV, text, binary, plot)
   - Delimiter, header presence, comment lines
   - Column names and types
   - Float precision (significant digits, scientific notation, trailing zeros)
   - NA/NaN/missing value representation
   - Sort order (if any)
   - Line endings (LF vs CRLF)
   - Thousands separators, locale-specific formatting
   - Header strings that downstream tools (especially MultiQC) use for format detection
7. **Identify downstream consumers**: For each output file, trace what reads it:
   - Does MultiQC parse it? Which columns/sections does it actually use?
   - Does another pipeline process consume it?
   - Are there summary statistics derived from it?
   - This determines which outputs are critical to match exactly vs. which can differ

### Step 3: Plan the Rust project

Reference `rust-project-scaffold.md` for the project structure.

Decide:
- **Binary name**: Short, descriptive (e.g., `rustqc`)
- **CLI structure**: Subcommands if multiple modes, flags for tool-specific options
- **Shared data structures**: Gene models, interval trees, BAM record processing
- **Output architecture**: Which files are generated, naming convention
- **Configuration**: YAML config file for toggling outputs, or just CLI flags

If bundling multiple tools into a single binary:
- Design a single-pass BAM iteration that feeds all analyses
- Plan shared annotation parsing (GTF/BED → gene model → per-tool indexes)
- Identify which tools need PE mate reconciliation vs single-read processing

### Step 4: Break into implementation phases

Each phase must produce testable output with clear acceptance criteria.

**Phase 1** (always first): Core algorithm + primary output file matching
- Parse inputs, implement main algorithm, produce the primary output
- Acceptance: primary output matches reference within tolerance

**Subsequent phases**: Additional outputs, secondary analyses, edge cases
- Each phase adds one or two output files
- Each phase has its own acceptance test

**Final phases**: Plots, HTML reports, MultiQC compatibility
- Plots are hardest — save for last. Budget significant time for visual matching.
- Ensure file names match what MultiQC expects (if applicable)
- Generate both SVG and PNG outputs. SVG is more portable for web/reports.

For each phase, estimate:
- Number of sessions (each session = 1-2 hours of agent time)
- Model recommendation (sonnet for mechanical coding, opus for complex algorithms)
- Risk level (low/medium/high)

**Estimation guidance from experience:**
- Simple I/O format matching: 1-2 sessions
- Complex numerical algorithms (statistical estimators, continued fractions, bootstrapping):
  multiply initial estimate by 3-5x. The bugs are in mathematical details.
- Float formatting to match upstream: can burn a full session on its own
- Plot generation to visually match upstream: 2-4 sessions minimum per plot type
- PE read handling with mate reconciliation: always medium or high risk

### Step 5: Plan test infrastructure

- Write validation scripts that compare outputs (save as `scripts/validate_*.py`)
- Plan CI: `cargo fmt --check`, `cargo clippy -- -D warnings`, `cargo test --release`
- Plan integration tests that run the binary and compare against reference outputs

### Step 6: Estimate cost and timeline

- Sessions × cost per session ($5-10 for Sonnet, $20-30 for Opus)
- Total estimated token cost
- Calendar time estimate (assuming N sessions per day)
- Flag phases that might need multiple iterations (complex algorithms, float matching)

### Step 7: Write the plan

Write to `TODO-rewrite-plan.md` using the template from `plan-template.md`.

### Memory budget estimation

For tools processing large files (BAM/CRAM), estimate per-accumulator memory budgets:
- Hash maps keyed by read name or fragment position can grow to 1-2 GB on large BAMs
- Per-gene/per-transcript arrays scale with annotation size
- Coverage arrays scale with genome/exome size
- Document expected memory for each accumulator in the plan
- Include warnings about OOM risks for specific data structures

## Key principles

- Source code URLs are the single most important artifact. Without them, everything is harder.
- Reference outputs are ground truth. Generate them BEFORE writing any Rust code.
- Small test data enables fast iteration. It should run in seconds, not minutes.
- Each phase must be independently testable. Don't design phases that require everything to work.
- Overestimate difficulty. Matching float formatting alone can burn a full session.
- Document upstream tool bugs if found — they must be replicated, not fixed.
- Check the pipeline's actual invocation of each tool before generating reference outputs.
  The pipeline may pass flags that change behavior significantly (PE mode, strandedness, etc.).
- Plan for cross-pipeline reuse from the start. Use a subcommand architecture so the same
  binary can serve multiple pipelines. Shared infrastructure (BAM iteration, CIGAR parsing,
  accumulator framework) should be modular.
- Separate "what to replicate" from "what to improve." First achieve exact output equivalence,
  then optimize. Never try to improve behavior and validate correctness simultaneously.
- For tools producing non-deterministic output (bootstrap CIs, random sampling), define
  correctness criteria upfront: deterministic metrics require exact match, stochastic outputs
  require statistical equivalence. Document which outputs will differ and why.
