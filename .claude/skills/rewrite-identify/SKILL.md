---
name: rewrite-identify
description: Analyze a Nextflow pipeline and rank tools by rewrite ROI (effort vs improvement). Use when starting a new rewrite project, evaluating which pipeline tools are worth rewriting, finding performance bottlenecks, or deciding what to rewrite first. Accepts trace files, Seqera Platform run IDs, or repo-only analysis.
model: opus
user-invocable: true
argument-hint: "<pipeline-path> [--trace <trace-file>] [--run-id <seqera-run-id>]"
---

# Identify Rewrite Targets

Analyze a Nextflow pipeline to find the highest-value tools to rewrite in a compiled language.
Rank candidates by effort:improvement ratio, identify bundling opportunities, and produce a
structured report.

## Arguments

- `$ARGUMENTS` — Path to the pipeline repository (required)
- `--trace <file>` — Nextflow execution trace file (`execution_trace.txt`)
- `--run-id <id>` — Seqera Platform run ID (uses Seqera MCP to fetch metadata)

## Procedure

### Step 1: Gather runtime data

Use ONE of these approaches depending on what's available:

**With Seqera MCP run ID** (best):
- Use `search_seqera_api` with query "list tasks for workflow run" to find the API
- Call `call_seqera_api` with the run ID to get task-level metadata
- Extract: task name, wall time, CPU usage, peak memory, read/write bytes
- This gives the most complete picture including resource usage patterns

**With trace file** (good):
- Parse the tab-separated trace file. Key columns:
  - `name`: process name (includes sample ID)
  - `realtime`: wall clock time per task
  - `%cpu`: CPU utilization
  - `peak_rss`: peak memory usage
  - `rchar` / `wchar`: bytes read/written
- Aggregate by process name (strip sample-specific parts)
- Calculate: total time per process, mean time per task, task count

**Repo-only** (acceptable):
- Inspect `main.nf`, `workflows/`, `modules/`, `subworkflows/`
- Identify tools from process definitions (container directives, script blocks)
- Estimate runtime from: tool type (aligner=slow, QC=medium, stats=fast), input data type (BAM=large, VCF=small), per-sample vs global

### Step 2: Inventory all tools

For each process/tool in the pipeline, document:
- Tool name and version
- Language (Python, R, shell, C/C++, Java, Perl)
- Source code availability (check GitHub, PyPI, CRAN, Bioconductor)
- Source code availability and hosting (GitHub/GitLab/Bitbucket/other — note if source is only available as compiled JARs or binaries, which significantly increases difficulty)
- Source code URL if found
- Input files (type and typical size)
- Output files (type and count)
- Whether the tool produces plots/visualizations (these require significant extra effort to match)
- What downstream tools consume the output (especially MultiQC — check its search patterns)
- Whether it runs per-sample or globally
- Whether it's a bottleneck (top 20% by runtime)

### Step 3: Score rewrite candidates

Rate each tool on these criteria:

**Difficulty** (low / medium / high):
- Low: Python/R tool with simple algorithm, source available, well-documented format
- Medium: Complex algorithm, multiple output files, PE read handling, coordinate math
- High: C/C++ tool (marginal speedup from rewrite), complex statistical model, GPU code

**Impact** (based on runtime data):
- Per-sample tools that run N times get N× multiplier on savings
- Tools processing large files (BAM, FASTQ) benefit most from single-pass architecture
- Tools with high I/O (read same file multiple times) benefit from bundling

**Bundling opportunities**:
- Tools reading the SAME input files (e.g., multiple BAM QC tools) can share a single pass
- Sequential tool chains where output of A feeds only into B can be merged
- Tools requiring the same annotation parsing (GTF/GFF) can share the parsed data structure

### Architecture considerations

- Single-pass architecture is the killer optimization for bundled tools: iterate the input file (BAM/CRAM) once, dispatching to all analysis accumulators simultaneously. This eliminates N-1 file passes.
- ~80-90% of wall time for QC tools is in BAM I/O. The individual computations are negligible. Architecture (single-pass) matters more than micro-optimization.
- Consider the accumulator pattern: each tool defines a struct with `process_read()` and `merge()` methods. All accumulators run in the same file pass, with per-thread instances merged after parallel processing.

### Step 4: Rank and recommend

Sort candidates by `(predicted_savings × task_count) / difficulty_score`.

Map difficulty to numeric values: low=1, medium=3, high=9.

For bundleable tools, calculate the bundled score: `sum(individual_savings) / max(individual_difficulties)` — bundling amortizes the hardest tool's difficulty across all bundled tools.

Group bundleable tools into rewrite targets (a single Rust binary can replace multiple tools).

For each recommended target:
- Predicted difficulty (low/medium/high)
- Predicted time savings (per task and total across pipeline run)
- Predicted token cost (rough estimate: low=$50-200, medium=$200-500, high=$500-2000)
- Risk factors (complex algorithms, poor documentation, no source code)

**Cost estimation caveats**:
- Initial estimates are typically 3-5x too optimistic for tools with complex numerical algorithms (statistical estimators, continued fractions, bootstrapping)
- Plot generation adds 1-3 sessions per tool
- Float formatting and output matching alone can consume a full session

### Step 5: Write report

Write the analysis to `TODO-rewrite-identify.md` using the template from `analysis-template.md`.

## Key principles

- Per-sample tools are almost always higher value than global tools (N× multiplier)
- I/O-bound tools benefit most from bundling (shared BAM iteration)
- Python/R tools with simple algorithms are the lowest hanging fruit
- Tools without available source code are significantly harder — flag as high risk
- Don't underestimate the value of eliminating staging/unstaging of large intermediate files
- A single Rust binary replacing 5 tools doesn't just save compute — it eliminates 4 task submissions, 4 container pulls, and 4 sets of file staging
