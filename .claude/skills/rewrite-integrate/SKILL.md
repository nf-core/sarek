---
name: rewrite-integrate
description: Wire the rewritten tool into a Nextflow pipeline with Wave module binaries. Use after the rewrite passes audit to integrate it into the pipeline with a toggleable process (original tool preserved as fallback). Covers cross-compilation, module structure, Wave config, testing, and MultiQC compatibility.
model: sonnet
user-invocable: true
argument-hint: "<pipeline-path> --tool <binary-name>"
---

# Pipeline Integration

Wire the rewritten tool into a Nextflow pipeline using Wave module binaries.
The rewritten tool runs by default, but the original tool remains available via a
parameter toggle.

**Arguments:**
- `$ARGUMENTS` — Path to the Nextflow pipeline repository
- `--tool <name>` — Name of the compiled binary

Read `nextflow-patterns.md` in this skill directory for all DSL2 code patterns
(module structure, process definition, toggle pattern, Wave config, versions YAML,
nf-test, preview/lint commands, resume strategy).

## Pre-Integration Checklist

Before starting integration:

1. **Binary exists and is tested:**
   ```bash
   ls ./target/release/<binary>
   file ./target/release/<binary>  # Should be ELF 64-bit for Linux
   ```

2. **Pipeline runs in current state:**
   ```bash
   cd <pipeline-path>
   nextflow run . -profile test -preview  # Cheap syntax check first
   ```

3. **Nextflow lint (if available):**
   ```bash
   nextflow lint .  # Check for common issues
   ```

## Cross-Compilation

The binary must run on Linux x86_64 (standard cloud executor target).
If developing on macOS, cross-compile:

### Option A: Docker build (recommended)
```bash
docker run --rm -v $(pwd):/src -w /src rust:latest \
  bash -c "apt-get update && apt-get install -y cmake libz-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev && cargo build --release"
```

### Option B: Cross crate
```bash
cargo install cross
cross build --release --target x86_64-unknown-linux-gnu
```

### Option C: Zig as linker
```bash
cargo install cargo-zigbuild
cargo zigbuild --release --target x86_64-unknown-linux-gnu
```

Verify the binary is statically linked or only depends on common system libraries:
```bash
# Inside a Linux environment
ldd ./target/release/<binary>
# Ideally shows only libc, libm, libpthread, libdl
```

## Module Structure

Create the module directory in the pipeline using the structure from `nextflow-patterns.md`.

Copy the binary:
```bash
mkdir -p <pipeline-path>/modules/local/<tool>/resources/usr/bin/
cp ./target/x86_64-unknown-linux-gnu/release/<binary> \
   <pipeline-path>/modules/local/<tool>/resources/usr/bin/
chmod +x <pipeline-path>/modules/local/<tool>/resources/usr/bin/<binary>
```

**Single-toggle pattern:** Use a single pipeline parameter (e.g., `--use_rewritten_tool`)
that auto-disables all replaced upstream processes. Don't require users to set multiple
flags (`--skip_tool_a`, `--skip_tool_b`, etc.).

## Process Definition (main.nf)

Write `modules/local/<tool>/main.nf` using the process definition pattern from
`nextflow-patterns.md`. Key requirements:

1. **Output channels must match the original process** so downstream works with either
2. **Use the binary name directly** — it's on PATH via moduleBinaries
3. **Include resource labels** appropriate for the tool
4. **Tag for logging** so pipeline logs show which sample is running

## Pipeline Configuration

Apply the Wave configuration and toggleable process pattern from `nextflow-patterns.md`.

**IMPORTANT:** Keep ALL original processes intact. Do not delete or modify them.
The toggle allows easy rollback and A/B comparison.

## Nextflow/Groovy Pitfalls

**String escaping in output declarations:**
Parentheses, brackets, and other special characters in Nextflow `path()` output
declarations are parsed as Groovy syntax, causing `LexerNoViableAltException`.
If upstream tools produce filenames with special characters
(e.g., `coverage_profile_(total).txt`), use glob patterns in all Nextflow-facing
declarations (`*total*`) and only emit the literal filename from the Rust binary.

**Boolean parameters with nf-schema:**
Default-false boolean params need explicit `"default": false` in the JSON schema.
The Nextflow CLI passes `--param_name` as the string `"true"`, which the JSON Schema
validator rejects as a type mismatch. Test boolean params via a params file, not bare
CLI flags.

**Conditional process execution:**
Using `if (params.flag) { PROCESS(...) }` in DSL2 causes `MissingPropertyException`
for `.out.channels` because Nextflow evaluates channel wiring at compile time. Instead,
always call every process but feed `Channel.empty()` when disabled:
```nextflow
ch_input = params.skip_tool ? Channel.empty() : ch_real_input
TOOL(ch_input)
```

**Disabling processes inside subworkflows:**
When replacing tools called inside subworkflows (e.g., stats tools inside a
markduplicates subworkflow), you can't gate them from the calling workflow. Use
Nextflow's `ext.when` config mechanism:
```nextflow
process {
    withName: 'SUBWORKFLOW:INNER_PROCESS' {
        ext.when = { !params.use_rewritten_tool }
    }
}
```

## Testing

### Step 1: Syntax validation (cheap, fast)
```bash
cd <pipeline-path>
nextflow run . -profile test -preview
```

**Always run `-preview` after every pipeline change.** It's nearly instant and catches
syntax, type, and wiring errors before a full run.

### Step 2: Full test run with rewritten tool
```bash
nextflow run . -profile test,docker -resume \
    --use_rewritten_<tool> true \
    --outdir results_rewritten
```

**Always use `-resume`** to cache upstream steps that haven't changed.
If previous run data is available (e.g., via Fusion mount), most steps will be cached
and only the new/changed process will actually run.

### Step 3: Full test run with original tool (for comparison)
```bash
nextflow run . -profile test,docker -resume \
    --use_rewritten_<tool> false \
    --outdir results_original
```

### Step 4: Compare outputs
```bash
# Compare all output files
diff -rq results_rewritten/ results_original/

# For numeric files, use the validation scripts from /rewrite-validate
python scripts/validate_outputs.py results_rewritten/ results_original/
```

### Step 5: MultiQC compatibility
```bash
# Run MultiQC on rewritten tool output
multiqc results_rewritten/ -o results_rewritten/multiqc/

# Compare with original
multiqc results_original/ -o results_original/multiqc/

# Check that the same modules are found and reports look correct
```

Ensure output file names and formats match what MultiQC modules expect.
If the rewritten tool produces files with different names, add symlinks or
rename in the process script (see output channel compatibility in `nextflow-patterns.md`).

**Critical MultiQC details:**
- Read MultiQC's source code search patterns for every tool being replaced. MultiQC
  detects file types by scanning for specific header strings (e.g.,
  `This file was produced by samtools stats`). Reproduce these trigger strings verbatim.
- Check which columns/sections MultiQC actually parses — it may only read columns 0
  and 1, ignoring the rest. This determines which parts of the output are critical to
  match exactly.
- Verify output file naming matches upstream conventions exactly (extensions, directory
  structure). MultiQC uses filename patterns for detection.

## Non-Wave Fallback

For environments without Wave, document manual Docker image build:

```dockerfile
FROM ubuntu:22.04
COPY modules/local/<tool>/resources/usr/bin/<binary> /usr/local/bin/
RUN chmod +x /usr/local/bin/<binary>
```

```bash
docker build -t <tool>:latest -f Dockerfile.rewritten .
```

Add to nextflow.config:
```nextflow
process {
    withName: 'REWRITTEN_TOOL' {
        container = '<tool>:latest'
    }
}
```

## Common Integration Issues

- **`.nf-test/`, `test_output/`, `work/` directories:** Add these to `.gitignore`
  immediately when setting up integration tests.
- **`versions.yml`:** nf-core modules require a `versions.yml` output. Emit your tool's
  version in the standard format.
- **`nf-core modules install` must be run sequentially:** The underlying git operations
  clash if parallelized. Install one module at a time.
- **Strandedness and PE logic:** These often live in the module's `main.nf` or templates
  (reading `meta.strandedness`, `meta.single_end`), NOT in `ext.args`. Trace the full
  parameter path when replacing a module.
- **NEVER make security-impacting changes** (repo visibility, access tokens) without
  explicit user approval.

## Output

Write `TODO-rewrite-integrate.md` documenting:
- Module structure created
- Config changes made
- Test run results (both tools)
- Output comparison results
- MultiQC compatibility status
- Any issues found and fixes applied
