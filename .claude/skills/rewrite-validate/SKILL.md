---
name: rewrite-validate
description: Exhaustive comparison of rewritten tool output against upstream reference. Use after implementing a phase to verify correctness, when outputs don't match and you need systematic diagnosis, or before benchmarking to confirm the code is correct. Compares file structure, numeric values, string formatting, and edge cases.
model: sonnet
user-invocable: true
argument-hint: "<output-dir> --reference <reference-dir>"
---

# /rewrite-validate — Verify Output Correctness

You are validating the output of a rewritten bioinformatics tool against the upstream
tool's reference output. Your job is to perform an exhaustive, rigorous comparison
and document every difference.

**Arguments:** `$ARGUMENTS` should specify the output directory and reference directory,
e.g., `output/ --reference benchmark/reference/`. If not provided, look for
`TODO-rewrite-plan.md` to find the expected paths.

## Core Principle

**The upstream tool output is ground truth. Your code is wrong until proven otherwise.**

NEVER rationalize a difference as an "upstream bug" without conclusive evidence from
the upstream source code. NEVER accept a difference as "close enough" unless the plan's
acceptance criteria explicitly allow it.

## Validation Process

### Step 1: File-Level Comparison

```bash
# List all files in both directories
diff <(cd "$OUTPUT_DIR" && find . -type f | sort) \
     <(cd "$REFERENCE_DIR" && find . -type f | sort)
```

Check:
- Same number of files
- Same file names (case-sensitive)
- No unexpected extra files in output
- No missing files from output
- File sizes are in the same ballpark (> 50% difference is suspicious)

### Step 2: Content Comparison by File Type

**TSV/CSV files:**
- Compare headers (column names, order, delimiter)
- Compare row counts
- Compare sort order (is the data sorted the same way?)
- Column-by-column comparison:
  - String columns: exact match
  - Integer columns: exact match
  - Float columns: numeric tolerance (see below)
- Check for trailing whitespace, empty trailing lines, BOM markers

**Text reports (e.g., `*_results.txt`):**
- Line-by-line comparison
- For lines with numbers: extract and compare numerically
- For label lines: exact string match
- Check section ordering

**Plot files (PNG/SVG/PDF):**
- File exists and is non-empty
- Approximate size comparison (within 2x)
- Verify plot dimensions match upstream tool defaults (e.g., R defaults to 7×7 inch square plots)
- For SVG: parse data values from `<path>`, `<polyline>`, and `<polygon>` elements; compare axis labels and tick values as text; verify color hex values (e.g., `fill="#FF0000"`) against the upstream source's palette definitions
- For PNG: pixel-perfect matching is unrealistic when using different rendering engines — visual side-by-side (A/B) comparison is acceptable
- Check that density estimation parameters (bandwidth, kernel), color palettes, and statistical computations behind the visualization match the upstream implementation
- Note: plot validation is primarily visual — automate data extraction where possible but expect manual review for final sign-off

**Binary files:**
- Checksum comparison (md5/sha256)

### Step 3: Numeric Tolerance

For floating-point comparisons, use this hierarchy:

1. **Exact string match** — Preferred. If the formatted string is identical, done.
2. **Absolute tolerance** — `|a - b| < 1e-10` — For values near zero.
3. **Relative tolerance** — `|a - b| / max(|a|, |b|) < 1e-6` — For larger values.
4. **Significant digits** — Compare first N significant digits (N from upstream format).

Apply the most restrictive tolerance that passes. Document which tolerance was needed.

### Step 4: Classify Differences

For each difference found, classify it:

| Category | Description | Action |
|----------|-------------|--------|
| **PASS** | Identical output | None |
| **FORMAT** | Whitespace, trailing zeros, line endings | Document, likely acceptable |
| **ROUND** | Within numeric tolerance | Document, check acceptance criteria |
| **VALUE** | Different numeric value beyond tolerance | **FIX REQUIRED** |
| **MISSING** | Data present in reference but not output | **FIX REQUIRED** |
| **EXTRA** | Data in output but not reference | **FIX REQUIRED** |
| **ORDER** | Same data, different sort order | Check if upstream sort is stable |

### Step 5: Non-Deterministic Output Handling

Some bioinformatics tools produce stochastic output (bootstrap confidence intervals,
random sampling, Monte Carlo methods, EM convergence). Validate these differently:

- **Different RNG implementations produce different sequences** even with the same seed.
  For example, Rust's `rand` crate defaults to ChaCha8 while C/C++ commonly uses
  Mersenne Twister. Do not expect bitwise-identical stochastic results.
- **For stochastic outputs:** validate statistical properties (median, CI width,
  distribution shape, convergence behavior) rather than exact values. Run multiple
  iterations if needed to confirm the distribution is equivalent.
- **For deterministic metrics in the same output file:** still require exact match.
  Many tools mix deterministic counts with stochastic estimates in a single file —
  validate each column/field according to its nature.
- **Document which outputs are deterministic vs stochastic** and what equivalence
  criteria apply for each. Record this in the validation report under
  "Known Acceptable Differences."

## Writing Validation Scripts

Write reusable validation scripts and save them in the project:

```python
#!/usr/bin/env python3
"""Compare rewritten tool output against upstream reference."""

import sys
import csv
import math

def compare_float(a, b, abs_tol=1e-10, rel_tol=1e-6):
    """Compare floats with absolute and relative tolerance."""
    if a == b:
        return True
    if math.isnan(a) and math.isnan(b):
        return True
    diff = abs(a - b)
    if diff < abs_tol:
        return True
    denom = max(abs(a), abs(b))
    if denom > 0 and diff / denom < rel_tol:
        return True
    return False

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", help="Directory with rewritten tool output")
    parser.add_argument("reference_dir", help="Directory with upstream tool output")
    args = parser.parse_args()

    # Compare each expected output file
    # Example: compare_tsv(args.output_dir + "/counts.tsv", args.reference_dir + "/counts.tsv")
    print("Validation complete")
```

These scripts should:
- Be self-contained (no external dependencies beyond Python stdlib)
- Print clear pass/fail per file
- Print specific discrepancies with line/column numbers
- Return non-zero exit code on failure
- Be saved in `scripts/` directory for reuse

## Common Validation Pitfalls

**Coordinate system mismatches:** GTF uses 1-based inclusive coordinates, BAM uses
0-based half-open, and some interval tree libraries use closed intervals. An off-by-one
in coordinate conversion can cause false-positive overlaps across the entire genome.
Always verify coordinate conventions on both sides before debugging count mismatches.

**PE fragment counting:** Verify which BAM flags the upstream tool uses for PE/SE
detection (`proper_pair` 0x2 vs `paired` 0x1), whether both mates are counted or only
one per fragment, and how orphaned mates (mate unmapped) are handled. These semantics
vary widely between tools and are a frequent source of 2× count differences.

**Chromosome naming:** BAM files may use `1, 2, X` while annotation files use
`chr1, chr2, chrX`. Mitochondrial naming varies (`MT`, `chrM`, `M`). Always check for
naming discrepancies first when counts don't match — a name mismatch silently produces
zero counts for the affected chromosomes.

**Different definitions of "total reads":** `samtools flagstat` counts ALL records
including secondary and supplementary alignments, while many tools exclude them.
Verify exact counting semantics per the upstream tool's source code before comparing
totals.

**Test data must exercise relevant code paths:** A test BAM with 0 mapped reads tells
you nothing about a counting algorithm. Verify test data is appropriate for the
functionality being validated — check that test inputs contain PE reads, multi-mappers,
overlapping features, etc. as needed.

**Stderr/stdout validation:** Some downstream tools or pipeline processes parse stderr.
Verify whether the upstream tool produces diagnostic output on stderr that consumers
depend on, and replicate it.

**Tests that silently pass on missing output are worse than no tests:** Always assert
file existence BEFORE comparing contents. Never use conditional guards (e.g.,
`if os.path.exists(f): compare(f)`) that skip assertions when output is missing —
this hides complete failures as passes.

## Validation Report

Write `TODO-rewrite-validate.md` with:

```markdown
# Validation Report — {Date}

## Summary
- Files compared: {N}
- Passed: {N}
- Failed: {N}
- Formatting differences: {N}

## Per-File Results

| File | Status | Notes |
|------|--------|-------|
| output.tsv | PASS | Exact match |
| summary.txt | FORMAT | Trailing whitespace on line 5 |
| counts.tsv | VALUE | Column 3, row 12: got 0.5432, expected 0.5431 |

## Specific Discrepancies

### {filename}
- Line {N}, Column {N}: expected `{value}`, got `{value}`
- Category: VALUE
- Upstream source reference: {file}:{line} — {explanation}
- Suggested fix: {description}

### Known Acceptable Differences
| Output | Difference | Reason | Evidence |
|--------|-----------|--------|----------|
| ci_lower.tsv | Values differ | Different RNG implementation | Statistical distribution equivalent |

## Acceptance Criteria Check
- [ ] All VALUE differences resolved
- [ ] All MISSING/EXTRA differences resolved
- [ ] FORMAT differences documented and accepted
- [ ] ROUND differences within plan's tolerance
```

## Reference Material

See `validation-checklist.md` in this skill directory for an exhaustive list of
things to check.
