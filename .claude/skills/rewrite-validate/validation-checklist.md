# Validation Checklist

Exhaustive list of things to check when comparing rewritten tool output against
upstream reference. Check ALL items, even the ones that seem trivial.

## File Structure
- [ ] Same number of output files
- [ ] Same file names (case-sensitive)
- [ ] Same directory structure (subdirectories match)
- [ ] No unexpected extra files
- [ ] No missing files
- [ ] File sizes within reasonable range (> 50% difference is suspicious)
- [ ] File permissions don't matter (but note if executable bit is set)
- [ ] Tests must assert file existence BEFORE comparing contents — tests that silently pass on missing output are worse than no tests

## Text Encoding
- [ ] UTF-8 encoding (no BOM marker `EF BB BF`)
- [ ] Line endings: LF (`\n`) not CRLF (`\r\n`)
- [ ] No null bytes
- [ ] No trailing whitespace on lines (unless upstream has it — then match it)
- [ ] Trailing newline at end of file (or not — match upstream)

## TSV/CSV Files
- [ ] Correct delimiter (tab vs comma vs space)
- [ ] Header line present/absent (match upstream)
- [ ] Column names match exactly (case-sensitive)
- [ ] Column order matches exactly
- [ ] Row count matches
- [ ] Sort order matches (check if sort is stable in upstream)
- [ ] Comment lines (lines starting with `#`) match

## Numeric Formatting
- [ ] Integer formatting: no leading zeros, no trailing decimal point
- [ ] Float significant digits: match upstream (R default: 15 sig digits)
- [ ] Scientific notation format: `1.23e+05` vs `1.23E5` vs `123000`
- [ ] Trailing zeros: `1.0` vs `1` vs `1.00` — match upstream exactly
- [ ] Negative zero: `-0` vs `0` — match upstream
- [ ] Very small numbers: check for underflow to zero
- [ ] Very large numbers: check for overflow or scientific notation switch
- [ ] NaN representation: `NA` vs `NaN` vs `nan` vs `N/A` vs empty string
- [ ] Infinity representation: `Inf` vs `inf` vs `Infinity`
- [ ] Percentage formatting: `50%` vs `50.0%` vs `0.5` — match upstream

## String Formatting
- [ ] Quoting: `"value"` vs `value` — match upstream
- [ ] Escaping: embedded quotes, tabs, newlines
- [ ] Padding/alignment: fixed-width columns if upstream uses them
- [ ] Case sensitivity: don't change case of any string values

## Report Files
- [ ] Section headers match
- [ ] Section order matches
- [ ] Spacing between sections matches
- [ ] Indentation matches
- [ ] Label text matches exactly
- [ ] Separator lines (e.g., `====` or `----`) match
- [ ] Stderr/stdout: some downstream tools (especially MultiQC) parse specific header strings. Verify trigger strings are reproduced verbatim.

## Plot Validation
- [ ] File exists and is non-empty
- [ ] Correct format (PNG/SVG/PDF as expected)
- [ ] Reasonable file size (not 0 bytes, not orders of magnitude different)
- [ ] Dimensions match upstream defaults (width/height in pixels or viewBox units)
- [ ] Color palette matches source code / upstream defaults
- [ ] Axis labels, tick marks, and titles are correct
- [ ] SVG: parse actual data values from `<path>`, `<polyline>`, `<rect>`, etc. elements and compare against expected numeric values
- [ ] PNG: visual side-by-side comparison (pixel-perfect match is unrealistic due to renderer differences)
- [ ] Density estimation parameters (bandwidth, kernel type) must match the upstream algorithm exactly — small differences cascade into visible plot divergence
- [ ] Regression testing: use A/B comparison with rated scoring (1-10 scale) for iterative improvement when exact match is infeasible

## Edge Cases to Test
- [ ] Empty input (0 reads/features)
- [ ] Single record input
- [ ] All-duplicate input
- [ ] Mixed chromosome input (if tool is chromosome-aware)
- [ ] Negative strand genes (if applicable)
- [ ] Overlapping features (if applicable)
- [ ] Very long sequences / very short sequences

## Stochastic Output Validation
Some tools produce non-deterministic output (bootstrap CIs, random sampling, Monte
Carlo methods, etc.). Different RNG implementations produce different exact values
even with the same seed.

- [ ] Identify which output fields/files are stochastic vs deterministic before writing tests
- [ ] Deterministic fields in the same file still require exact match
- [ ] For stochastic fields, validate statistical properties: median, confidence interval width, distribution shape
- [ ] Define equivalence criteria upfront (e.g., "95% CI width within 10% of reference", "median within 1 SD")
- [ ] Run multiple iterations to confirm variance is within expected bounds
- [ ] Document which outputs are stochastic and what the acceptance criteria are in test comments

## Common Root Causes of Mismatches
When values don't match, check these frequent sources of divergence:

- [ ] **Coordinate system off-by-one**: 1-based GTF vs 0-based BAM vs closed-interval trees — a single ±1 error propagates to every count
- [ ] **PE/SE detection flags**: `proper_pair` vs `paired` vs CLI-specified `--paired` vs auto-detected from BAM header — these are not equivalent
- [ ] **Secondary vs supplementary alignment filtering**: tools differ on which SAM flag bits they use to exclude reads (0x100 vs 0x800 vs both)
- [ ] **Chromosome naming**: `chr1` vs `1`, `chrM` vs `MT` vs `M` — mismatches cause silent zero counts, not errors
- [ ] **Union vs intersection strategy**: paired-end gene assignment can use the union of both mates' overlaps or the intersection — different tools default differently
- [ ] **Fragment vs read counting**: per-read counts double-count PE fragments; per-fragment counts require mate deduplication
- [ ] **Genomic span vs exonic coordinates**: coverage/length computed on full genomic span (introns included) vs exonic bases only produces very different values
- [ ] **Upstream tool bugs that must be reproduced**: if the upstream tool has a known bug that affects output, the rewrite must reproduce the buggy behavior for compatibility — do not "fix" it
