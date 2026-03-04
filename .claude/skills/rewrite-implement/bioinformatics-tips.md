# Bioinformatics Implementation Tips

Lessons learned from rewriting bioinformatics tools (dupRadar, featureCounts, RSeQC,
Qualimap, preseq, samtools) in Rust. These are domain-specific gotchas that will save
you hours of debugging.

## BAM File Processing

### Coordinate Systems
- BAM uses **0-based, half-open** coordinates: `[start, end)`
- GTF uses **1-based, inclusive** coordinates: `[start, end]`
- Picard/htsjdk uses **1-based, inclusive**: `[start, end]`
- BED uses **0-based, half-open**: `[start, end)`
- Convert GTF to 0-based half-open: `(gtf_start - 1, gtf_end)`
- ALWAYS document which coordinate system each variable uses

### CIGAR Operations
- M (match/mismatch): consumes both query and reference
- I (insertion): consumes query only
- D (deletion): consumes reference only
- N (reference skip / intron): consumes reference only
- S (soft clip): consumes query only
- H (hard clip): consumes neither (but length is in CIGAR)
- Some upstream tools have CIGAR bugs — e.g., advancing reference position for
  insertion (I) or soft-clip (S) operations. You MUST replicate these bugs.

### PE Read Handling
- Buffer mates by read name for paired-end reconciliation
- Handle edge cases: orphaned mates (mate filtered/unmapped), cross-chromosome pairs,
  supplementary alignments sharing read names
- `record.is_proper_pair()` checks the SAM flag, not actual pair validity
- NH tag indicates multi-mapping: NH > 1 means multi-mapped
- Secondary alignments (flag 0x100) typically have NH > 1 but check your data
- preseq-style fragment counting uses `proper_pair` flag (0x2) not `paired` flag (0x1)
  for PE/SE mode detection — tools vary in which flags they check
- When computing fragment start from TLEN for reverse-strand mates, use
  `cigar_end + tlen`, NOT `pos + tlen`. This is a common subtle bug.
- Tools may count both mates (2 reads per fragment) or count fragments (1 per pair).
  Verify counting semantics from the upstream source.
- For multi-mapped PE reads, mate buffer keys must include position to disambiguate
  alignment locations: `(read_name, (min(pos, mpos), max(pos, mpos)))`
- Orphaned proper-pair mates (whose partner was filtered) should typically be DISCARDED,
  not processed as single-end reads — verify against upstream behavior
- UNION vs INTERSECTION for PE gene hits can differ significantly: if R1 overlaps
  genes {A, B} and R2 overlaps {A}, intersection assigns to A while union marks as
  ambiguous

### Multi-mapper Handling
- Check NH:i tag for multi-mapping status
- Default behavior varies by tool: some include, some exclude multi-mappers
- Some tools use NH tag, others use MAPQ (MAPQ=0 in some aligners means multi-mapped)
- Document exactly which reads your rewrite includes/excludes

### Single-Pass Architecture
- The #1 performance win for multi-tool rewrites: iterate the BAM file ONCE
- Collect all analyses in parallel during the same pass
- Buffer PE mates efficiently (HashMap by read name)
- Flush buffers at chromosome boundaries for sorted BAMs

## GTF/GFF Annotation Parsing

### Attribute Extraction
- GTF attributes: `gene_id "ENSG00000..."; gene_name "TP53";`
- GFF3 attributes: `ID=gene:ENSG00000...;Name=TP53`
- GENCODE GTFs use `gene_type` for biotype
- Ensembl GTFs use `gene_biotype` for biotype
- Auto-detect: check which attribute exists and fall back gracefully

### Transcript Structure
- A gene has multiple transcripts, each with its own set of exons
- Exons within a transcript may overlap (rare but possible in some annotations)
- Merged gene model: union of all exons across all transcripts, then merge overlapping
- Per-transcript model: each transcript keeps its own exon list
- Different tools use different models — read the upstream source to know which

### Interval Operations
- **Enclosure**: interval A encloses B if `A.start <= B.start && A.end >= B.end`
- **Overlap**: A overlaps B if `A.start < B.end && B.start < A.end` (half-open)
- **Abutting**: A abuts B if `A.end == B.start` (half-open) — whether to merge these
  depends on the upstream tool. Some do, some don't.
- Use `coitrees` crate for cache-oblivious interval trees (fast queries)
- COITree uses **closed** interval semantics `[start, last]` internally. Passing 0-based
  half-open intervals directly extends each interval by 1 base, causing false-positive
  overlaps. Convert to closed: `(start, end - 1)`.
- Coverage arrays for gene body coverage should use EXONIC coordinates (sum of exon
  lengths), not genomic span. Intronic zeros dilute mean coverage.

### Gzip Transparency
- GTF/GFF/BED files may be gzip-compressed (`.gz` extension)
- Detect by magic bytes (`1f 8b`), not file extension
- Use `flate2` crate for transparent decompression

## Float Formatting and Output Matching

### R Output Conventions
- R prints up to 15 significant digits by default
- R uses `NA` (not `NaN` or `nan`) for missing values
- R's `format()` and `sprintf()` have specific rounding behavior
- Trailing zeros: R may or may not include them depending on context
- Scientific notation: R uses `e+00` format (e.g., `1.23e+05`)

### Matching Float Output
- Use `format!("{:.15e}", value)` as a starting point, then adjust
- Write a formatting function that matches the upstream tool's exact output
- Test with edge cases: 0.0, very small numbers, very large numbers, negative numbers
- NaN handling: check what the upstream tool outputs for NaN (could be "NA", "NaN", "nan", "N/A")
- Infinity: check upstream behavior
- Tools vary widely in float formatting: samtools stats uses `{:.6e}` for error rates,
  1-decimal for percentages; other tools use comma-separated thousands
- Float formatting alone can consume a full implementation session. Budget time for it.

### Numeric Comparison
- For validation, use relative tolerance: `|a - b| / max(|a|, |b|) < 1e-6`
- For exact string matching: compare formatted strings directly
- Watch for: sort instability affecting output order of equal-valued rows

## Common Upstream Bugs to Replicate

These are real bugs found in widely-used tools that you MUST replicate to match output:

1. **CIGAR offset bug**: Some tools advance the reference position for ALL CIGAR
   operations (including I, S, H, P which should NOT advance reference). This shifts
   M-block positions for reads with insertions or soft clips before the first M block.

2. **Off-by-one in coverage**: Picard's `addCoverageCounts(start, end)` treats `end`
   as exclusive, but Qualimap passes 1-based inclusive end coordinates — missing the
   last base of each alignment block.

3. **TreeMap duplicate key loss**: Java's TreeMap silently drops entries with duplicate
   keys. If keyed by a floating-point mean coverage, transcripts with identical means
   are silently lost.

4. **Junction counting outside guard**: Some tools increment a junction counter for
   every N CIGAR operation, even when the junction motif extraction fails (e.g., near
   read boundaries where there aren't enough bases to extract the motif).

5. **Sort instability**: When sorting by a floating-point value with ties, different
   implementations may produce different orders. Match the upstream tool's sort behavior.

## Source Code Analysis

- If the upstream source is only available as compiled bytecode (Java `.class`,
  Python `.pyc`), decompiling is acceptable as a last resort. Use CFR for Java,
  uncompyle6 for Python.
- Original tools often have actual bugs (e.g., CIGAR parsing errors) that must be
  REPRODUCED to match output. Document with `// BUG COMPAT:` comments.
- Check the target pipeline's actual invocation flags in the module's `main.nf` —
  flags like `-pe`, strandedness mode, and annotation parameters change output
  significantly.

## Testing and Validation

- Different RNG implementations (ChaCha8 vs Mersenne Twister vs mt19937) produce
  different sequences. For stochastic outputs, validate statistical properties
  (median, CI width) not exact values.
- Hash-based counting approaches (HashMap of positions) are often more correct than
  sort-dependent sequential counting, but will produce different output ordering.
  Document the improvement.

## Performance Patterns

### Memory
- Pre-allocate HashMap/Vec capacity when size is known
- Use `&[u8]` instead of `String` for sequence data in hot loops
- For large interval trees, build once and query many times
- Flush PE mate buffers at chromosome boundaries

### Parallelism
- `rayon` for data-parallel iteration (per-chromosome, per-transcript)
- For BAM processing, use `rust-htslib`'s thread pool for decompression
- Be careful with shared mutable state — prefer per-thread accumulators merged at end

### I/O
- Buffer output writes (`BufWriter`)
- For BAM: single pass is 10-100x faster than multiple passes
- Read annotation files once at startup, build indexes, then process BAM
