# tabix/tabix + tabix/bgziptabix → htslib/bgziptabix migration

**Branch:** `topic/cnvkit` (PR #2243) — folded into the existing cnvkit/versions-topic cleanup work in this branch, not a separate PR.

## Background

`nf-core/modules` PR [#11571](https://github.com/nf-core/modules/commit/54e41f4ed) (merged 2026-05-22, "Add HTSLIB/BGZIPTABIX and deprecate redundant modules") deprecated both `tabix/tabix` and `tabix/bgziptabix`. Both modules' upstream `main.nf` now unconditionally `assert false` with a deprecation message pointing at `htslib/bgziptabix` — not a subtle channel-shape break, a hard crash on the next `nf-core modules update`. Sarek is currently pinned to old, working shas for both (`tabix/tabix` @ `f2cfcf9d3`, `tabix/bgziptabix` @ `23004c9c6`), so nothing is broken today, but this is a landmine that any future full-modules update would trigger.

This was discovered while investigating why `prepare_genome`'s `versions` channel still had live `.mix()` calls after the sentieon/tnscope/consensus/annotate dead-versions-plumbing cleanup earlier in this session — `prepare_genome`'s only remaining version producers are 5 aliased `TABIX_TABIX` calls, which are the only thing in that subworkflow still on the classic (non-topic) versions pattern.

## Scope

15 total call sites across 11 files:

- `tabix/tabix`: 9 declared aliases, but one (`TABIX_GERMLINE_VCFS_CONCAT_SORT` in `vcf_concatenate_germline`) is a dead include — never actually called in that file. **8 real call sites.**
- `tabix/bgziptabix`: 6 aliases, one of which (`vcf_annotate_snpeff`) is a vendored `subworkflows/nf-core` file that nf-core/modules has _already_ fixed upstream. **5 real local call sites** to hand-rewire + 1 vendored subworkflow to pull via `nf-core subworkflows update`.

Total real work: 13 mechanical call-site rewires + 1 dead-include deletion + 1 subworkflow update.

## Module setup

1. `nf-core modules install htslib/bgziptabix` — standard tooling, registers `main.nf`/`meta.yml`/`environment.yml` + `modules.json` entry.
2. `nf-core subworkflows update vcf_annotate_snpeff` — pulls upstream's already-fixed version. Upstream's fixed file is confirmed identical to sarek's vendored copy except for the tabix→htslib swap:
   ```groovy
   include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix'
   ...
   HTSLIB_BGZIPTABIX(
       SNPEFF_SNPEFF.out.vcf.map { meta, vcf -> [ meta, vcf, [], [] ] },
       "compress",
       true,
       "vcf"
   )
   ch_vcf_tbi = HTSLIB_BGZIPTABIX.out.output.join(HTSLIB_BGZIPTABIX.out.index)
   ```
   This is the canonical call pattern replicated at the other 13 sites.
3. Once all 13 local sites are rewired and the dead include deleted: `nf-core modules remove tabix/tabix` and `nf-core modules remove tabix/bgziptabix`. Both are fully removed — nothing should reference upstream-dead modules going forward.

## `htslib/bgziptabix` interface

```groovy
input:
tuple val(meta), path(infile), path(infile_tbi), path(regions)
val action       // "compress" | "decompress"
val make_index
val out_ext

output:
tuple val(meta), path("${outfile}"), emit: output
tuple val(meta), path("${outfile}.{tbi,csi}"), emit: index, optional: true
tuple val("${task.process}"), val('htslib'), eval("bgzip --version | ..."), topic: versions, emit: versions_htslib
tuple val("${task.process}"), val('xz'), eval("xz --version | ..."), topic: versions, emit: versions_xz
```

`outfile = action == "compress" ? (out_ext ? "${prefix}.${out_ext}.gz" : "${prefix}.gz") : ...`, `prefix = task.ext.prefix ?: "${meta.id}"`. If the input is already BGZF-compressed and `action == "compress"`, the script symlinks input→outfile instead of re-compressing, then indexes if `make_index`.

Every sarek call site uses `action = 'compress'`, `make_index = true`, `infile_tbi = []`, `regions = []` — none currently pass a pre-existing index or region filter.

**Key semantic gap:** the old `tabix/tabix` module ignores `task.ext.prefix` entirely — its script just runs `tabix $tab`, producing `<original-filename>.tbi` unconditionally. `htslib/bgziptabix` names output from `task.ext.prefix ?: meta.id` instead. So the 8 former-`tabix/tabix` sites need `ext.prefix` **added** (not carried over) to preserve today's exact filenames — this wasn't a config concern before because the old module never read it.

## Call-site mapping

Uniform pattern for all 13 sites: input tuple becomes `[meta, file, [], []]`, call gets `'compress', true, <out_ext>` appended, `.out.tbi` accessor becomes `.out.index`, `.out.gz_index` (old 3-tuple) becomes `.out.output.join(.out.index)` (reconstructing the 3-tuple from two 2-tuples).

| File                                                   | Alias(es)                                                                                                                     | out_ext | ext.prefix change                                                                                                                                                              | Versions cascade                                                                                                                                                                                                                |
| ------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `subworkflows/local/bam_variant_calling_freebayes`     | `TABIX_VC_FREEBAYES`, `TABIX_VC_FREEBAYES_FILT`                                                                               | `''`    | add `{ infile.baseName }` (`conf/modules/freebayes.config`, both blocks)                                                                                                       | subworkflow's `versions` becomes fully dead (was 100% tabix-fed) → full teardown; cascades to `bam_variant_calling_germline_all`/`somatic_all`/`tumor_only_all` (each loses one `.mix()` line, stays alive via other producers) |
| `subworkflows/local/vcf_concatenate_germline`          | delete dead `TABIX_GERMLINE_VCFS_CONCAT_SORT` include; rewire `TABIX_EXT_VCF`                                                 | `''`    | rename `input.baseName`→`infile.baseName` in the single shared `conf/modules/post_variant_calling.config` block (also covers `vcf_normalization`'s instance of the same alias) | no change — this subworkflow only ever mixed `ADD_INFO_TO_VCF.out.versions`, never `TABIX_EXT_VCF`'s                                                                                                                            |
| `subworkflows/local/vcf_varlociraptor_somatic`         | `TABIX_GERMLINE`, `TABIX_SOMATIC`                                                                                             | `''`    | add `{ infile.baseName }` (`conf/modules/varlociraptor.config`, shared block)                                                                                                  | subworkflow's `ch_versions` becomes fully dead (was 100% tabix-fed) → full teardown; cascades to `post_variantcalling` (loses one `.mix()` line, stays alive via `NORMALIZE_VCFS`/`CONCATENATE_GERMLINE_VCFS`)                  |
| `subworkflows/local/prepare_genome`                    | `TABIX_BCFTOOLS_ANNOTATIONS`, `TABIX_DBSNP`, `TABIX_GERMLINE_RESOURCE`, `TABIX_KNOWN_INDELS`, `TABIX_KNOWN_SNPS`, `TABIX_PON` | `''`    | add `{ infile.baseName }` to all 6 blocks in `conf/modules/prepare_genome.config`                                                                                              | subworkflow's `versions` becomes fully dead (was 100% tabix-fed — the exact channel the sentieon investigation was about) → full teardown; cascades to `main.nf`                                                                |
| `subworkflows/local/prepare_intervals`                 | `TABIX_BGZIPTABIX_INTERVAL_SPLIT`, `TABIX_BGZIPTABIX_INTERVAL_COMBINED`                                                       | `'bed'` | keep existing `{"${meta.id}"}` in `conf/modules/prepare_intervals.config` unchanged                                                                                            | no change (was never mixed in this subworkflow)                                                                                                                                                                                 |
| `subworkflows/local/vcf_normalization`                 | `TABIX_EXT_VCF` (shares config with `vcf_concatenate_germline`'s instance)                                                    | `''`    | (see above, shared block)                                                                                                                                                      | no change                                                                                                                                                                                                                       |
| `subworkflows/local/bam_variant_calling_single_tiddit` | `TABIX_BGZIP_TIDDIT_SV`                                                                                                       | `'vcf'` | keep existing meta.id-based prefixes in `conf/modules/tiddit.config` (base + normal/tumor overrides) unchanged                                                                 | no change (this file has no versions tracking at all)                                                                                                                                                                           |
| `subworkflows/nf-core/vcf_annotate_snpeff`             | via `nf-core subworkflows update`                                                                                             | —       | —                                                                                                                                                                              | already fully dead from the earlier cleanup pass, no further change                                                                                                                                                             |

Also: `conf/modules/annotate.config`'s wildcard `withName: '...(TABIX_BGZIPTABIX|TABIX_TABIX)'` with `ext.prefix = { input.name - '.vcf' }` needs `input`→`infile` renamed; the `TABIX_TABIX` alternative in that regex was already dead (nothing under `VCF_ANNOTATE_ALL` ever used it) and can be dropped.

## Testing & verification

- **Static:** `nextflow lint` after each edit (hook-enforced already), full `pre-commit run --all-files` once all 13 sites + deletion + subworkflow update land. Lint catches syntax/arity but not the tuple-shape change (`gz_index` 3-tuple → `output`+`index` 2-tuples), which is the primary correctness risk here — more so than most of this migration's prior module bumps.
- **Blast radius:** unlike cnvkit's own migration, these sites span `freebayes`, `varlociraptor`, `prepare_genome`/`prepare_intervals` (exercised by nearly every test profile), `tiddit`, and germline VCF concatenation. This needs the full nf-test matrix on CI, not a narrow cnvkit-scoped subset.
- **Expected snapshot fallout:**
  - Version-dict entries change shape: `htslib/bgziptabix` emits topic entries keyed `htslib` (bgzip version) and `xz` (xz version), replacing whatever `tabix`/`bgzip` keys existed for these process names in every `.snap` exercising them — same mechanical find/replace pattern as the `CNVKIT_BATCH` → `samtools` entry added earlier this session.
  - Container changes (`community.wave.seqera.io/library/htslib_xz:...` vs old `htslib:1.21`) — worth checking whether the effective htslib version itself bumps.
  - File content (md5s) should be unaffected — same underlying `bgzip`/`tabix` binaries — _provided_ every site's `ext.prefix` is correct. A wrong prefix produces a wrong/missing output filename (a real bug), not just a snapshot diff — that's the top regression risk to watch for in CI logs specifically, not just diff away.
  - Snapshot regeneration goes through the CI-based `regenerate-nf-test-snapshots` mechanism (ARM Mac / sandbox can't run these containers locally), same as the rest of this migration.
- **Order:** rewire all 13 sites + subworkflow update in the working tree → lint clean → pre-commit clean → push → full CI run → pull real md5/version deltas from CI → apply to snapshots → re-verify.

## Out of scope

- Any further nf-core/modules deprecations beyond `tabix/tabix`/`tabix/bgziptabix` (e.g. `tabix/bgzip`, `samtools/bgzip`, also deprecated in the same upstream PR but not used anywhere in sarek).
- `modules/local/*` migrations (add_info_to_vcf, create_intervals_bed, samtools/reindex_bam) — that's `topic/local` PR #2244's scope.
- `bam_variant_calling_indexcov` and `utils_nfcore_sarek_pipeline` — confirmed not dead, or pre-existing template boilerplate respectively; left untouched per the earlier cleanup pass.
