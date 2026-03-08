# Rewrite Target Analysis

## Pipeline Summary

- **Pipeline**: `nf-core/sarek` (commit `3880e26b` on `evanfloden/sarek`)
- **Repository**: `https://github.com/nf-core/sarek`
- **Data source**: Seqera Platform run metadata (`run_id=4rk1M6nCGI6cYF`) and task export (`/tmp/seqera_tasks_4rk1M6nCGI6cYF.json`)
- **Run name**: `sarek-germline-opt-v14-trial-0`
- **Run status**: `SUCCEEDED`
- **Total pipeline wall time**: ~1.5 hours
- **Total tasks**: 117
- **Total task-minutes**: 378.9
- **Environment**: AWS Batch (`r5ad` family), `eu-west-1`
- **Config context**: dragmap aligner, mpileup-only calling, GRCh38, 1 germline WGS sample split into 21 intervals

### Scoring Method Used

- Difficulty mapping: `low=1`, `medium=3`, `high=9`
- ROI formula: `(predicted_savings_per_task_minutes x task_count) / difficulty_score`
- Ranking below uses measured runtime from this run only
- A separate `N=20` projection (10 tumor-normal pairs) is included and clearly labeled as projection

## Tool Inventory

### Processes Executed In This Run (real runtime data)

| Tool / Process | Language | Per-sample? | Mean runtime (min) | Total runtime (min) | Source URL | Notes |
|---|---|---:|---:|---:|---|---|
| `DRAGMAP_ALIGN` | C++ | Yes | 11.4 | 136.7 | `https://github.com/Illumina/dragmap` | Dominant wall time; already highly optimized native aligner |
| `BCFTOOLS_MPILEUP` | C | Yes (interval-scattered) | 4.6 | 96.2 | `https://github.com/samtools/bcftools` | Reads 163.97 GB; I/O-heavy |
| `GATK4_BASERECALIBRATOR` | Java | Yes (interval-scattered) | 4.2 | 87.9 | `https://github.com/broadinstitute/gatk` | BQSR modeling stage |
| `GATK4_APPLYBQSR` | Java | Yes (interval-scattered) | 1.5 | 31.8 | `https://github.com/broadinstitute/gatk` | BQSR application stage |
| `MOSDEPTH` | Rust | Yes | 3.0 | 5.9 | `https://github.com/brentp/mosdepth` | Already replaced by `ironqc` option in this repo |
| `GATK4_MARKDUPLICATES` | Java | Yes | 5.0 | 5.0 | `https://github.com/broadinstitute/gatk` | Single task in this run |
| `SAMTOOLS_STATS` | C | Yes | 2.4 | 4.9 | `https://github.com/samtools/samtools` | Already replaced by `ironqc` option in this repo |
| `FASTQC` | Java | Yes | 2.0 | 2.0 | `https://www.bioinformatics.babraham.ac.uk/projects/fastqc/` | Small runtime in this config |
| `FASTP` | C++ | Yes | 0.6 | 0.6 | `https://github.com/OpenGene/fastp` | Small runtime in this config |
| `MULTIQC` | Python | Global | 0.7 | 0.7 | `https://github.com/MultiQC/MultiQC` | Aggregation only |

### Not Exercised In This Run (important caveat)

This run did not execute tumor-normal, SV, CNV, annotation, or multi-caller paths (`Mutect2`, `HaplotypeCaller`, `DeepVariant`, `Manta`, `Strelka`, `SnpEff`, `VEP`, `CNVkit`, `Control-FREEC`, `ASCAT`, etc.). These remain in scope for future ranking once measured runtime is available from those modes.

## Completed Rewrites

### `ironqc` (COMPLETED)

- **Status**: Completed and integrated via `params.use_ironqc` in `subworkflows/local/cram_qc_mosdepth_samtools/main.nf`
- **Bundle delivered**: `SAMTOOLS_STATS` + `MOSDEPTH` + `GOLEFT_INDEXCOV` in one Rust binary (`/Users/evan/projects/sarek-investigate/ironqc/`)
- **Architecture**: single-pass BAM/CRAM iteration with `noodles` + `rayon` accumulators
- **Current run impact (measured baseline replaced)**: `SAMTOOLS_STATS (4.9 min)` + `MOSDEPTH (5.9 min)` = **10.8 min**
- **Share of total run**: `10.8 / 378.9 = 2.85%`
- **N=20 projection**: ~216 task-minutes removed from baseline path (before accounting for second-order staging/container overhead)

## Ranked Rewrite Targets

### Target 1: GATK BQSR Chain Fusion (`BASERECALIBRATOR` + `APPLYBQSR`) (Priority: HIGH)

- **Tools to bundle**: `GATK4_BASERECALIBRATOR`, `GATK4_APPLYBQSR` (optionally keep `GATK4_GATHERBQSRREPORTS` as compatibility shim)
- **Current total runtime**: `119.7 min` (`87.9 + 31.8`), plus `2.4 min` gather
- **Predicted savings**: `~47.8 min` total (`~40%` on the 119.7 min chain; conservative for high-complexity Java statistical model)
- **Difficulty**: `high` (score `9`)
- **Estimated token cost**: `$900-2000`
- **Shared inputs**: recalibration BAM/CRAM, known-sites VCF, reference
- **Source code**: `https://github.com/broadinstitute/gatk`
- **Risk factors**: covariate table parity, floating-point reproducibility, read-group edge cases, strict report-format compatibility
- **Bundling rationale**: same data path and tightly coupled stages; fusing can remove repeated serialization and task boundaries
- **ROI score**: per-task savings `1.90 min` x `21` tasks / `9` = **4.43**
- **N=20 projection**: `~956 min` saved (same % reduction applied to measured single-sample baseline)

### Target 2: BCFTOOLS mpileup Modernization (`MPILEUP` + merge shim) (Priority: MEDIUM)

- **Tools to bundle**: `BCFTOOLS_MPILEUP`, `MERGE_BCFTOOLS_MPILEUP`, and optional trivial QC emitters
- **Current total runtime**: `96.4 min` (`96.2 + 0.2`)
- **Predicted savings**: `~14.4 min` total (`~15%`, realistic upper bound for C-to-Rust rewrite on already optimized htslib stack)
- **Difficulty**: `high` (score `9`)
- **Estimated token cost**: `$700-1600`
- **Shared inputs**: interval BAM/CRAM + reference + mpileup params
- **Source code**: `https://github.com/samtools/bcftools`, `https://github.com/samtools/htslib`
- **Risk factors**: mpileup edge-case parity, BAQ behavior, pileup filtering semantics, deep compatibility burden
- **Bundling rationale**: little arithmetic speed gain expected; value is mostly reduced process handoff and minor I/O pipeline cleanup
- **ROI score**: per-task savings `0.69 min` x `21` tasks / `9` = **1.61**
- **N=20 projection**: `~288 min` saved

### Target 3: MarkDuplicates + CRAM Finalization Chain (Priority: LOW-MEDIUM)

- **Tools to bundle**: `GATK4_MARKDUPLICATES`, `MERGE_CRAM`, `INDEX_CRAM`
- **Current total runtime**: `6.6 min` (`5.0 + 1.5 + 0.1`)
- **Predicted savings**: `~2.0 min` total (`~30%`, conservative)
- **Difficulty**: `high` (score `9`)
- **Estimated token cost**: `$600-1400`
- **Shared inputs**: post-alignment BAM/CRAM stream
- **Source code**: `https://github.com/broadinstitute/gatk`, `https://github.com/samtools/samtools`
- **Risk factors**: duplicate-marking equivalence and metrics parity are non-trivial; modest upside in this run configuration
- **Bundling rationale**: mostly a DAG simplification play, not an algorithmic speedup play
- **ROI score**: per-task savings `2.00 min` x `1` task / `9` = **0.22**
- **N=20 projection**: `~40 min` saved

### Target 4: FASTQ QC/Trim Consolidation (`FASTP` + `FASTQC`) (Priority: LOW)

- **Tools to bundle**: `FASTP`, `FASTQC`
- **Current total runtime**: `2.6 min`
- **Predicted savings**: `~1.0 min` total (`~38%`)
- **Difficulty**: `medium` (score `3`)
- **Estimated token cost**: `$250-700`
- **Shared inputs**: FASTQ reads
- **Source code**: `https://github.com/OpenGene/fastp`, `https://www.bioinformatics.babraham.ac.uk/projects/fastqc/`
- **Risk factors**: plot/report compatibility overhead dominates implementation effort
- **Bundling rationale**: can remove one pass and one task, but impact is small in current profile
- **ROI score**: per-task savings `1.00 min` x `1` task / `3` = **0.33**
- **N=20 projection**: `~20 min` saved

### Explicitly Not Recommended (despite high runtime)

- **`DRAGMAP_ALIGN` rewrite**: huge runtime share, but C++ aligners are already near hardware limits; rewrite difficulty is extreme and expected gains are marginal.
- **Already-Rust replacements**: do not target `varlociraptor` family or completed `ironqc` bundle.

## Bundling Groups

| Group | Tools | Shared input | Single-pass feasible? |
|---|---|---|---|
| `BAM-QC-COMPLETED` | `SAMTOOLS_STATS`, `MOSDEPTH`, `GOLEFT_INDEXCOV` via `IRONQC` | BAM/CRAM + index (+ optional BED) | Yes (already implemented) |
| `BQSR-CHAIN` | `GATK4_BASERECALIBRATOR`, `GATK4_APPLYBQSR` | recalibration BAM/CRAM + known sites | Partial (same read traversal logic, but statistical/state model is complex) |
| `MPILEUP-CHAIN` | `BCFTOOLS_MPILEUP`, `MERGE_BCFTOOLS_MPILEUP` | interval BAM/CRAM + reference | Partial (shared code path, limited upside) |
| `DEDUP-FINALIZE` | `GATK4_MARKDUPLICATES`, `MERGE_CRAM`, `INDEX_CRAM` | post-alignment BAM/CRAM | Partial |
| `FASTQ-QC` | `FASTP`, `FASTQC` | FASTQ | Yes |

## Recommended Next Target

`GATK BQSR Chain Fusion` - it is the best remaining measured ROI after `ironqc`, and unlike aligner/mpileup rewrites it has meaningful JVM overhead and workflow-boundary overhead that can plausibly be reduced despite high algorithmic complexity.

## Data Limitations And Ranking Sensitivity

- This dataset is **single-sample germline**; per-sample tasks only ran once.
- Only **mpileup** calling was enabled; no HaplotypeCaller/Mutect2/DeepVariant/Manta/Strelka.
- No annotation stack executed (`SnpEff`, `VEP`, `SnpSift`).
- No CNV stack executed (`CNVkit`, `Control-FREEC`, `ASCAT`).
- No SV stack executed (`Manta`, `Strelka`, `TIDDIT`, `SVDB`).
- In full tumor-normal multi-caller runs, VCF post-processing and annotation bundles would likely move up significantly in absolute savings; this report does not assign them measured ROI without runtime evidence.

## Methodology Notes

- Runtime numbers are from the provided Seqera run summary and consistent with `/tmp/seqera_tasks_4rk1M6nCGI6cYF.json` task metadata.
- Rankings use measured per-process totals from this run and conservative speedup assumptions by language class:
  - C/C++ to Rust rewrites: generally `10-30%` (often closer to lower end for mature htslib-class tools)
  - Java to Rust rewrites: potentially `2-5x` in favorable cases, but capped conservatively here due to statistical complexity and parity risk
- Difficulty reflects implementation and parity burden, not just code volume.
- Cost ranges include realism caveats: complex numerical parity, float/text formatting parity, and plot reproduction often expand effort by `3-5x` from optimistic first estimates.

## Source Availability And Hosting Risk (preserved reference)

- **Open-source, straightforward access (GitHub or project docs)**: bwa/bwamem2, samtools/bcftools/tabix, mosdepth, goleft, fastp, fastqc, cnvkit, freebayes, manta, strelka, tiddit, deepvariant, snpeff/snpsift, VEP, vcflib, vcftools, varlociraptor, rust-bio-tools.
- **Closed/proprietary binaries (high rewrite risk due unavailable source or licensing constraints)**: Sentieon (`sentieon/*`), NVIDIA Parabricks (`parabricks/fq2bam`).
- **Mixed script ecosystems (higher parity burden due plotting/format details)**: ASCAT (R), Control-FREEC ancillary scripts, MuSE, NGSCheckMate.

## MultiQC Compatibility (preserved reference)

From `assets/multiqc_config.yml`, primary modules consumed are: `fastqc`, `fastp`, `bbmap`, `picard`, `samtools`, `mosdepth`, `gatk`, `bcftools`, `vcftools`, `snpeff`, `vep`, plus custom `dedup_metrics` table.

Implication: rewrites must preserve report file names and headers for these producers:

- `FASTQC`, `FASTP`, `BBMAP_BBSPLIT`
- `GATK4_MARKDUPLICATES` or `GATK4_ESTIMATELIBRARYCOMPLEXITY`
- `SAMTOOLS_STATS`, `MOSDEPTH`
- `BCFTOOLS_STATS`, `VCFTOOLS`
- `SNPEFF_SNPEFF`, `ENSEMBLVEP_VEP`
- `SENTIEON_DEDUP` custom metrics TSV

## Appendix: Exhaustive Process-to-Module Mapping (preserved)

| Process | Module path |
|---|---|
| `ADD_INFO_TO_VCF` | `modules/local/add_info_to_vcf/main.nf` |
| `CONSENSUS_FROM_SITES` | `modules/local/consensus_from_sites/main.nf` |
| `CREATE_INTERVALS_BED` | `modules/local/create_intervals_bed/main.nf` |
| `SAMTOOLS_REINDEX_BAM` | `modules/local/samtools/reindex_bam/main.nf` |
| `ASCAT` | `modules/nf-core/ascat/main.nf` |
| `BBMAP_BBSPLIT` | `modules/nf-core/bbmap/bbsplit/main.nf` |
| `BCFTOOLS_ANNOTATE` | `modules/nf-core/bcftools/annotate/main.nf` |
| `BCFTOOLS_CONCAT` | `modules/nf-core/bcftools/concat/main.nf` |
| `BCFTOOLS_ISEC` | `modules/nf-core/bcftools/isec/main.nf` |
| `BCFTOOLS_MERGE` | `modules/nf-core/bcftools/merge/main.nf` |
| `BCFTOOLS_MPILEUP` | `modules/nf-core/bcftools/mpileup/main.nf` |
| `BCFTOOLS_NORM` | `modules/nf-core/bcftools/norm/main.nf` |
| `BCFTOOLS_SORT` | `modules/nf-core/bcftools/sort/main.nf` |
| `BCFTOOLS_STATS` | `modules/nf-core/bcftools/stats/main.nf` |
| `BCFTOOLS_VIEW` | `modules/nf-core/bcftools/view/main.nf` |
| `BWA_INDEX` | `modules/nf-core/bwa/index/main.nf` |
| `BWA_MEM` | `modules/nf-core/bwa/mem/main.nf` |
| `BWAMEM2_INDEX` | `modules/nf-core/bwamem2/index/main.nf` |
| `BWAMEM2_MEM` | `modules/nf-core/bwamem2/mem/main.nf` |
| `CAT_CAT` | `modules/nf-core/cat/cat/main.nf` |
| `CAT_FASTQ` | `modules/nf-core/cat/fastq/main.nf` |
| `CNVKIT_ANTITARGET` | `modules/nf-core/cnvkit/antitarget/main.nf` |
| `CNVKIT_BATCH` | `modules/nf-core/cnvkit/batch/main.nf` |
| `CNVKIT_CALL` | `modules/nf-core/cnvkit/call/main.nf` |
| `CNVKIT_EXPORT` | `modules/nf-core/cnvkit/export/main.nf` |
| `CNVKIT_GENEMETRICS` | `modules/nf-core/cnvkit/genemetrics/main.nf` |
| `CNVKIT_REFERENCE` | `modules/nf-core/cnvkit/reference/main.nf` |
| `CONTROLFREEC_ASSESSSIGNIFICANCE` | `modules/nf-core/controlfreec/assesssignificance/main.nf` |
| `CONTROLFREEC_FREEC` | `modules/nf-core/controlfreec/freec/main.nf` |
| `CONTROLFREEC_FREEC2BED` | `modules/nf-core/controlfreec/freec2bed/main.nf` |
| `CONTROLFREEC_FREEC2CIRCOS` | `modules/nf-core/controlfreec/freec2circos/main.nf` |
| `CONTROLFREEC_MAKEGRAPH2` | `modules/nf-core/controlfreec/makegraph2/main.nf` |
| `DEEPVARIANT_RUNDEEPVARIANT` | `modules/nf-core/deepvariant/rundeepvariant/main.nf` |
| `DRAGMAP_ALIGN` | `modules/nf-core/dragmap/align/main.nf` |
| `DRAGMAP_HASHTABLE` | `modules/nf-core/dragmap/hashtable/main.nf` |
| `ENSEMBLVEP_DOWNLOAD` | `modules/nf-core/ensemblvep/download/main.nf` |
| `ENSEMBLVEP_VEP` | `modules/nf-core/ensemblvep/vep/main.nf` |
| `FASTP` | `modules/nf-core/fastp/main.nf` |
| `FASTQC` | `modules/nf-core/fastqc/main.nf` |
| `FGBIO_CALLMOLECULARCONSENSUSREADS` | `modules/nf-core/fgbio/callmolecularconsensusreads/main.nf` |
| `FGBIO_COPYUMIFROMREADNAME` | `modules/nf-core/fgbio/copyumifromreadname/main.nf` |
| `FGBIO_FASTQTOBAM` | `modules/nf-core/fgbio/fastqtobam/main.nf` |
| `FGBIO_GROUPREADSBYUMI` | `modules/nf-core/fgbio/groupreadsbyumi/main.nf` |
| `FREEBAYES` | `modules/nf-core/freebayes/main.nf` |
| `GATK4_APPLYBQSR` | `modules/nf-core/gatk4/applybqsr/main.nf` |
| `GATK4_APPLYVQSR` | `modules/nf-core/gatk4/applyvqsr/main.nf` |
| `GATK4_BASERECALIBRATOR` | `modules/nf-core/gatk4/baserecalibrator/main.nf` |
| `GATK4_CALCULATECONTAMINATION` | `modules/nf-core/gatk4/calculatecontamination/main.nf` |
| `GATK4_CNNSCOREVARIANTS` | `modules/nf-core/gatk4/cnnscorevariants/main.nf` |
| `GATK4_CREATESEQUENCEDICTIONARY` | `modules/nf-core/gatk4/createsequencedictionary/main.nf` |
| `GATK4_ESTIMATELIBRARYCOMPLEXITY` | `modules/nf-core/gatk4/estimatelibrarycomplexity/main.nf` |
| `GATK4_FILTERMUTECTCALLS` | `modules/nf-core/gatk4/filtermutectcalls/main.nf` |
| `GATK4_FILTERVARIANTTRANCHES` | `modules/nf-core/gatk4/filtervarianttranches/main.nf` |
| `GATK4_GATHERBQSRREPORTS` | `modules/nf-core/gatk4/gatherbqsrreports/main.nf` |
| `GATK4_GATHERPILEUPSUMMARIES` | `modules/nf-core/gatk4/gatherpileupsummaries/main.nf` |
| `GATK4_GENOMICSDBIMPORT` | `modules/nf-core/gatk4/genomicsdbimport/main.nf` |
| `GATK4_GENOTYPEGVCFS` | `modules/nf-core/gatk4/genotypegvcfs/main.nf` |
| `GATK4_GETPILEUPSUMMARIES` | `modules/nf-core/gatk4/getpileupsummaries/main.nf` |
| `GATK4_HAPLOTYPECALLER` | `modules/nf-core/gatk4/haplotypecaller/main.nf` |
| `GATK4_INTERVALLISTTOBED` | `modules/nf-core/gatk4/intervallisttobed/main.nf` |
| `GATK4_LEARNREADORIENTATIONMODEL` | `modules/nf-core/gatk4/learnreadorientationmodel/main.nf` |
| `GATK4_MARKDUPLICATES` | `modules/nf-core/gatk4/markduplicates/main.nf` |
| `GATK4_MERGEMUTECTSTATS` | `modules/nf-core/gatk4/mergemutectstats/main.nf` |
| `GATK4_MERGEVCFS` | `modules/nf-core/gatk4/mergevcfs/main.nf` |
| `GATK4_MUTECT2` | `modules/nf-core/gatk4/mutect2/main.nf` |
| `GATK4_VARIANTRECALIBRATOR` | `modules/nf-core/gatk4/variantrecalibrator/main.nf` |
| `GATK4SPARK_APPLYBQSR` | `modules/nf-core/gatk4spark/applybqsr/main.nf` |
| `GATK4SPARK_BASERECALIBRATOR` | `modules/nf-core/gatk4spark/baserecalibrator/main.nf` |
| `GATK4SPARK_MARKDUPLICATES` | `modules/nf-core/gatk4spark/markduplicates/main.nf` |
| `GAWK` | `modules/nf-core/gawk/main.nf` |
| `GOLEFT_INDEXCOV` | `modules/nf-core/goleft/indexcov/main.nf` |
| `GUNZIP` | `modules/nf-core/gunzip/main.nf` |
| `LOFREQ_CALLPARALLEL` | `modules/nf-core/lofreq/callparallel/main.nf` |
| `MANTA_GERMLINE` | `modules/nf-core/manta/germline/main.nf` |
| `MANTA_SOMATIC` | `modules/nf-core/manta/somatic/main.nf` |
| `MANTA_TUMORONLY` | `modules/nf-core/manta/tumoronly/main.nf` |
| `MOSDEPTH` | `modules/nf-core/mosdepth/main.nf` |
| `MSISENSOR2_MSI` | `modules/nf-core/msisensor2/msi/main.nf` |
| `MSISENSORPRO_MSISOMATIC` | `modules/nf-core/msisensorpro/msisomatic/main.nf` |
| `MSISENSORPRO_SCAN` | `modules/nf-core/msisensorpro/scan/main.nf` |
| `MULTIQC` | `modules/nf-core/multiqc/main.nf` |
| `MUSE_CALL` | `modules/nf-core/muse/call/main.nf` |
| `MUSE_SUMP` | `modules/nf-core/muse/sump/main.nf` |
| `NGSCHECKMATE_NCM` | `modules/nf-core/ngscheckmate/ncm/main.nf` |
| `PARABRICKS_FQ2BAM` | `modules/nf-core/parabricks/fq2bam/main.nf` |
| `RBT_VCFSPLIT` | `modules/nf-core/rbt/vcfsplit/main.nf` |
| `SAMTOOLS_BAM2FQ` | `modules/nf-core/samtools/bam2fq/main.nf` |
| `SAMTOOLS_COLLATEFASTQ` | `modules/nf-core/samtools/collatefastq/main.nf` |
| `SAMTOOLS_CONVERT` | `modules/nf-core/samtools/convert/main.nf` |
| `SAMTOOLS_FAIDX` | `modules/nf-core/samtools/faidx/main.nf` |
| `SAMTOOLS_INDEX` | `modules/nf-core/samtools/index/main.nf` |
| `SAMTOOLS_MERGE` | `modules/nf-core/samtools/merge/main.nf` |
| `SAMTOOLS_MPILEUP` | `modules/nf-core/samtools/mpileup/main.nf` |
| `SAMTOOLS_STATS` | `modules/nf-core/samtools/stats/main.nf` |
| `SAMTOOLS_VIEW` | `modules/nf-core/samtools/view/main.nf` |
| `SENTIEON_APPLYVARCAL` | `modules/nf-core/sentieon/applyvarcal/main.nf` |
| `SENTIEON_BWAMEM` | `modules/nf-core/sentieon/bwamem/main.nf` |
| `SENTIEON_DEDUP` | `modules/nf-core/sentieon/dedup/main.nf` |
| `SENTIEON_DNAMODELAPPLY` | `modules/nf-core/sentieon/dnamodelapply/main.nf` |
| `SENTIEON_DNASCOPE` | `modules/nf-core/sentieon/dnascope/main.nf` |
| `SENTIEON_GVCFTYPER` | `modules/nf-core/sentieon/gvcftyper/main.nf` |
| `SENTIEON_HAPLOTYPER` | `modules/nf-core/sentieon/haplotyper/main.nf` |
| `SENTIEON_TNSCOPE` | `modules/nf-core/sentieon/tnscope/main.nf` |
| `SENTIEON_VARCAL` | `modules/nf-core/sentieon/varcal/main.nf` |
| `SNPEFF_DOWNLOAD` | `modules/nf-core/snpeff/download/main.nf` |
| `SNPEFF_SNPEFF` | `modules/nf-core/snpeff/snpeff/main.nf` |
| `SNPSIFT_ANNMEM` | `modules/nf-core/snpsift/annmem/main.nf` |
| `SNPSIFT_ANNMEMCREATE` | `modules/nf-core/snpsift/annmemcreate/main.nf` |
| `SPRING_DECOMPRESS` | `modules/nf-core/spring/decompress/main.nf` |
| `STRELKA_GERMLINE` | `modules/nf-core/strelka/germline/main.nf` |
| `STRELKA_SOMATIC` | `modules/nf-core/strelka/somatic/main.nf` |
| `SVDB_MERGE` | `modules/nf-core/svdb/merge/main.nf` |
| `TABIX_BGZIPTABIX` | `modules/nf-core/tabix/bgziptabix/main.nf` |
| `TABIX_TABIX` | `modules/nf-core/tabix/tabix/main.nf` |
| `TIDDIT_SV` | `modules/nf-core/tiddit/sv/main.nf` |
| `UNTAR` | `modules/nf-core/untar/main.nf` |
| `UNZIP` | `modules/nf-core/unzip/main.nf` |
| `VARLOCIRAPTOR_CALLVARIANTS` | `modules/nf-core/varlociraptor/callvariants/main.nf` |
| `VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES` | `modules/nf-core/varlociraptor/estimatealignmentproperties/main.nf` |
| `VARLOCIRAPTOR_PREPROCESS` | `modules/nf-core/varlociraptor/preprocess/main.nf` |
| `VCFLIB_VCFFILTER` | `modules/nf-core/vcflib/vcffilter/main.nf` |
| `VCFTOOLS` | `modules/nf-core/vcftools/main.nf` |
| `YTE` | `modules/nf-core/yte/main.nf` |
