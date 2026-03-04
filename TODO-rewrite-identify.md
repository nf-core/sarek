# Rewrite Target Analysis

## Pipeline Summary

- **Pipeline**: `nf-core/sarek` (tag `3.8.1`)
- **Repository**: `https://github.com/nf-core/sarek`
- **Data source**: repo inspection only (`main.nf`, `workflows/`, `subworkflows/`, `modules/`)
- **Total pipeline runtime**: unknown (repo-only analysis)
- **Sample count**: unknown (assumptions below use a representative 20-sample cohort = 10 tumor-normal pairs)

### Runtime Assumptions Used For Scoring (Repo-only)

Because no trace file or Seqera run metadata was provided, ROI scoring uses normalized estimates:

- `per-sample` task multiplier: **20**
- `per-pair` task multiplier: **10**
- `global` task multiplier: **1**
- Difficulty mapping: **low=1, medium=3, high=9**
- ROI formula: **(predicted_savings_per_task_minutes x task_count) / difficulty_score**

## Tool Inventory

| Tool | Language | Per-sample? | Mean runtime | Total runtime | Source URL | Notes |
|------|----------|-------------|-------------|---------------|-----------|-------|
| `bwa/bwamem2/dragmap/sentieon_bwamem/parabricks_fq2bam` | C/C++ or proprietary binaries | Yes | very high (120-300m) | very high | `http://bio-bwa.sourceforge.net/`, `https://github.com/bwa-mem2/bwa-mem2`, `https://github.com/Illumina/dragmap` | FASTQ -> BAM/CRAM alignment; bottleneck class; high rewrite difficulty; sentieon/parabricks closed-source binaries |
| `gatk4*` (23 module processes) | Java | mixed (sample/pair/global) | low to very high (2-180m) | very high | `https://gatk.broadinstitute.org/` | Calling + recalibration + filtering; open-source but algorithmically complex; many scatter/gather tasks |
| `samtools*` (9 module processes + local reindex) | C + shell wrappers | mostly per-sample | low to medium (1-25m) | high | `http://www.htslib.org/` | BAM/CRAM conversion/index/stats/mpileup; strong bundling potential with mosdepth/indexcov |
| `mosdepth` | Rust | per-sample | medium (10-25m) | high | `https://github.com/brentp/mosdepth` | Coverage/QC; MultiQC-compatible outputs (`*.summary.txt`, dist files) |
| `goleft indexcov` | Go | per-sample | medium (8-20m) | medium-high | `https://github.com/brentp/goleft` | Coverage QC from BAM index; can be fused with QC pass |
| `fastqc` | Java | per-sample (mapping step) | medium (5-15m) | medium | `https://www.bioinformatics.babraham.ac.uk/projects/fastqc/` | MultiQC-compatible; rewriting likely low ROI vs established tool |
| `fastp` | C++ | per-sample (mapping step) | medium (5-20m) | medium-high | `https://github.com/OpenGene/fastp` | Trimming/QC, emits HTML/JSON plots; rewrite risk from parity requirements |
| `bbmap/bbsplit` | Java | per-sample if contamination removal enabled | medium-high (10-45m) | medium-high | `https://jgi.doe.gov/data-and-tools/software-tools/bbtools/` | Produces QC text used by MultiQC bbmap module |
| `cnvkit*` (6 module processes) | Python | per-sample + global reference prep | medium-high (20-90m) | high when enabled | `https://cnvkit.readthedocs.io/` | Batch/call/export chain; Python stack = rewrite-friendly compared to C/Java callers |
| `controlfreec*` (5 module processes) | C++ + Perl scripts | per-sample/per-pair | medium-high (20-90m) | medium-high | `http://boevalab.inf.ethz.ch/FREEC` | Includes plot generation (`makegraph2` PNG outputs); mixed language; medium complexity |
| `ascat` | R | per-pair | high (40-120m) | medium-high | `https://github.com/VanLoo-lab/ascat` | Statistical model + plots; difficult output parity |
| `deepvariant` | C++/TensorFlow | per-sample | very high (60-240m) | high | `https://github.com/google/deepvariant` | ML-heavy caller; high rewrite difficulty |
| `freebayes` | C++ | sample/pair | medium-high (20-90m) | medium-high | `https://github.com/freebayes/freebayes` | Often followed by filter/merge/tabix |
| `manta/strelka/tiddit/svdb` | C++ | sample/pair | medium-high to high | high | `https://github.com/Illumina/manta`, `https://github.com/Illumina/strelka`, `https://github.com/SciLifeLab/TIDDIT` | SV calling chain, difficult to faithfully reimplement |
| `mutect2 auxiliary` (`getpileupsummaries`, `gatherpileupsummaries`, `calculatecontamination`, `learnreadorientationmodel`, `filtermutectcalls`) | Java | pair/tumor-only | medium (5-45m) | high combined | `https://gatk.broadinstitute.org/` | Strong chain/bundling candidate but high algorithmic risk |
| `bcftools*` (annotate/concat/isec/merge/mpileup/norm/sort/stats/view) | C | sample/pair/global | low-medium (1-30m) | high in aggregate | `http://samtools.github.io/bcftools/` | Heavy VCF post-processing fan-out/fan-in; excellent workflow-level bundling target |
| `vcftools` | C++ | sample/pair | low-medium (2-15m) | medium | `http://vcftools.sourceforge.net/` | VCF QC stats for MultiQC |
| `snpeff/snpsift/vep` | Java + Java + Perl | sample/pair/global | medium-high (10-90m) | high if all enabled | `https://pcingola.github.io/SnpEff/`, `https://www.ensembl.org/info/docs/tools/vep/` | Annotation stack; formatting parity and plugin behavior are risks |
| `varlociraptor*` + `yte` + `rbt_vcfsplit` | Rust + template engine | sample/pair | medium-high (15-90m) | medium-high | `https://varlociraptor.github.io/`, `https://github.com/rust-bio/rust-bio-tools` | Core caller already Rust (rewrite generally not needed) |
| `msisensor2/msisensorpro` | C/C++ | sample/pair | medium (10-45m) | medium | `https://github.com/niu-lab/msisensor2`, `https://github.com/xjtu-omics/msisensor-pro` | Specialized MSI callers; limited benefit from rewrite |
| `muse` | C++ | per-pair | medium-high (20-70m) | medium | `https://github.com/wwylab/MuSE` | Two-step call+sump chain |
| `ngscheckmate` | Perl/R + binaries | per-pair/group | low-medium (5-20m) | medium | `https://github.com/parklab/NGSCheckMate` | Identity check from pileups/VCFs |
| `utility processes` (`untar`, `unzip`, `gunzip`, `cat`, `tabix`, `gawk`, local awk helpers) | shell/C utilities | mostly global or per-file | low (seconds to minutes) | medium aggregated | mixed | Individually small, but many containerized steps add scheduling/staging overhead |

### Complete Process Inventory (all process definitions)

Total process definitions found: **123** (`modules/**/main.nf`).
No `process` blocks are defined in `workflows/**` or `subworkflows/**` (these compose modules).

- **GATK family (23)**: `GATK4_APPLYBQSR`, `GATK4_APPLYVQSR`, `GATK4_BASERECALIBRATOR`, `GATK4_CALCULATECONTAMINATION`, `GATK4_CNNSCOREVARIANTS`, `GATK4_CREATESEQUENCEDICTIONARY`, `GATK4_ESTIMATELIBRARYCOMPLEXITY`, `GATK4_FILTERMUTECTCALLS`, `GATK4_FILTERVARIANTTRANCHES`, `GATK4_GATHERBQSRREPORTS`, `GATK4_GATHERPILEUPSUMMARIES`, `GATK4_GENOMICSDBIMPORT`, `GATK4_GENOTYPEGVCFS`, `GATK4_GETPILEUPSUMMARIES`, `GATK4_HAPLOTYPECALLER`, `GATK4_INTERVALLISTTOBED`, `GATK4_LEARNREADORIENTATIONMODEL`, `GATK4_MARKDUPLICATES`, `GATK4_MERGEMUTECTSTATS`, `GATK4_MERGEVCFS`, `GATK4_MUTECT2`, `GATK4_VARIANTRECALIBRATOR`, `GATK4SPARK_*` trio.
- **SAMtools family (10 including local wrapper)**: `SAMTOOLS_BAM2FQ`, `SAMTOOLS_COLLATEFASTQ`, `SAMTOOLS_CONVERT`, `SAMTOOLS_FAIDX`, `SAMTOOLS_INDEX`, `SAMTOOLS_MERGE`, `SAMTOOLS_MPILEUP`, `SAMTOOLS_STATS`, `SAMTOOLS_VIEW`, `SAMTOOLS_REINDEX_BAM`.
- **Sentieon family (9)**: `SENTIEON_APPLYVARCAL`, `SENTIEON_BWAMEM`, `SENTIEON_DEDUP`, `SENTIEON_DNAMODELAPPLY`, `SENTIEON_DNASCOPE`, `SENTIEON_GVCFTYPER`, `SENTIEON_HAPLOTYPER`, `SENTIEON_TNSCOPE`, `SENTIEON_VARCAL`.
- **BCFtools family (9)**: `BCFTOOLS_ANNOTATE`, `BCFTOOLS_CONCAT`, `BCFTOOLS_ISEC`, `BCFTOOLS_MERGE`, `BCFTOOLS_MPILEUP`, `BCFTOOLS_NORM`, `BCFTOOLS_SORT`, `BCFTOOLS_STATS`, `BCFTOOLS_VIEW`.
- **CNVkit family (6)**: `CNVKIT_ANTITARGET`, `CNVKIT_BATCH`, `CNVKIT_CALL`, `CNVKIT_EXPORT`, `CNVKIT_GENEMETRICS`, `CNVKIT_REFERENCE`.
- **Control-FREEC family (5)**: `CONTROLFREEC_FREEC`, `CONTROLFREEC_ASSESSSIGNIFICANCE`, `CONTROLFREEC_FREEC2BED`, `CONTROLFREEC_FREEC2CIRCOS`, `CONTROLFREEC_MAKEGRAPH2`.
- **FGBIO family (4 including copy UMI)**: `FGBIO_FASTQTOBAM`, `FGBIO_GROUPREADSBYUMI`, `FGBIO_CALLMOLECULARCONSENSUSREADS`, `FGBIO_COPYUMIFROMREADNAME`.
- **Alignment family**: `BWA_INDEX`, `BWA_MEM`, `BWAMEM2_INDEX`, `BWAMEM2_MEM`, `DRAGMAP_HASHTABLE`, `DRAGMAP_ALIGN`, `PARABRICKS_FQ2BAM`.
- **SV/variant callers**: `FREEBAYES`, `DEEPVARIANT_RUNDEEPVARIANT`, `LOFREQ_CALLPARALLEL`, `MANTA_GERMLINE`, `MANTA_SOMATIC`, `MANTA_TUMORONLY`, `STRELKA_GERMLINE`, `STRELKA_SOMATIC`, `TIDDIT_SV`, `MUSE_CALL`, `MUSE_SUMP`, `ASCAT`, `MSISENSOR2_MSI`, `MSISENSORPRO_SCAN`, `MSISENSORPRO_MSISOMATIC`, `GOLEFT_INDEXCOV`, `SVDB_MERGE`, `NGSCHECKMATE_NCM`.
- **Annotation and post-processing**: `SNPEFF_DOWNLOAD`, `SNPEFF_SNPEFF`, `SNPSIFT_ANNMEMCREATE`, `SNPSIFT_ANNMEM`, `ENSEMBLVEP_DOWNLOAD`, `ENSEMBLVEP_VEP`, `VCFTOOLS`, `VCFLIB_VCFFILTER`, `VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES`, `VARLOCIRAPTOR_PREPROCESS`, `VARLOCIRAPTOR_CALLVARIANTS`, `RBT_VCFSPLIT`, `YTE`.
- **Utility and local glue**: `FASTQC`, `FASTP`, `MOSDEPTH`, `BBMAP_BBSPLIT`, `CAT_CAT`, `CAT_FASTQ`, `TABIX_TABIX`, `TABIX_BGZIPTABIX`, `GUNZIP`, `UNTAR`, `UNZIP`, `SPRING_DECOMPRESS`, `GAWK`, `MULTIQC`, `ADD_INFO_TO_VCF`, `CONSENSUS_FROM_SITES`, `CREATE_INTERVALS_BED`.

### Source availability and hosting risk

- **Open-source, straightforward access (GitHub or project docs)**: bwa/bwamem2, samtools/bcftools/tabix, mosdepth, goleft, fastp, fastqc, cnvkit, freebayes, manta, strelka, tiddit, deepvariant, snpeff/snpsift, VEP, vcflib, vcftools, varlociraptor, rust-bio-tools.
- **Closed/proprietary binaries (high rewrite risk due unavailable source or licensing constraints)**: Sentieon (`sentieon/*`), NVIDIA Parabricks (`parabricks/fq2bam`).
- **Mixed script ecosystems (higher parity burden due plotting/format details)**: ASCAT (R), Control-FREEC ancillary scripts, MuSE, NGSCheckMate.

### MultiQC compatibility (important for rewrite effort)

From `assets/multiqc_config.yml`, primary modules consumed are: `fastqc`, `fastp`, `bbmap`, `picard`, `samtools`, `mosdepth`, `gatk`, `bcftools`, `vcftools`, `snpeff`, `vep`, plus custom `dedup_metrics` table.

Implication: rewrites must preserve report file names/headers for these producers:

- `FASTQC`, `FASTP`, `BBMAP_BBSPLIT`
- `GATK4_MARKDUPLICATES` or `GATK4_ESTIMATELIBRARYCOMPLEXITY`
- `SAMTOOLS_STATS`, `MOSDEPTH`
- `BCFTOOLS_STATS`, `VCFTOOLS`
- `SNPEFF_SNPEFF`, `ENSEMBLVEP_VEP`
- `SENTIEON_DEDUP` custom metrics TSV

## Ranked Rewrite Targets

### Target 1: BAM/CRAM Single-pass QC Bundle (Priority: HIGH)

- **Tools to bundle**: `SAMTOOLS_STATS`, `MOSDEPTH`, `SAMTOOLS_REINDEX_BAM`, `GOLEFT_INDEXCOV` (and optionally a light `flagstat`-style summary emitter)
- **Current total runtime**: ~35 min per sample (repo-only estimate)
- **Predicted savings**: ~22 min per sample (~63%); ~440 min across 20-sample run
- **Difficulty**: `medium` (difficulty score 3)
- **Estimated token cost**: `$250-600`
- **Shared inputs**: BAM/CRAM + index (+ optional intervals)
- **Source code**: `http://www.htslib.org/`, `https://github.com/brentp/mosdepth`, local wrappers in `modules/local/samtools/reindex_bam/main.nf`
- **Risk factors**: matching mosdepth distributions exactly; preserving MultiQC-compatible output names/columns; CRAM reference edge cases
- **Bundling rationale**: all stages reread the same alignment files; single-pass accumulator architecture removes repeated BAM I/O and several task/container boundaries
- **ROI score**: `(22 x 20) / 3 = 146.7`

### Target 2: VCF Post-processing Fusion (Priority: HIGH)

- **Tools to bundle**: `BCFTOOLS_VIEW`, `BCFTOOLS_NORM`, `BCFTOOLS_SORT`, `BCFTOOLS_CONCAT`, `BCFTOOLS_ISEC`, `TABIX_BGZIPTABIX`, local `ADD_INFO_TO_VCF`, local `CONSENSUS_FROM_SITES`
- **Current total runtime**: ~12 min per VCF stream
- **Predicted savings**: ~6 min per VCF stream (~50%); ~720 min across ~120 VCF streams in 20-sample cohort (multi-caller)
- **Difficulty**: `medium` (difficulty score 3)
- **Estimated token cost**: `$300-700`
- **Shared inputs**: VCF/BCF (+ reference for normalization)
- **Source code**: `http://samtools.github.io/bcftools/`, local scripts under `modules/local/`
- **Risk factors**: exact VCF normalization semantics and header ordering; tabix index parity; consensus-calling corner cases
- **Bundling rationale**: mostly sequential transforms with repeated compress/decompress/index cycles; merge into one streaming VCF graph
- **ROI score**: `(6 x 120) / 3 = 240.0`

### Target 3: CNVkit Python Stack Rewrite (Priority: MEDIUM)

- **Tools to bundle**: `CNVKIT_BATCH`, `CNVKIT_CALL`, `CNVKIT_EXPORT`, `CNVKIT_GENEMETRICS` (+ optionally `CNVKIT_ANTITARGET`/`CNVKIT_REFERENCE` prep)
- **Current total runtime**: ~45 min per sample
- **Predicted savings**: ~20 min per sample (~44%); ~400 min across 20 samples
- **Difficulty**: `medium` (difficulty score 3)
- **Estimated token cost**: `$250-650`
- **Shared inputs**: BAM + target/antitarget BED + CNV reference
- **Source code**: `https://cnvkit.readthedocs.io/`
- **Risk factors**: segmentation/calling numerical parity; export format compatibility used downstream
- **Bundling rationale**: Python-heavy multi-step chain with repeated I/O over same per-sample CNV intermediates
- **ROI score**: `(20 x 20) / 3 = 133.3`

### Target 4: Control-FREEC Postprocessing Bundle (Priority: MEDIUM)

- **Tools to bundle**: `CONTROLFREEC_ASSESSSIGNIFICANCE`, `CONTROLFREEC_FREEC2BED`, `CONTROLFREEC_FREEC2CIRCOS`, `CONTROLFREEC_MAKEGRAPH2` (keep `CONTROLFREEC_FREEC` as upstream caller initially)
- **Current total runtime**: ~30 min per tumor or pair output set
- **Predicted savings**: ~12 min per task (~40%); ~360 min across ~30 tasks
- **Difficulty**: `medium` (difficulty score 3)
- **Estimated token cost**: `$200-500` (higher if plot replication required)
- **Shared inputs**: FREEC ratio/CNV result tables
- **Source code**: `http://boevalab.inf.ethz.ch/FREEC`
- **Risk factors**: plot generation parity (PNG); undocumented script assumptions; coordinate formatting quirks
- **Bundling rationale**: pure sequential post-processing chain over same FREEC outputs
- **ROI score**: `(12 x 30) / 3 = 120.0`

### Target 5: Interval and Staging Utility Consolidation (Priority: MEDIUM)

- **Tools to bundle**: `CREATE_INTERVALS_BED`, `GAWK` interval splitter calls, `TABIX_BGZIPTABIX` for interval artifacts, selected `CAT_*` glue
- **Current total runtime**: low per task, but many small tasks and container startups
- **Predicted savings**: ~15 min per run from orchestration overhead reduction (compute + scheduling + staging)
- **Difficulty**: `low` (difficulty score 1)
- **Estimated token cost**: `$80-220`
- **Shared inputs**: interval lists/BED files
- **Source code**: local modules + GNU coreutils/gawk
- **Risk factors**: minimal algorithmic risk; mostly file naming and exact split behavior
- **Bundling rationale**: removes many tiny utility tasks and simplifies DAG fan-out
- **ROI score**: `(15 x 1) / 1 = 15.0`

## Bundling Groups

| Group | Tools | Shared input | Single-pass feasible? |
|-------|-------|-------------|----------------------|
| BAM-QC-1 | `SAMTOOLS_STATS`, `MOSDEPTH`, `SAMTOOLS_REINDEX_BAM`, `GOLEFT_INDEXCOV` | BAM/CRAM + index | Yes |
| VCF-POST-1 | `BCFTOOLS_VIEW/NORM/SORT/CONCAT/ISEC`, `TABIX_*`, `ADD_INFO_TO_VCF`, `CONSENSUS_FROM_SITES` | VCF/BCF streams | Yes (streaming transform pipeline) |
| CNV-1 | `CNVKIT_BATCH/CALL/EXPORT/GENEMETRICS` | BAM + target BED + CNV reference | Partially (shared parse/cache; not pure one-pass BAM) |
| FREEC-POST-1 | `FREEC2BED`, `FREEC2CIRCOS`, `ASSESSSIGNIFICANCE`, `MAKEGRAPH2` | FREEC ratio/CNV outputs | Yes |
| VARLOCIRAPTOR-CHAIN | `RBT_VCFSPLIT`, `VARLOCIRAPTOR_PREPROCESS`, `VARLOCIRAPTOR_CALLVARIANTS`, `BCFTOOLS_CONCAT/SORT` | VCF/BCF chunks + alignment properties | Partial (already Rust-heavy; lower rewrite need) |

## Recommended First Target

`BAM/CRAM Single-pass QC Bundle` - highest practical payoff with moderate difficulty and immediate reduction of repeated BAM I/O, task submissions, container pulls, and staging overhead while preserving existing callers.

## Methodology notes

- Runtime data was inferred from repository structure and tool classes (no trace/run metadata available).
- Process inventory was built by scanning all `process` definitions under `modules/**/main.nf`; workflow/subworkflow files contain orchestration only.
- Bottleneck designation is estimated from expected WGS/WES behavior: aligners and large-file callers dominate; BAM/VCF utility chains dominate cumulative overhead through repetition.
- Cost estimates include caveats from the skill guidance: numerical parity and plotting often take 3-5x longer than initial optimistic estimates.
- Tools already in compiled languages were still considered for **bundling architecture gains** (single-pass + fewer tasks), but pure algorithm rewrites are deprioritized unless orchestration savings are large.
- Closed-source/proprietary tools (`sentieon`, `parabricks`) are high-risk and were not selected as initial rewrite targets.

## Appendix: Exhaustive Process-to-Module Mapping

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

