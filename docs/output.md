# nf-core/sarek: Output <!-- omit in toc -->

## Introduction <!-- omit in toc -->

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Directory Structure](#directory-structure)
- [Preprocessing](#preprocessing)
  - [Preparation of input files (FastQ or (u)BAM)](#preparation-of-input-files-fastq-or-ubam)
    - [Clip and filter read length](#clip-and-filter-read-length)
    - [Trim adapters](#trim-adapters)
    - [Split FastQ files](#split-fastq-files)
    - [UMI consensus](#umi-consensus)
  - [Map to Reference](#map-to-reference)
    - [BWA](#bwa)
    - [BWA-mem2](#bwa-mem2)
    - [DragMap](#dragmap)
    - [Sentieon BWA mem](#sentieon-bwa-mem)
  - [Mark Duplicates](#mark-duplicates)
    - [GATK MarkDuplicates (Spark)](#gatk-markduplicates-spark)
  - [Sentieon LocusCollector and Dedup](#sentieon-locuscollector-and-dedup)
  - [Base Quality Score Recalibration](#base-quality-score-recalibration)
    - [GATK BaseRecalibrator (Spark)](#gatk-baserecalibrator-spark)
    - [GATK ApplyBQSR (Spark)](#gatk-applybqsr-spark)
  - [CSV files](#csv-files)
- [Variant Calling](#variant-calling)
  - [SNVs and small indels](#snvs-and-small-indels)
    - [bcftools](#bcftools)
    - [DeepVariant](#deepvariant)
    - [FreeBayes](#freebayes)
    - [GATK HaplotypeCaller](#gatk-haplotypecaller)
      - [GATK Germline Single Sample Variant Calling](#gatk-germline-single-sample-variant-calling)
      - [GATK Joint Germline Variant Calling](#gatk-joint-germline-variant-calling)
    - [GATK Mutect2](#gatk-mutect2)
    - [Sentieon DNAscope](#sentieon-dnascope)
      - [Sentieon DNAscope joint germline variant calling](#sentieon-dnascope-joint-germline-variant-calling)
    - [Sentieon Haplotyper](#sentieon-haplotyper)
      - [Sentieon Haplotyper joint germline variant calling](#sentieon-haplotyper-joint-germline-variant-calling)
    - [Strelka](#strelka)
    - [Lofreq](#lofreq)
  - [Structural Variants](#structural-variants)
    - [Indexcov](#indexcov)
    - [Manta](#manta)
    - [TIDDIT](#tiddit)
  - [Sample heterogeneity, ploidy and CNVs](#sample-heterogeneity-ploidy-and-cnvs)
    - [ASCAT](#ascat)
    - [CNVKit](#cnvkit)
    - [Control-FREEC](#control-freec)
  - [Microsatellite instability (MSI)](#microsatellite-instability-msi)
    - [MSIsensorPro](#msisensorpro)
  - [Concatenation](#concatenation)
- [Variant annotation](#variant-annotation)
  - [snpEff](#snpeff)
  - [VEP](#vep)
  - [BCFtools annotate](#bcftools-annotate)
- [Quality control and reporting](#quality-control-and-reporting)
  - [Quality control](#quality-control)
    - [FastQC](#fastqc)
    - [FastP](#fastp)
    - [Mosdepth](#mosdepth)
    - [NGSCheckMate](#ngscheckmate)
    - [GATK MarkDuplicates reports](#gatk-markduplicates-reports)
    - [Sentieon Dedup reports](#sentieon-dedup-reports)
    - [samtools stats](#samtools-stats)
    - [bcftools stats](#bcftools-stats)
    - [VCFtools](#vcftools)
    - [snpEff reports](#snpeff-reports)
    - [VEP reports](#vep-reports)
  - [Reporting](#reporting)
    - [MultiQC](#multiqc)
  - [Pipeline information](#pipeline-information)
- [Reference files](#reference-files)

## Directory Structure

The default directory structure is as follows

```
{outdir}
├── csv
├── multiqc
├── pipeline_info
├── preprocessing
│   ├── markduplicates
│       └── <sample>
│   ├── recal_table
│       └── <sample>
│   └── recalibrated
│       └── <sample>
├── reference
└── reports
    ├── <tool1>
    └── <tool2>
work/
.nextflow.log
```

## Preprocessing

Sarek pre-processes raw FastQ files or unmapped BAM files, based on [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

### Preparation of input files (FastQ or (u)BAM)

[FastP](https://github.com/OpenGene/fastp) is a tool designed to provide all-in-one preprocessing for FastQ files and as such is used for trimming and splitting. By default, these files are not published. However, if publishing is enabled, please be aware that these files are only published once, meaning if trimming and splitting is enabled, then the resulting files will be sharded FastQ files with trimmed reads. If only one of them is enabled then the files contain either trimmed or split reads, respectively.

#### Clip and filter read length

[FastP](https://github.com/OpenGene/fastp) enables efficient clipping of reads from either the 5' end (`--clip_r1`, `--clip_r2`) or the 3' end (`--three_prime_clip_r1`, `--three_prime_clip_r2`). Additionally, FastP allows the filtering of reads based on insert size by specifying a minimum required length with the `--length_required` parameter (default: 15bp). It is recommended to optimize these parameters according to the specific characteristics of your data.

#### Trim adapters

[FastP](https://github.com/OpenGene/fastp) supports global trimming, which means it trims all reads in the front or the tail. This function is useful since sometimes you want to drop some cycles of a sequencing run. In the current implementation in Sarek
`--detect_adapter_for_pe` is set by default which enables auto-detection of adapter sequences. For more information on how to fine-tune adapter trimming, take a look into the parameter docs.

The resulting files are intermediate and by default not kept in the final files delivered to users. Set `--save_trimmed` to enable publishing of the files in:

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/fastp/<sample>`**

- `<sample>_<lane>_{1,2}.fastp.fastq.gz>`
  - Bgzipped FastQ file

</details>

#### Split FastQ files

[FastP](https://github.com/OpenGene/fastp) supports splitting of one FastQ file into multiple files allowing parallel alignment of sharded FastQ file. To enable splitting, the number of reads per output can be specified. For more information, take a look into the parameter `--split_fastq`in the parameter docs.

These files are intermediate and by default not placed in the output-folder kept in the final files delivered to users. Set `--save_split` to enable publishing of these files to:

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/fastp/<sample>/`**

- `<sample_lane_{1,2}.fastp.fastq.gz>`
  - Bgzipped FastQ file

</details>

#### UMI consensus

Sarek can process UMI-reads, using [fgbio](http://fulcrumgenomics.github.io/fgbio/tools/latest/) tools.

These files are intermediate and by default not placed in the output-folder kept in the final files delivered to users. Set `--save_split` to enable publishing of these files to:

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/umi/<sample>/`**

- `<sample_lane_{1,2}.umi-consensus.bam>`

**Output directory: `{outdir}/reports/umi/`**

- `<sample_lane_{1,2}_umi_histogram.txt>`

</details>

### Map to Reference

#### BWA

[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome. The aligned reads are then coordinate-sorted (or name-sorted if [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html).

#### BWA-mem2

[BWA-mem2](https://github.com/bwa-mem2/bwa-mem2) is a software package for mapping low-divergent sequences against a large reference genome.The aligned reads are then coordinate-sorted (or name-sorted if [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html).

#### DragMap

[DragMap](https://github.com/Illumina/dragmap) is an open-source software implementation of the DRAGEN mapper, which the Illumina team created so that we would have an open-source way to produce the same results as their proprietary DRAGEN hardware. The aligned reads are then coordinate-sorted (or name-sorted if [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html).

These files are intermediate and by default not placed in the output-folder kept in the final files delivered to users. Set `--save_mapped` to enable publishing, furthermore add the flag `save_output_as_bam` for publishing in BAM format.

#### Sentieon BWA mem

Sentieon [bwa mem](https://support.sentieon.com/manual/usages/general/#bwa-mem-syntax) is a subroutine for mapping low-divergent sequences against a large reference genome. It is part of the proprietary software package [DNAseq](https://www.sentieon.com/detailed-description-of-pipelines/#dnaseq) from [Sentieon](https://www.sentieon.com/).

The aligned reads are coordinate-sorted with Sentieon.

<details markdown="1">
<summary>Output files for all mappers and samples</summary>

The alignment files (BAM or CRAM) produced by the chosen aligner are not published by default. CRAM output files will not be saved in the output-folder (`outdir`), unless the flag `--save_mapped` is used. BAM output can be selected by setting the flag `--save_output_as_bam`.

**Output directory: `{outdir}/preprocessing/mapped/<sample>/`**

- if `--save_mapped`: `<sample>.sorted.cram` and `<sample>.sorted.cram.crai`

  - CRAM file and index

- if `--save_mapped --save_output_as_bam`: `<sample>.sorted.bam` and `<sample>.sorted.bam.bai`
  - BAM file and index
  </details>

### Mark Duplicates

During duplicate marking, read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artificial processes are identified. These are considered to be non-independent observations, so all but a single read pair within each set of duplicates are marked, causing the marked pairs to be ignored by default during the variant discovery process.

For further reading and documentation see the [data pre-processing for variant discovery from the GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).

#### GATK MarkDuplicates (Spark)

By default, Sarek will use [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/5358880192027-MarkDuplicates-Picard-).

To use the corresponding spark implementation [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark), please specify `--use_gatk_spark markduplicates`. The resulting files are converted to CRAM with either [samtools](https://www.htslib.org/doc/samtools.html), when GATK MarkDuplicates is used, or, implicitly, by GATK MarkDuplicatesSpark.

The resulting CRAM files are delivered to the users.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/markduplicates/<sample>/`**

- `<sample>.md.cram` and `<sample>.md.cram.crai`
  - CRAM file and index
- if `--save_output_as_bam`:
  - `<sample>.md.bam` and `<sample>.md.bam.bai`

</details>

### Sentieon LocusCollector and Dedup

The subroutines LocusCollector and Dedup are part of Sentieon DNAseq packages with speedup versions of the standard GATK tools, and together those two subroutines correspond to GATK's MarkDuplicates.

The subroutine [LocusCollector](https://support.sentieon.com/manual/usages/general/#driver-algorithm-syntax) collects read information that will be used for removing or tagging duplicate reads; its output is the score file indicating which reads are likely duplicates.

The subroutine [Dedup](https://support.sentieon.com/manual/usages/general/#dedup-algorithm) marks or removes duplicate reads based on the score file supplied by LocusCollector, and produces a BAM or CRAM file.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/sentieon_dedup/<sample>/`**

- `<sample>.dedup.cram` and `<sample>.dedup.cram.crai`
  - CRAM file and index
- if `--save_output_as_bam`:
  - `<sample>.dedup.bam` and `<sample>.dedup.bam.bai`

</details>

### Base Quality Score Recalibration

During Base Quality Score Recalibration, systematic errors in the base quality scores are corrected by applying machine learning to detect and correct for them. This is important for evaluating the correct call of a variant during the variant discovery process. However, this is not needed for all combinations of tools in Sarek. Notably, this should be turned off when having UMI tagged reads or using DragMap (see [here](https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode)) as mapper.

For further reading and documentation see the [technical documentation by GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-).

#### GATK BaseRecalibrator (Spark)

[GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360042477672-BaseRecalibrator) generates a recalibration table based on various co-variates.

To use the corresponding spark implementation [GATK BaseRecalibratorSpark](https://gatk.broadinstitute.org/hc/en-us/articles/5358896138011-BaseRecalibrator), please specify `--use_gatk_spark baserecalibrator`.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/recal_table/<sample>/`**

- `<sample>.recal.table`
  - Recalibration table associated to the duplicates-marked CRAM file.

</details>

#### GATK ApplyBQSR (Spark)

[GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/5358826654875-ApplyBQSR) recalibrates the base qualities of the input reads based on the recalibration table produced by the [GATK BaseRecalibrator](#gatk-baserecalibrator) tool.

Specify `--use_gatk_spark baserecalibrator` to use [GATK ApplyBQSRSpark](https://gatk.broadinstitute.org/hc/en-us/articles/5358898266011-ApplyBQSRSpark-BETA-) instead, the respective spark implementation.

The resulting recalibrated CRAM files are delivered to the user. Recalibrated CRAM files are usually 2-3 times larger than the duplicate-marked CRAM files.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/recalibrated/<sample>/`**

- `<sample>.recal.cram` and `<sample>.recal.cram.crai`
  - CRAM file and index
- if `--save_output_as_bam`:
  - `<sample>.recal.bam` and `<sample>.recal.bam.bai` - BAM file and index
  </details>

### CSV files

The CSV files are auto-generated and can be used by Sarek for further processing and/or variant calling.

See the [`input`](usage#input-sample-sheet-configurations) section in the usage documentation for further reading and documentation on how to make the most of them.

<details markdown="1">
<summary>Output files:</summary>

**Output directory: `{outdir}/preprocessing/csv`**

- `mapped.csv`
  - if `--save_mapped`
  - CSV containing an entry for each sample with the columns `patient,sample,sex,status,bam,bai`
- `markduplicates_no_table.csv`
  - CSV containing an entry for each sample with the columns `patient,sample,sex,status,cram,crai`
- `markduplicates.csv`
  - CSV containing an entry for each sample with the columns `patient,sample,sex,status,cram,crai,table`
- `recalibrated.csv`
  - CSV containing an entry for each sample with the columns`patient,sample,sex,status,cram,crai`
- `variantcalled.csv`
  - CSV containing an entry for each sample with the columns `patient,sample,vcf`
  </details>

## Variant Calling

The results regarding variant calling are collected in `{outdir}/variantcalling/`.
If some results from a variant caller do not appear here, please check out the `--tools` section in the parameter [documentation](https://nf-co.re/sarek/latest/parameters).

(Recalibrated) CRAM files can used as an input to start the variant calling.

### SNVs and small indels

For single nucleotide variants (SNVs) and small indels, multiple tools are available for normal (germline), tumor-only, and tumor-normal (somatic) paired data. For a list of the appropriate tool(s) for the data and sequencing type at hand, please check [here](usage#which-tool).

#### bcftools

[bcftools mpileup](https://samtools.github.io/bcftools/bcftools.html#mpileup) generates pileup of a CRAM file, followed by [bcftools call](https://samtools.github.io/bcftools/bcftools.html#call) and filtered with `-i 'count(GT==\"RR\")==0`.
For further reading and documentation see the [bcftools manual](https://samtools.github.io/bcftools/howtos/variant-calling.html).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/variantcalling/bcftools/<sample>/`**

- `<sample>.bcftools.vcf.gz` and `<sample>.bcftools.vcf.gz.tbi`
  - VCF with tabix index

</details>

#### DeepVariant

[DeepVariant](https://github.com/google/deepvariant) is a deep learning-based variant caller that takes aligned reads, produces pileup image tensors from them, classifies each tensor using a convolutional neural network and finally reports the results in a standard VCF or gVCF file. For further documentation take a look [here](https://github.com/google/deepvariant/tree/r1.4/docs).

<details markdown="1">
<summary>Output files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/deepvariant/<sample>/`**

- `<sample>.deepvariant.vcf.gz` and `<sample>.deepvariant.vcf.gz.tbi`
  - VCF with tabix index
- `<sample>.deepvariant.g.vcf.gz` and `<sample>.deepvariant.g.vcf.gz.tbi`
  - gVCF with tabix index
  </details>

#### FreeBayes

[FreeBayes](https://github.com/ekg/freebayes) is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs, indels, MNPs, and complex events smaller than the length of a short-read sequencing alignment. For further reading and documentation see the [FreeBayes manual](https://github.com/ekg/freebayes/blob/master/README.md#user-manual-and-guide).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/variantcalling/freebayes/{sample,normalsample_vs_tumorsample}/`**

- `<sample>.freebayes.vcf.gz` and `<sample>.freebayes.vcf.gz.tbi`
  - VCF with tabix index

</details>

#### GATK HaplotypeCaller

[GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller) calls germline SNPs and indels via local re-assembly of haplotypes.

<details markdown="1">
<summary>Output files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/haplotypecaller/<sample>/`**

- `<sample>.haplotypecaller.vcf.gz` and `<sample>.haplotypecaller.vcf.gz.tbi`
  - VCF with tabix index

</details>

##### GATK Germline Single Sample Variant Calling

[GATK Single Sample Variant Calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)
uses HaplotypeCaller in its default single-sample mode to call variants. The VCF that HaplotypeCaller emits errors on the side of sensitivity, therefore they are filtered by first running the [CNNScoreVariants](https://gatk.broadinstitute.org/hc/en-us/articles/5358904862107-CNNScoreVariants) tool. This tool annotates each variant with a score indicating the model's prediction of the quality of each variant. To apply filters based on those scores run the [FilterVariantTranches](https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches) tool with SNP and INDEL sensitivity tranches appropriate for your task.

If the haplotype-called VCF files are not filtered, then Sarek should be run with at least one of the options `--dbsnp` or `--known_indels`.

<details markdown="1">
<summary>Output files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/haplotypecaller/<sample>/`**

- `<sample>.haplotypecaller.filtered.vcf.gz` and `<sample>.haplotypecaller.filtered.vcf.gz.tbi`
  - VCF with tabix index

</details>

##### GATK Joint Germline Variant Calling

[GATK Joint germline Variant Calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) uses Haplotypecaller per sample in `gvcf` mode. Next, the gVCFs are consolidated from multiple samples into a [GenomicsDB](https://gatk.broadinstitute.org/hc/en-us/articles/5358869876891-GenomicsDBImport) datastore. After joint [genotyping](https://gatk.broadinstitute.org/hc/en-us/articles/5358906861083-GenotypeGVCFs), [VQSR](https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator) is applied for filtering to produce the final multisample callset with the desired balance of precision and sensitivity.

<details markdown="1">
<summary>Output files from joint germline variant calling</summary>

**Output directory: `{outdir}/variantcalling/haplotypecaller/<sample>/`**

- `<sample>.haplotypecaller.g.vcf.gz` and `<sample>.haplotypecaller.g.vcf.gz.tbi`
  - gVCF with tabix index

**Output directory: `{outdir}/variantcalling/haplotypecaller/joint_variant_calling/`**

- `joint_germline.vcf.gz` and `joint_germline.vcf.gz.tbi`
  - VCF with tabix index
- `joint_germline_recalibrated.vcf.gz` and `joint_germline_recalibrated.vcf.gz.tbi`
  - variant recalibrated VCF with tabix index (if VQSR is applied)

</details>

#### GATK Mutect2

[GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2) calls somatic SNVs and indels via local assembly of haplotypes.
When `--joint_mutect2` is used, Mutect2 subworkflow outputs will be saved in a subfolder named with the patient ID and `{patient}.mutect2.vcf.gz` file will contain variant calls from all of the normal and tumor samples of the patient.
For further reading and documentation see the [Mutect2 manual](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132).
It is not required, but recommended to have a [panel of normals (PON)](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON) using at least 40 normal samples to get filtered somatic calls. When using `--genome GATK.GRCh38`, a panel-of-normals file is available. However, it is _highly_ recommended to create one matching your tumor samples. Creating your own panel-of-normals is currently not natively supported by the pipeline. See [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132) for how to create one manually.

<details markdown="1">
<summary>Output files for tumor-only and tumor/normal paired samples</summary>

**Output directory: `{outdir}/variantcalling/mutect2/{sample,tumorsample_vs_normalsample,patient}/`**

Files created:

- `{sample,tumorsample_vs_normalsample,patient}.mutect2.vcf.gz` and `{sample,tumorsample_vs_normalsample,patient}.mutect2.vcf.gz.tbi`
  - unfiltered (raw) Mutect2 calls VCF with tabix index
- `{sample,tumorsample_vs_normalsample,patient}.mutect2.vcf.gz.stats`
  - a stats file generated during calling of raw variants (needed for filtering)
- `{sample,tumorsample_vs_normalsample}.mutect2.contamination.table`
  - table calculating the fraction of reads coming from cross-sample contamination
- `{sample,tumorsample_vs_normalsample}.mutect2.segmentation.table`
  - table containing segmentation of the tumor by minor allele fraction
- `{sample,tumorsample_vs_normalsample,patient}.mutect2.artifactprior.tar.gz`
  - prior probabilities for read orientation artifacts
- `{sample,tumorsample,normalsample}.mutect2.pileups.table`
  - tabulates pileup metrics for inferring contamination
- `{sample,tumorsample_vs_normalsample,patient}.mutect2.filtered.vcf.gz` and `{sample,tumorsample_vs_normalsample,patient}.mutect2.filtered.vcf.gz.tbi`
  - filtered Mutect2 calls VCF with tabix index based on the probability that a variant is somatic
- `{sample,tumorsample_vs_normalsample,patient}.mutect2.filtered.vcf.gz.filteringStats.tsv`
  - a stats file generated during the filtering of Mutect2 called variants

</details>

#### Sentieon DNAscope

[Sentieon DNAscope](https://support.sentieon.com/appnotes/dnascope_ml/#dnascope-germline-variant-calling-with-a-machine-learning-model) is a variant-caller which aims at outperforming GATK's Haplotypecaller in terms of both speed and accuracy. DNAscope allows you to use a machine learning model to perform variant calling with higher accuracy by improving the candidate detection and filtering.

<details markdown="1">
<summary>Unfiltered VCF-files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/sentieon_dnascope/<sample>/`**

- `<sample>.dnascope.unfiltered.vcf.gz` and `<sample>.dnascope.unfiltered.vcf.gz.tbi`
  - VCF with tabix index

</details>

The output from Sentieon's DNAscope can be controlled through the option `--sentieon_dnascope_emit_mode` for Sarek, see [Basic usage of Sentieon functions](#basic-usage-of-sentieon-functions).

Unless `dnascope_filter` is listed under `--skip_tools` in the nextflow command, Sentieon's [DNAModelApply](https://support.sentieon.com/manual/usages/general/#dnamodelapply-algorithm) is applied to the unfiltered VCF-files in order to obtain filtered VCF-files.

<details markdown="1">
<summary>Filtered VCF-files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/sentieon_dnascope/<sample>/`**

- `<sample>.dnascope.filtered.vcf.gz` and `<sample>.dnascope.filtered.vcf.gz.tbi`
  - VCF with tabix index

</details>

##### Sentieon DNAscope joint germline variant calling

In Sentieon's package DNAscope, joint germline variant calling is done by first running Sentieon's Dnacope in emit-mode `gvcf` for each sample and then running Sentieon's [GVCFtyper](https://support.sentieon.com/manual/usages/general/#gvcftyper-algorithm) on the set of gVCF-files. See [Basic usage of Sentieon functions](#basic-usage-of-sentieon-functions) for information on how joint germline variant calling can be done in Sarek using Sentieon's DNAscope.

<details markdown="1">
<summary>Output files from joint germline variant calling</summary>

**Output directory: `{outdir}/variantcalling/sentieon_dnascope/<sample>/`**

- `<sample>.dnascope.g.vcf.gz` and `<sample>.dnascope.g.vcf.gz.tbi`
  - VCF with tabix index

**Output directory: `{outdir}/variantcalling/sentieon_dnascope/joint_variant_calling/`**

- `joint_germline.vcf.gz` and `joint_germline.vcf.gz.tbi`
  - VCF with tabix index

</details>

#### Sentieon Haplotyper

[Sentieon Haplotyper](https://support.sentieon.com/manual/usages/general/#haplotyper-algorithm) is Sention's speedup version of GATK's Haplotypecaller (see above).

<details markdown="1">
<summary>Unfiltered VCF-files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/sentieon_haplotyper/<sample>/`**

- `<sample>.haplotyper.unfiltered.vcf.gz` and `<sample>.haplotyper.unfiltered.vcf.gz.tbi`
  - VCF with tabix index

</details>

The output from Sentieon's Haplotyper can be controlled through the option `--sentieon_haplotyper_emit_mode` for Sarek, see [Basic usage of Sentieon functions](#basic-usage-of-sentieon-functions).

Unless `haplotyper_filter` is listed under `--skip_tools` in the nextflow command, GATK's CNNScoreVariants and FilterVariantTranches (see above) is applied to the unfiltered VCF-files in order to obtain filtered VCF-files.

<details markdown="1">
<summary>Filtered VCF-files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/sentieon_haplotyper/<sample>/`**

- `<sample>.haplotyper.filtered.vcf.gz` and `<sample>.haplotyper.filtered.vcf.gz.tbi`
  - VCF with tabix index

</details>

##### Sentieon Haplotyper joint germline variant calling

In Sentieon's package DNAseq, joint germline variant calling is done by first running Sentieon's Haplotyper in emit-mode `gvcf` for each sample and then running Sentieon's [GVCFtyper](https://support.sentieon.com/manual/usages/general/#gvcftyper-algorithm) on the set of gVCF-files. See [Basic usage of Sentieon functions](#basic-usage-of-sentieon-functions) for information on how joint germline variant calling can be done in Sarek using Sentieon's DNAseq. After joint genotyping, Sentieon's version of VQSR ([VarCal](https://support.sentieon.com/manual/usages/general/#varcal-algorithm) and [ApplyVarCal](https://support.sentieon.com/manual/usages/general/#applyvarcal-algorithm)) is applied for filtering to produce the final multisample callset with the desired balance of precision and sensitivity.

<details markdown="1">
<summary>Output files from joint germline variant calling</summary>

**Output directory: `{outdir}/variantcalling/sentieon_haplotyper/<sample>/`**

- `<sample>.haplotyper.g.vcf.gz` and `<sample>.haplotyper.g.vcf.gz.tbi`
  - VCF with tabix index

**Output directory: `{outdir}/variantcalling/sentieon_haplotyper/joint_variant_calling/`**

- `joint_germline.vcf.gz` and `joint_germline.vcf.gz.tbi`
  - VCF with tabix index
- `joint_germline_recalibrated.vcf.gz` and `joint_germline_recalibrated.vcf.gz.tbi`
  - variant recalibrated VCF with tabix index (if VarCal is applied)

</details>

#### Strelka

[Strelka](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. For further reading and documentation see the [Strelka user guide](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md). If [Strelka](https://github.com/Illumina/strelka) is used for somatic variant calling and [Manta](https://github.com/Illumina/manta) is also specified in tools, the output candidate indels from [Manta](https://github.com/Illumina/manta) are used according to [Strelka Best Practices](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic-configuration-example).
For further downstream analysis, take a look [here](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#interpreting-the-germline-multi-sample-variants-vcf).

<details markdown="1">
<summary>Output files for single samples (normal)</summary>

**Output directory: `{outdir}/variantcalling/strelka/<sample>/`**

- `<sample>.strelka.genome.vcf.gz` and `<sample>.strelka.genome.vcf.gz.tbi`
  - genome VCF with tabix index
- `<sample>.strelka.variants.vcf.gz` and `<sample>.strelka.variants.vcf.gz.tbi`
  - VCF with tabix index with all potential variant loci across the sample. Note this file includes non-variant loci if they have a non-trivial level of variant evidence or contain one or more alleles for which genotyping has been forced.
  </details>

<details markdown="1">
<summary>Output files for tumor/normal paired samples</summary>

**Output directory: `{outdir}/variantcalling/strelka/<tumorsample_vs_normalsample>/`**

- `<tumorsample_vs_normalsample>.strelka.somatic_indels.vcf.gz` and `<tumorsample_vs_normalsample>.strelka.somatic_indels.vcf.gz.tbi`
  - VCF with tabix index with all somatic indels inferred in the tumor sample.
- `<tumorsample_vs_normalsample>.strelka.somatic_snvs.vcf.gz` and `<tumorsample_vs_normalsample>.strelka.somatic_snvs.vcf.gz.tbi`
  - VCF with tabix index with all somatic SNVs inferred in the tumor sample.

</details>

#### Lofreq

[Lofreq](https://github.com/CSB5/lofreq) is a fast and sensitive variant-caller for inferring SNVs and indels from next-generation sequencing data. It makes full use of base-call qualities and other sources of errors inherent in sequencing, which are usually ignored by other methods or only used for filtering. For further reading and documentation see the [Lofreq user guide](https://csb5.github.io/lofreq/).

<details markdown = "1">
<summary>Output files for tumor-only samples</summary>

**Output directory: `{outdir}/variant_calling/lofreq/<sample>/`**

-`<tumorsample>.vcf.gz`
-VCF which provides a detailed description of the detected genetic variants.

  </details>

### Structural Variants

#### indexcov

[indexcov](https://github.com/brentp/goleft/tree/master/indexcov) quickly estimate coverage from a whole-genome bam or cram index.
A bam index has 16KB resolution and it is used as a coverage estimate .
The output is scaled to around 1. So a long stretch with values of 1.5 would be a heterozygous duplication. This is useful as a quick QC to get coverage values across the genome.

**Output directory: `{outdir}/variantcalling/indexcov/`**

In addition to the interactive HTML files, `indexcov` outputs a number of text files:

- `<sample>-indexcov.ped`: a .ped/.fam file with the inferred sex in the appropriate column if the sex chromosomes were found.
  the CNX and CNY columns indicating the floating-point estimate of copy-number for those chromosomes.
  `bins.out`: how many bins had a coverage value outside of (0.85, 1.15). high values can indicate high-bias samples.
  `bins.lo`: number of bins with value < 0.15. high values indicate missing data.
  `bins.hi`: number of bins with value > 1.15.
  `bins.in`: number of bins with value inside of (0.85, 1.15)
  `p.out`: `bins.out/bins.in`
  `PC1...PC5`: PCA projections calculated with depth of autosomes.

- `<sample>-indexcov.roc`: tab-delimited columns of chrom, scaled coverage cutoff, and $n_samples columns where each indicates the
  proportion of 16KB blocks at or above that scaled coverage value.
- `<sample>-indexcov.bed.gz`: a bed file with columns of chrom, start, end, and a column per sample where the values indicate there
  scaled coverage for that sample in that 16KB chunk.

#### Manta

[Manta](https://github.com/Illumina/manta) calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
[Manta](https://github.com/Illumina/manta) provides a candidate list for small indels that can be fed to [Strelka](https://github.com/Illumina/strelka) following [Strelka Best Practices](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic-configuration-example). For further reading and documentation see the [Manta user guide](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md).

<details markdown="1">
<summary>Output files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/manta/<sample>/`**

- `<sample>.manta.diploid_sv.vcf.gz` and `<sample>.manta.diploid_sv.vcf.gz.tbi`
  - VCF with tabix index containing SVs and indels scored and genotyped under a diploid model for the sample.
  </details>

<details markdown="1">
<summary>Output files for tumor-only samples</summary>

**Output directory: `{outdir}/variantcalling/manta/<sample>/`**

- `<sample>.manta.tumor_sv.vcf.gz` and `<sample>.manta.tumor_sv.vcf.gz.tbi`
  - VCF with tabix index containing a subset of the candidateSV.vcf.gz file after removing redundant candidates and small indels less than the minimum scored variant size (50 by default). The SVs are not scored, but include additional details: (1) paired and split read supporting evidence counts for each allele (2) a subset of the filters from the scored tumor-normal model are applied to the single tumor case to improve precision.
  </details>

<details markdown="1">
<summary>Output files for tumor/normal paired samples</summary>

**Output directory: `{outdir}/variantcalling/manta/<tumorsample_vs_normalsample>/`**

- `<tumorsample_vs_normalsample>.manta.diploid_sv.vcf.gz` and `<tumorsample_vs_normalsample>.manta.diploid_sv.vcf.gz.tbi`
  - VCF with tabix index containing SVs and indels scored and genotyped under a diploid model for the sample. In the case of a tumor/normal subtraction, the scores in this file do not reflect any information from the tumor sample.
- `<tumorsample_vs_normalsample>.manta.somatic_sv.vcf.gz` and `<tumorsample_vs_normalsample>.manta.somatic_sv.vcf.gz.tbi`
  - VCF with tabix index containing SVs and indels scored under a somatic variant model.
  </details>

#### TIDDIT

[TIDDIT](https://github.com/SciLifeLab/TIDDIT) identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions. For further reading and documentation see the [TIDDIT manual](https://github.com/SciLifeLab/TIDDIT/blob/master/README.md).

<details markdown="1">
<summary>Output files for normal and tumor-only samples</summary>

**Output directory: `{outdir}/variantcalling/tiddit/<sample>/`**

- `<sample>.tiddit.vcf.gz` and `<sample>.tiddit.vcf.gz.tbi`
  - VCF with tabix index containing SV calls
- `<sample>.tiddit.ploidies.tab`
  - tab file describing the estimated ploidy and coverage across each contig

</details>

<details markdown="1">
<summary>Output files for tumor/normal paired samples</summary>

**Output directory: `{outdir}/variantcalling/tiddit/<tumorsample_vs_normalsample>/`**

- `<tumorsample_vs_normalsample>.tiddit.normal.vcf.gz` and `<tumorsample_vs_normalsample>.tiddit.normal.vcf.gz.tbi`
  - VCF with tabix index containing SV calls
- `<tumorsample_vs_normalsample>.tiddit.tumor.vcf.gz` and `<tumorsample_vs_normalsample>.tiddit.tumor.vcf.gz.tbi`
  - VCF with tabix index containing SV calls
- `<tumorsample_vs_normalsample>_sv_merge.tiddit.vcf.gz` and `<tumorsample_vs_normalsample>_sv_merge.tiddit.vcf.gz.tbi`
  - merged tumor/normal VCF with tabix index
- `<tumorsample_vs_normalsample>.tiddit.ploidies.tab`
  - tab file describing the estimated ploidy and coverage across each contig

</details>

### Sample heterogeneity, ploidy and CNVs

#### ASCAT

[ASCAT](https://github.com/VanLoo-lab/ascat) is a software for performing allele-specific copy number analysis of tumor samples and for estimating tumor ploidy and purity (normal contamination).
It infers tumor purity and ploidy and calculates whole-genome allele-specific copy number profiles.
The [ASCAT](https://github.com/VanLoo-lab/ascat) process gives several images as output, described in detail in this [book chapter](http://www.ncbi.nlm.nih.gov/pubmed/22130873).
Running ASCAT on NGS data requires that the BAM files are converted into BAF and LogR values.
This is done internally using the software [AlleleCount](https://github.com/cancerit/alleleCount). For further reading and documentation see the [ASCAT manual](https://www.crick.ac.uk/research/labs/peter-van-loo/software).

<details markdown="1">
<summary>Output files for tumor/normal paired samples</summary>

**Output directory: `{outdir}/variantcalling/ascat/<tumorsample_vs_normalsample>/`**

- `<tumorsample_vs_normalsample>.tumour.ASCATprofile.png`
  - image with information about allele-specific copy number profile
- `<tumorsample_vs_normalsample>.tumour.ASPCF.png`
  - image with information about allele-specific copy number segmentation
- `<tumorsample_vs_normalsample>.before_correction_Tumour.<tumorsample_vs_normalsample>.tumour.png`
  - image with information about raw profile of tumor sample of logR and BAF values before GC correction
- `<tumorsample_vs_normalsample>.before_correction_Tumour.<tumorsample_vs_normalsample>.germline.png`
  - image with information about raw profile of normal sample of logR and BAF values before GC correction
- `<tumorsample_vs_normalsample>.after_correction_GC_Tumour.<tumorsample_vs_normalsample>.tumour.png`
  - image with information about GC and RT corrected logR and BAF values of tumor sample after GC correction
- `<tumorsample_vs_normalsample>.after_correction_GC_Tumour.<tumorsample_vs_normalsample>.germline.png`
  - image with information about GC and RT corrected logR and BAF values of normal sample after GC correction
- `<tumorsample_vs_normalsample>.tumour.sunrise.png`
  - image visualising the range of ploidy and tumor percentage values
- `<tumorsample_vs_normalsample>.metrics.txt`
  - file with information about different metrics from ASCAT profiles
- `<tumorsample_vs_normalsample>.cnvs.txt`
  - file with information about CNVS
- `<tumorsample_vs_normalsample>.purityploidy.txt`
  - file with information about purity and ploidy
- `<tumorsample_vs_normalsample>.segments.txt`
  - file with information about copy number segments
- `<tumorsample_vs_normalsample>.tumour_tumourBAF.txt` and `<tumorsample_vs_normalsample>.tumour_normalBAF.txt`
  - file with beta allele frequencies
- `<tumorsample_vs_normalsample>.tumour_tumourLogR.txt` and `<tumorsample_vs_normalsample>.tumour_normalLogR.txt`
  - file with total copy number on a logarithmic scale

The text file `<tumorsample_vs_normalsample>.cnvs.txt` contains predictions about copy number state for all the segments.
The output is a tab delimited text file with the following columns:

- _chr_: chromosome number
- _startpos_: start position of the segment
- _endpos_: end position of the segment
- _nMajor_: number of copies of one of the allels (for example the chromosome inherited of one parent)
- _nMinor_: number of copies of the other allele (for example the chromosome inherited of the other parent)

The file `<tumorsample_vs_normalsample>.cnvs.txt` contains all segments predicted by ASCAT, both those with normal copy number (nMinor = 1 and nMajor =1) and those corresponding to copy number aberrations.

</details>

#### CNVKit

[CNVKit](https://cnvkit.readthedocs.io/en/stable/) is a toolkit to infer and visualize copy number from high-throughput DNA sequencing data. It is designed for use with hybrid capture, including both whole-exome and custom target panels, and short-read sequencing platforms such as Illumina. For further reading and documentation, see the [CNVKit Documentation](https://cnvkit.readthedocs.io/en/stable/plots.html)

<details markdown="1">
<summary>Output files for normal and tumor-only samples</summary>

**Output directory: `{outdir}/variantcalling/cnvkit/<sample>/`**

- `<sample>.antitargetcoverage.cnn`
  - file containing coverage information
- `<sample>.targetcoverage.cnn`
  - file containing coverage information
- `<sample>-diagram.pdf`
  - file with plot of copy numbers or segments on chromosomes
- `<sample>-scatter.png`
  - file with plot of bin-level log2 coverages and segmentation calls
- `<sample>.bintest.cns`
  - file containing copy number segment information
- `<sample>.cnr`
  - file containing copy number ratio information
- `<sample>.cns`
  - file containing copy number segment information
- `<sample>.call.cns`
  - file containing copy number segment information
- `<sample>.genemetrics.tsv`
  - file containing per gene copy number information (if input files are annotated)
  </details>

<details markdown="1">
<summary>Output files for tumor/normal samples</summary>

**Output directory: `{outdir}/variantcalling/cnvkit/<tumorsample_vs_normalsample>/`**

- `<normalsample>.antitargetcoverage.cnn`
  - file containing coverage information
- `<normalsample>.targetcoverage.cnn`
  - file containing coverage information
- `<tumorsample>.antitargetcoverage.cnn`
  - file containing coverage information
- `<tumorsample>.targetcoverage.cnn`
  - file containing coverage information
- `<tumorsample>.bintest.cns`
  - file containing copy number segment information
- `<tumorsample>-scatter.png`
  - file with plot of bin-level log2 coverages and segmentation calls
- `<tumorsample>-diagram.pdf`
  - file with plot of copy numbers or segments on chromosomes
- `<tumorsample>.cnr`
  - file containing copy number ratio information
- `<tumorsample>.cns`
  - file containing copy number segment information
- `<tumorsample>.call.cns`
  - file containing copy number segment information
- `<tumorsample>.genemetrics.tsv`
  - file containing per gene copy number information (if input files are annotated)
  </details>

#### Control-FREEC

[Control-FREEC](https://github.com/BoevaLab/FREEC) is a tool for detection of copy-number changes and allelic imbalances (including loss of heterozygoity (LOH)) using deep-sequencing data.
[Control-FREEC](https://github.com/BoevaLab/FREEC) automatically computes, normalizes, segments copy number and beta allele frequency profiles, then calls copy number alterations and LOH.
It also detects subclonal gains and losses and evaluates the most likely average ploidy of the sample. For further reading and documentation see the [Control-FREEC Documentation](http://boevalab.inf.ethz.ch/FREEC/tutorial.html).

<details markdown="1">
<summary>Output files for tumor-only and tumor/normal paired samples</summary>

**Output directory: `{outdir}/variantcalling/controlfreec/{tumorsample,tumorsample_vs_normalsample}/`**

- `config.txt`
  - Configuration file used to run Control-FREEC
- `<tumorsample>_BAF.png` and `<tumorsample_vs_normalsample>_BAF.png`
  - image of BAF plot
- `<tumorsample>_ratio.log2.png` and `<tumorsample_vs_normalsample>_ratio.log2.png`
  - image of ratio log2 plot
- `<tumorsample>_ratio.png` and `<tumorsample_vs_normalsample>_ratio.png`
  - image of ratio plot
- `<tumorsample>.bed` and `<tumorsample_vs_normalsample>.bed`
  - translated output to a .BED file (so to view it in the UCSC Genome Browser)
- `<tumorsample>.circos.txt` and `<tumorsample_vs_normalsample>.circos.txt`
  - translated output to the Circos format
- `<tumorsample>.p.value.txt` and `<tumorsample_vs_normalsample>.p.value.txt`
  - CNV file containing p_values for each call
- `<tumorsample>_BAF.txt` and `<tumorsample_vs_normalsample>.mpileup.gz_BAF.txt`
  - file with beta allele frequencies for each possibly heterozygous SNP position
- `<tumorsample_vs_normalsample>.tumor.mpileup.gz_CNVs`
  - file with coordinates of predicted copy number alterations
- `<tumorsample>_info.txt` and `<tumorsample_vs_normalsample>.tumor.mpileup.gz_info.txt`
  - parsable file with information about FREEC run
- ` <tumorsample>_ratio.BedGraph` and `<tumorsample_vs_normalsample>.tumor.mpileup.gz_ratio.BedGraph `
  - file with ratios in BedGraph format for visualization in the UCSC genome browser. The file contains tracks for normal copy number, gains and losses, and copy neutral LOH (\*).
- `<tumorsample>_ratio.txt` and `<tumorsample_vs_normalsample>.tumor.mpileup.gz_ratio.txt`
  - file with ratios and predicted copy number alterations for each window
- `<tumorsample>_sample.cpn` and `<tumorsample_vs_normalsample>.tumor.mpileup.gz_sample.cpn`
  - files with raw copy number profiles for the tumor sample
- `<tumorsample_vs_normalsample>.normal.mpileup.gz_control.cpn`
  - files with raw copy number profiles for the control sample
- `<GC_profile.<tumorsample>.cpn>`
  - file with GC-content profile

</details>

### Microsatellite instability (MSI)

[Microsatellite instability](https://en.wikipedia.org/wiki/Microsatellite_instability) is a genetic condition associated with deficiencies in the mismatch repair (MMR) system which causes a tendency to accumulate a high number of mutations (SNVs and indels).
An altered distribution of microsatellite length is associated with a missed replication slippage which would be corrected under normal MMR conditions.

#### MSIsensorPro

[MSIsensorPro](https://github.com/xjtu-omics/msisensor-pro) is a tool to detect the MSI status of a tumor scanning the length of the microsatellite regions.
It requires a normal sample for each tumour to differentiate the somatic and germline cases. For further reading see the [MSIsensor paper](https://www.ncbi.nlm.nih.gov/pubmed/24371154).

<details markdown="1">
<summary>Output files for tumor/normal paired samples</summary>

**Output directory: `{outdir}/variantcalling/msisensor/<tumorsample_vs_normalsample>/`**

- `<tumorsample_vs_normalsample>`
  - MSI score output, contains information about the number of somatic sites.
- `<tumorsample_vs_normalsample>_dis`
  - The normal and tumor length distribution for each microsatellite position.
- `<tumorsample_vs_normalsample>_germline`
  - Somatic sites detected.
- `<tumorsample_vs_normalsample>_somatic`
  - Germline sites detected.
  </details>

### Concatenation

Germline VCFs from `DeepVariant`, `FreeBayes`, `HaplotypeCaller`, `Haplotyper`, `Manta`, `bcftools mpileup`, `Strelka`, or `Tiddit` are concatenated with `bcftools concat`. The field `SOURCE` is added to the VCF header to report the variant caller.

<details markdown="1">
<summary>Concatenated VCF-files for normal samples</summary>

**Output directory: `{outdir}/variantcalling/concat/<sample>/`**

- `<sample>.germline.vcf.gz` and `<sample>.germline.vcf.gz.tbi`
  - VCF with tabix index

</details>

## Variant annotation

This directory contains results from the final annotation steps: two tools are used for annotation, [snpEff](http://snpeff.sourceforge.net/) and [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html). Both results can also be combined by setting `--tools merge`.
All variants present in the called VCF files are annotated. For some variant callers this can mean that the variants are already filtered by `PASS`, for some this needs to be done during post-processing.

### snpEff

[snpeff](http://snpeff.sourceforge.net/) is a genetic variant annotation and effect prediction toolbox.
It annotates and predicts the effects of variants on genes (such as amino acid changes) using multiple databases for annotations.
The generated VCF header contains the software version and the used command line. For further reading and documentation see the [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#outputSummary).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/annotation/{sample,tumorsample_vs_normalsample}`**

- `{sample,tumorsample_vs_normalsample}.<variantcaller>_snpEff.ann.vcf.gz` and `{sample,tumorsample_vs_normalsample}.<variantcaller>_snpEff.ann.vcf.gz.tbi`
  - VCF with tabix index
  </details>

### VEP

[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html), based on `Ensembl`, is a tool to determine the effects of all sorts of variants, including SNPs, indels, structural variants, CNVs.
The generated VCF header contains the software version, also the version numbers for additional databases like [Clinvar](https://www.ncbi.nlm.nih.gov/clinvar/) or [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) used in the [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) line.
The format of the [consequence annotations](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html) is also in the VCF header describing the `INFO` field.
For further reading and documentation see the [VEP manual](https://www.ensembl.org/info/docs/tools/vep/index.html).

Currently, it contains:

- _Consequence_: impact of the variation, if there is any
- _Codons_: the codon change, i.e. cGt/cAt
- _Amino_acids_: change in amino acids, i.e. R/H if there is any
- _Gene_: ENSEMBL gene name
- _SYMBOL_: gene symbol
- _Feature_: actual transcript name
- _EXON_: affected exon
- _PolyPhen_: prediction based on [PolyPhen](http://genetics.bwh.harvard.edu/pph2/)
- _SIFT_: prediction by [SIFT](http://sift.bii.a-star.edu.sg/)
- _Protein_position_: Relative position of amino acid in protein
- _BIOTYPE_: Biotype of transcript or regulatory feature

plus any additional filed selected via the plugins: [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP), [LOFTEE](https://github.com/konradjk/loftee), [SpliceAI](https://spliceailookup.broadinstitute.org/), [SpliceRegion](https://www.ensembl.info/2018/10/26/cool-stuff-the-vep-can-do-splice-site-variant-annotation/).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/annotation/{sample,tumorsample_vs_normalsample}`**

- `{sample,tumorsample_vs_normalsample}.<variantcaller>_VEP.ann.vcf.gz` and `{sample,tumorsample_vs_normalsample}.<variantcaller>_VEP.ann.vcf.gz.tbi`
  - VCF with tabix index

</details>

### BCFtools annotate

[BCFtools annotate](https://samtools.github.io/bcftools/bcftools.html#annotate) is used to add annotations to VCF files. The annotations are added to the INFO column of the VCF file. The annotations are added to the VCF header and the VCF header is updated with the new annotations. For further reading and documentation see the [BCFtools annotate manual](https://samtools.github.io/bcftools/bcftools.html#annotate).

<details markdown="1">
<summary>Output files for all samples</summary>

- `{sample,tumorsample_vs_normalsample}.<variantcaller>_bcf.ann.vcf.gz` and `{sample,tumorsample_vs_normalsample}.<variantcaller>_bcf.ann.vcf.gz.tbi`
  - VCF with tabix index

</details>

## Quality control and reporting

### Quality control

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

The plots display:

- Sequence counts for each sample.
- Sequence Quality Histograms: The mean quality value across each base position in the read.
- Per Sequence Quality Scores: The number of reads with average quality scores. Shows if a subset of reads has poor quality.
- Per Base Sequence Content: The proportion of each base position for which each of the four normal DNA bases has been called.
- Per Sequence GC Content: The average GC content of reads. Normal random library typically have a roughly normal distribution of GC content.
- Per Base N Content: The percentage of base calls at each position for which an N was called.
- Sequence Length Distribution.
- Sequence Duplication Levels: The relative level of duplication found for each sequence.
- Overrepresented sequences: The total amount of overrepresented sequences found in each library.
- Adapter Content: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/fastqc/<sample-lane>`**

- `<sample-lane_1>_fastqc.html` and `<sample-lane_2>_fastqc.html`
  - [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) report containing quality metrics for your untrimmed raw FastQ files
- `<sample-lane_1>_fastqc.zip` and `<sample-lane_2>_fastqc.zip`
  - Zip archive containing the [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) report, tab-delimited data file and plot images

:::note
The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads.
They may contain adapter sequence and potentially regions with low quality.
:::

</details>

#### FastP

[FastP](https://github.com/OpenGene/fastp) is a tool designed to provide all-in-one preprocessing for FastQ files and is used for trimming and splitting. The tool then determines QC metrics for the processed reads.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/fastp/<sample>`**

- `<sample-lane>_fastp.html`
  - report in HTML format
- `<sample-lane>_fastp.json`
  - report in JSON format
- `<sample-lane>_fastp.log`
  - FastQ log file

</details>

#### Mosdepth

[Mosdepth](https://github.com/brentp/mosdepth) reports information for the evaluation of the quality of the provided alignment data.
In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced.
For further reading and documentation see the [Mosdepth documentation](https://github.com/brentp/mosdepth).

Plots will show:

- cumulative coverage distribution
- absolute coverage distribution
- average coverage per contig/chromosome

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/mosdepth/<sample>`**

- `<sample>.{sorted,md,recal}.mosdepth.global.dist.txt`
  - file used by [MultiQC](https://multiqc.info/), if `.region` file does not exist
- `<sample>.{sorted,md,recal}.mosdepth.region.dist.txt`
  - file used by [MultiQC](https://multiqc.info/)
- `<sample>.{sorted,md,recal}.mosdepth.summary.txt`
  -A summary of mean depths per chromosome and within specified regions per chromosome.
- `<sample>.{sorted,md,recal}.{per-base,regions}.bed.gz`
  - per-base depth for targeted data, per-window (500bp) depth of WGS
- `<sample>.{sorted,md,recal}.regions.bed.gz.csi`
  - CSI index for per-base depth for targeted data, per-window (500bp) depth of WGS
  </details>

#### NGSCheckMate

[NGSCheckMate](https://github.com/parklab/NGSCheckMate) is a tool for determining whether samples come from the same genetic individual, using a set of commonly heterozygous SNPs. This enables for the detecting of sample mislabelling events. The output includes a text file indicating whether samples have matched or not according to the algorithm, as well as a dendrogram visualising these results.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/ngscheckmate/`**

- `ngscheckmate_all.txt`
  - Tab delimited text file listing all the comparisons made, whether they were considered as a match, with the correlation and a normalised depth.
- `ngscheckmate_matched.txt`
  - Tab delimited text file listing only the comparison that were considered to match, with the correlation and a normalised depth.
- `ngscheckmate_output_corr_matrix.txt`
  - Tab delimited text file containing a matrix of all correlations for all comparisons made.
- `vcfs/<sample>.vcf.gz`
  - Set of vcf files for each sample. Contains calls for the set of SNP positions used to calculate sample relatedness.
  </details>

#### GATK MarkDuplicates reports

More information in the [GATK MarkDuplicates section](#gatk-markduplicates)

Duplicates can arise during sample preparation _e.g._ library construction using PCR.
Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.
These duplication artifacts are referred to as optical duplicates. If [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/5358880192027-MarkDuplicates-Picard-) is used, the metrics file generated by the tool is used, if [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) is used the report is generated by [GATK4 EstimateLibraryComplexity](https://gatk.broadinstitute.org/hc/en-us/articles/5358838684187-EstimateLibraryComplexity-Picard-) on the mapped BAM files.
For further reading and documentation see the [MarkDuplicates manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_markduplicates_MarkDuplicates.php).

The plot will show:

- duplication statistics

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/markduplicates/<sample>`**

- `<sample>.md.cram.metrics`
  - file used by [MultiQC](https://multiqc.info/)
  </details>

#### Sentieon Dedup reports

Sentieon's DNAseq subroutine Dedup produces a metrics report much like the one produced by GATK's MarkDuplicates. The Dedup metrics are imported into MultiQC as custom content and displayed in a table.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/sentieon_dedup/<sample>`**

- `<sample>.dedup.cram.metrics`
  - file used by [MultiQC](https://multiqc.info/).
  </details>

#### samtools stats

[samtools stats](https://www.htslib.org/doc/samtools.html) collects statistics from CRAM files and outputs in a text format.
For further reading and documentation see the [`samtools` manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS).

The plots will show:

- Alignment metrics.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/samtools/<sample>`**

- `<sample>.{sorted,md,recal}.samtools.stats.out`
  - Raw statistics used by `MultiQC`

</details>

#### bcftools stats

[bcftools stats](https://samtools.github.io/bcftools/bcftools.html#stats) produces a statistics text file which is suitable for machine processing and can be plotted using plot-vcfstats.
For further reading and documentation see the [bcftools stats manual](https://samtools.github.io/bcftools/bcftools.html#stats).

Plots will show:

- Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.
- Note: When using [Strelka](https://github.com/Illumina/strelka), there will be no depth distribution plot, as Strelka does not report the INFO/DP field

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/bcftools/`**

- `<sample>.<variantcaller>.bcftools_stats.txt`
  - Raw statistics used by `MultiQC`
  </details>

#### VCFtools

[VCFtools](https://vcftools.github.io/) is a program package designed for working with VCF files. For further reading and documentation see the [VCFtools manual](https://vcftools.github.io/man_latest.html#OUTPUT%20OPTIONS).

Plots will show:

- the summary counts of each type of transition to transversion ratio for each `FILTER` category.
- the transition to transversion ratio as a function of alternative allele count (using only bi-allelic SNPs).
- the transition to transversion ratio as a function of SNP quality threshold (using only bi-allelic SNPs).

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/vcftools/`**

- `<sample>.<variantcaller>.FILTER.summary`
  - Raw statistics used by `MultiQC` with a summary of the number of SNPs and Ts/Tv ratio for each FILTER category
- `<sample>.<variantcaller>.TsTv.count`
  - Raw statistics used by `MultiQC` with the Transition / Transversion ratio as a function of alternative allele count. Only uses bi-allelic SNPs.
- `<sample>.<variantcaller>.TsTv.qual`
  - Raw statistics used by `MultiQC` with Transition / Transversion ratio as a function of SNP quality threshold. Only uses bi-allelic SNPs.
  </details>

#### snpEff reports

[snpeff](http://snpeff.sourceforge.net/) is a genetic variant annotation and effect prediction toolbox.
It annotates and predicts the effects of variants on genes (such as amino acid changes) using multiple databases for annotations. For further reading and documentation see the [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#outputSummary).

The plots will show:

- locations of detected variants in the genome and the number of variants for each location.
- the putative impact of detected variants and the number of variants for each impact.
- the effect of variants at protein level and the number of variants for each effect type.
- the quantity as function of the variant quality score.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/SnpEff/{sample,tumorsample_vs_normalsample}/<variantcaller>/`**

- `<sample>.<variantcaller>_snpEff.csv`
  - Raw statistics used by [MultiQC](http://multiqc.info)
- `<sample>.<variantcaller>_snpEff.html`
  - Statistics to be visualised with a web browser
- `<sample>.<variantcaller>_snpEff.genes.txt`
  - TXT (tab separated) summary counts for variants affecting each transcript and gene
  </details>

#### VEP reports

[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html), based on `Ensembl`, is a tool to determine the effects of all sorts of variants, including SNPs, indels, structural variants, CNVs. For further reading and documentation see the [VEP manual](https://www.ensembl.org/info/docs/tools/vep/index.html)

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/EnsemblVEP/{sample,tumorsamplt_vs_normalsample}/<variantcaller>/`**

- `<sample>.<variantcaller>_VEP.summary.html`
  - Summary of the VEP run to be visualised with a web browser
  </details>

### Reporting

#### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project.
Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.
Results generated by MultiQC collect pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.
  </details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report_<timestamp>.html`, `execution_timeline_<timestamp>.html`, `execution_trace_<timestamp>.txt`, `pipeline_dag_<timestamp>.dot`/`pipeline_dag_<timestamp>.svg` and `manifest_<timestamp>.bco.json`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params_<timestamp>.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Reference files

Contains reference folders generated by the pipeline. These files are only published, if `--save_reference` is set.

<details markdown="1">
<summary>Output files</summary>

- `bwa/`
  - Index corresponding to the [BWA](https://github.com/lh3/bwa) aligner
- `bwamem2/`
  - Index corresponding to the [BWA-mem2](https://github.com/bwa-mem2/bwa-mem2) aligner
- `cnvkit/`
  - Reference files generated by [CNVKit](https://cnvkit.readthedocs.io/en/stable/)
- `dragmap/`
  - Index corresponding to the [DragMap](https://github.com/Illumina/dragmap) aligner
- `dbsnp/`
  - Tabix index generated by [Tabix](http://www.htslib.org/doc/tabix.html) from the given dbsnp file
- `dict/`
  - Sequence dictionary generated by [GATK4 CreateSequenceDictionary](https://gatk.broadinstitute.org/hc/en-us/articles/5358872471963-CreateSequenceDictionary-Picard-) from the given fasta
- `fai/`
  - Fasta index generated with [samtools faidx](http://www.htslib.org/doc/samtools-faidx.html) from the given fasta
- `germline_resource/`
  - Tabix index generated by [Tabix](http://www.htslib.org/doc/tabix.html) from the given gernline resource file
- `intervals/`
  - Bed files in various stages: .bed, .bed.gz, .bed per chromosome, .bed.gz per chromsome
- `known_indels/`
  - Tabix index generated by [Tabix](http://www.htslib.org/doc/tabix.html) from the given known indels file
- `msi/`
  - [MSIsensorPro](https://github.com/xjtu-omics/msisensor-pro) scan of the reference genome to get microsatellites information
- `pon/`
  - Tabix index generated by [Tabix](http://www.htslib.org/doc/tabix.html) from the given panel-of-normals file
  </details>
