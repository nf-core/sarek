# nf-core/sarek: Output <!-- omit in toc -->

## Introduction <!-- omit in toc -->

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

The directories listed below will be created in the results directory after the pipeline has finished.
All paths are relative to the top-level results directory.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [Prepare input files](#prepare-input)
    - [Trim adapters](#trimming)
    - [Split fastq files](#split)
  - [Map to Reference](#map-to-reference)
    - [BWA](#bwa)
    - [BWA-mem2](#bwa-mem2)
    - [DragMap](#dragmap)
  - [Duplicate Marking](#duplicate-marking)
    - [GATK MarkDuplicates (Spark)](#gatk-markduplicates)
  - [Base Quality Score Recalibration](#base-quality-score-recalibration)
    - [GATK BaseRecalibrator (Spark)](#gatk-baserecalibrator)
    - [GATK ApplyBQSR (Spark)](#gatk-applybqsr)
  - [CSV files](#csv-files)
- [Variant Calling](#variant-calling)
  - [SNVs and small indels](#snvs-and-small-indels)
    - [DeepVariant](#deepvariant)
    - [FreeBayes](#freebayes)
    - [GATK HaplotypeCaller](#gatk-haplotypecaller)
    - [GATK Mutect2](#gatk-mutect2)
    - [samtools mpileup](#samtools-mpileup)
    - [Strelka2](#strelka2)
  - [Structural Variants](#structural-variants)
    - [Manta](#manta)
    - [TIDDIT](#tiddit)
  - [Sample heterogeneity, ploidy and CNVs](#sample-heterogeneity-ploidy-and-cnvs)
    - [ASCAT](#ascat)
    - [Control-FREEC](#control-freec)
    - [CNVKit](#cnvkit)
  - [MSI status](#msi-status)
    - [MSIsensorPro](#msisensorpro)
- [Variant annotation](#variant-annotation)
  - [snpEff](#snpeff)
  - [VEP](#vep)
- [QC and reporting](#qc-and-reporting)
  - [QC](#qc)
    - [FastQC](#fastqc)
    - [FastP](#fastp)
    - [GATK MarkDuplicates reports](#gatk-markduplicates-reports)
    - [Mosdepth](#mosdepth)
    - [samtools stats](#samtools-stats)
    - [bcftools stats](#bcftools-stats)
    - [VCFtools](#vcftools)
    - [snpEff reports](#snpeff-reports)
    - [VEP reports](#vep-reports)
  - [Reporting](#reporting)
    - [MultiQC](#multiqc)
  - [Pipeline information](#pipeline-information)
- [Reference files](#reference-files)

## Preprocessing

`Sarek` pre-processes raw `FASTQ` files or `unmapped BAM` files, based on [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

### Preparation of input files (FastQ or (u)BAM)

[FastP](https://github.com/OpenGene/fastp) is a tool designed to provide all-in-one preprocessing for FastQ files and as such is used for trimming and splitting. By default, these files are not published. However, if publishing is enabled, please be aware that these files are only published once, meaning if trimming and splitting is enabled, then the resulting files will be sharded fastq files with trimmed reads. If only one of them is enabled then the files contain either trimmed or split reads respectively.

#### Trim adapters

[FastP](https://github.com/OpenGene/fastp) supports global trimming, which means trim all reads in the front or the tail. This function is useful since sometimes you want to drop some cycles of a sequencing run. In the current implementation in `sarek`
`--detect_adapter_for_pe` is set by default which enables auto-detection of adapter sequences. For more information on how to fine-tune adapter trimming, take a look into the parameter docs.

The resulting files are intermediate and by default not kept in the final files delivered to users. Set `--save_trimmed`to enable publishing of the files in `{outdir}/preprocessing/<sample>/fastp/<sample_lane_{1,2}.fastp.fastq.gz>`

#### Split fastq files

[FastP](https://github.com/OpenGene/fastp) supports splitting of FastQ into multiple files which enables parallel alignment. To enable splitting, the number of reads per output can be specified. For more information, take a look into the parameter `--split_fastq`in the parameter docs .

These files are intermediate and by default not kept in the final files delivered to users. Set `--save_split` to enable publishing of these files to `{outdir}/preprocessing/<sample>/fastp/<sample_lane_{1,2}.fastp.fastq.gz>`

### Map to Reference

#### BWA

[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome. The aligned reads are then `coordinate` sorted (or `name` sorted if MarkDuplicates Spark is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html)

These files are intermediate and by default not kept in the final files delivered to users. Set `--save_bam_mapped` to enable publishing.

#### BWA-mem2

[BWA-mem2](https://github.com/bwa-mem2/bwa-mem2) is a software package for mapping low-divergent sequences against a large reference genome.The aligned reads are then `coordinate` sorted (or `name` sorted if MarkDuplicates Spark is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html)

These files are intermediate and by default not kept in the final files delivered to users.Set `--save_bam_mapped` to enable publishing.

#### DragMap

[DragMap](https://github.com/Illumina/dragmap) is an open-source software implementation of the DRAGEN mapper, which the Illumina team created so that we would have an open-source way to produce the same results as their proprietary DRAGEN hardware. The aligned reads are then `coordinate` sorted (or `name` sorted if MarkDuplicates Spark is used for duplicate marking) with [samtools](https://www.htslib.org/doc/samtools.html)

These files are intermediate and by default not kept in the final files delivered to users. Set `--save_bam_mapped` to enable publishing.

For all mappers and samples:

**Output directory: `{outdir}/preprocessing/<sample>/mapped`**

- `<sample>.bam` and `<sample>.bam.bai`
  - `BAM` file and index

### Mark Duplicates

During duplicate marking read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artifactual processes are identified. These are considered to be non-independent observations, so all but a single read pair within each set of duplicates, causing the marked pairs to be ignored by default during the variant discovery process

For further reading and documentation see the [data pre-processing for variant discovery from the GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).

#### GATK MarkDuplicates (Spark)

By default, `Sarek` will use [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/5358880192027-MarkDuplicates-Picard-).

Specify `--use_gatk_spark markduplicates` to use [`GATK MarkDuplicatesSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) instead, the `Spark` implementation. The resulting files are converted to `CRAM` either with [samtools](https://www.htslib.org/doc/samtools.html) when `GATK MarkDuplicates` is used or implicitly by `GATK MarkDuplicatesSpark`

The resulting `CRAM` files are delivered to the users.

For all samples:

**Output directory: `{outdir}/preprocessing/<sample>/markduplicates`**

- `<sample>.md.cram` and `<sample>.md.cram.crai`
  - `CRAM` file and index
- if `--save_output_as_bam`:
  - `<sample>.md.bam` and `<sample>.md.bam.bai`

### Base (Quality Score) Recalibration

During Base Quality Score Recalibration systematic errors in the base quality scores are corrected by applying machine learning to detect and correct for them. This is important, for evaluating the correct call of a variant during the variant discovery process. However, this is not needed for all combinations of tools in sarek. Notably, this should be turned of when having UMI tagged reads or using DragMap (see [here](https://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode)) as mapper.

For further reading and documentation see the [technical documentation by GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-).

#### GATK BaseRecalibrator (Spark)

[GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360042477672-BaseRecalibrator) generates a recalibration table based on various co-variates.

Specify `--use_gatk_spark baserecalibrator` to use [`GATK BaseRecalibratorSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358896138011-BaseRecalibrator) instead, the respective `Spark` implementation.

For all samples:

**Output directory: `results/preprocessing/<sample>/recal_table`**

- `<sample>.recal.table`
  - Recalibration table associated to the `duplicates-marked CRAM` file.

#### GATK ApplyBQSR (Spark)

[GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/5358826654875-ApplyBQSR) recalibrates the base qualities of the input reads based on the recalibration table produced by the [GATK BaseRecalibrator](#gatk-baserecalibrator) tool.

Specify `--use_gatk_spark baserecalibrator` to use [`GATK ApplyBQSRSpark`](https://gatk.broadinstitute.org/hc/en-us/articles/5358898266011-ApplyBQSRSpark-BETA-) instead, the respective `Spark` implementation.

The resulting `recalibrated CRAM` files are delivered to the user. `Recalibrated CRAM` files are usually 2-3 times larger than the `duplicate-marked CRAM` files.
To re-generate `recalibrated CRAM` files you have to apply the recalibration table delivered to the `recal_table/` folder either using `Sarek` ( [`--step recalibrate`](usage.md#step-recalibrate) ) , or doing this recalibration yourself.

For all samples:

**Output directory: `results/Preprocessing/<sample>/recalibrated`**

- `<sample>.recal.cram` and `<sample>.recal.cram.crai`
  - `CRAM` file and index
- if `--save_output_as_bam`:
  - `<sample>.recal.bam` and `<sample>.recal.bam.bai`

### CSV files

The `CSV` files are auto-generated and can be used by `Sarek` for further processing and/or variant calling.

For further reading and documentation see the [`--input`](usage.md#--input) section in the usage documentation.

For all samples:

**Output directory: `results/Preprocessing/TSV`**

- `duplicates_marked_no_table.tsv`, `duplicates_marked.tsv` and `recalibrated.tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps.
- `duplicates_marked_no_table_[SAMPLE].tsv`, `duplicates_marked_[SAMPLE].tsv` and `recalibrated_[SAMPLE].tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps for a specific sample.

> **WARNING** Only with [`--skip_markduplicates`](usage.md#--skip_markduplicates)

For all samples:

**Output directory: `results/Preprocessing/TSV`**

- `mapped.tsv`, `mapped_no_duplicates_marked.tsv` and `recalibrated.tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps.
- `mapped_[SAMPLE].tsv`, `mapped_no_duplicates_marked_[SAMPLE].tsv` and `recalibrated_[SAMPLE].tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps for a specific sample.

## Variant Calling

The results regarding variant calling are collected in `{outdir}/variantcalling/`.
If some results from a variant caller do not appear here, please check out the [`--tools`]section in the parameter [documentation](https://nf-co.re/sarek/3.0.0/parameters).

`(Recalibrated) CRAM` files can used as an input to start the Variant Calling.

### SNVs and small indels

For single nucleotide variants (SNVs) and small indels, multiple tools are available for normal, tumor-only, and tumor-normal paired data. For a list of the appropriate tool(s) for the data and sequencing type at hand, please check [here](usage.md#Which tool can be used for for which data type?)

#### DeepVariant

[DeepVariant](https://github.com/google/deepvariant) is a deep learning-based variant caller that takes aligned reads, produces pileup image tensors from them, classifies each tensor using a convolutional neural network, and finally reports the results in a standard VCF or gVCF file.

For normal samples:

**Output directory: `results/variantcalling/<sample>/deepvariant`**

- `<sample>.vcf.gz` and `<sample>.vcf.gz.tbi`
  - `VCF` with Tabix index
- `<sample>.g.vcf.gz` and `<sample>.g.vcf.gz.tbi`
  - `.g.VCF` with Tabix index

#### FreeBayes

[FreeBayes](https://github.com/ekg/freebayes) is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs, indels, MNPs, and complex events smaller than the length of a short-read sequencing alignment. For further reading and documentation see the [FreeBayes manual](https://github.com/ekg/freebayes/blob/master/README.md#user-manual-and-guide).

For all samples:

**Output directory: `results/variantcalling/{sample,normalsample_vs_tumorsample}/freebayes`**

- `<sample>.vcf.gz` and `<sample>.vcf.gz.tbi`
  - `VCF` with Tabix index

#### GATK HaplotypeCaller

[GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller) calls germline SNPs and indels via local re-assembly of haplotypes.

Germline calls are provided for all samples, to enable comparison of both, tumor and normal, for possible mixup.

For normal samples:

**Output directory: `results/VariantCalling/[SAMPLE]/HaploTypeCaller`**

- `HaplotypeCaller_[SAMPLE].vcf.gz` and `HaplotypeCaller_[SAMPLE].vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [HaplotypeCaller manual](https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller).

##### GATK GenotypeGVCFs

[GATK GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360042914991-GenotypeGVCFs) performs joint genotyping on one or more samples pre-called with HaplotypeCaller.

Germline calls are provided for all samples, to enable comparison of both, tumor and normal, for possible mixup.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/HaplotypeCallerGVCF`**

- `HaplotypeCaller_[SAMPLE].g.vcf.gz` and `HaplotypeCaller_[SAMPLE].g.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [GenotypeGVCFs manual](https://gatk.broadinstitute.org/hc/en-us/articles/360042914991-GenotypeGVCFs).

#### GATK Mutect2

[GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360042477952-Mutect2) calls somatic SNVs and indels via local assembly of haplotypes.

For further reading and documentation see the [Mutect2 manual](https://gatk.broadinstitute.org/hc/en-us/articles/360042477952-Mutect2).
It is recommended to have [panel of normals (PON)](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON) for this version of `GATK Mutect2` using at least 40 normal samples.
Additionally, you can add your `PON` file to get filtered somatic calls.

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Mutect2`**

Files created:

- `Mutect2_unfiltered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz` and `Mutect2_unfiltered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.tbi`
  - unfiltered (raw) Mutect2 calls `VCF` with Tabix index
- `Mutect2_filtered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz` and `Mutect2_filtered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.tbi`
  - filtered Mutect2 calls `VCF` with Tabix index: these entries have a `PASS` filter, you can get these when supplying a panel of normals using the `--pon` option
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.stats`
  - a stats file generated during calling of raw variants (needed for filtering)
- `[TUMORSAMPLE]_contamination.table`
  - a text file exported when panel-of-normals about sample contamination are provided

#### samtools mpileup

[samtools mpileup](https://www.htslib.org/doc/samtools.html) generates pileup of a `CRAM` file.
For further reading and documentation see the [samtools manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS).

For all samples:

**Output directory: `results/variantcalling/<sample>/mpileup`**

- `<sample>.pileup.gz`
  - The pileup format is a text-based format for summarizing the base calls of aligned reads to a reference sequence. Alignment records are grouped by sample (`SM`) identifiers in `@RG` header lines.

#### Strelka2

[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. For further reading and documentation see the [Strelka2 user guide](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md). If [Strelka2](https://github.com/Illumina/strelka) is used for somatic variant calling and [Manta](https://github.com/Illumina/manta) is also specified in tools, the output candidate indels from [Manta](https://github.com/Illumina/manta) are used according to [Strelka Best Practices](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic-configuration-example).

For single samples (normal-only or tumor-only):

**Output directory: `results/variantcalling/<sample>/strelka`**

- `<sample>.genome.vcf.gz` and `<sample>.genome.vcf.gz.tbi`
  - `VCF` with Tabix index
- `<sample>.variants.vcf.gz` and `<sample>.variants.vcf.gz.tbi`
  - `VCF` with Tabix index

For a Tumor/Normal pair:

**Output directory: `results/variantcalling/<tumorsample_vs_normalsample>/strelka`**

- `<tumorsample_vs_normalsample>.somatic_indels.vcf.gz` and `<tumorsample_vs_normalsample>.somatic_indels.vcf.gz.tbi`
  - `VCF` with Tabix index
- `<tumorsample_vs_normalsample>.somatic_snvs.vcf.gz` and `<tumorsample_vs_normalsample>.somatic_snvs.vcf.gz.tbi`
  - `VCF` with Tabix index

### Structural Variants

#### Manta

[Manta](https://github.com/Illumina/manta) calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
`Manta` provides a candidate list for small indels that can be fed to `Strelka` following [Strelka Best Practices](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic-configuration-example). For further reading and documentation see the [Manta user guide](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md).

For normal samples:

**Output directory: `results/variantcalling/<sample>/manta`**

- `<sample>.diploid_sv.vcf.gz` and `<sample>.diploid_sv.vcf.gz.tbi`
  - `VCF` with Tabix index

For a Tumor-only samples:

- `<sample>.tumor_sv.vcf.gz` and `<sample>.tumor_sv.vcf.gz.tbi`
  - `VCF` with Tabix index

For a Tumor/Normal pair:

**Output directory: `results/variantcalling/<tumorsample_vs_normalsample>/manta`**

- `<tumorsample_vs_normalsample>.diploid_sv.vcf.gz` and `<tumorsample_vs_normalsample>.diploid_sv.vcf.gz.tbi`
  - `VCF` with Tabix index
- `<tumorsample_vs_normalsample>.somatic_sv.vcf.gz` and `<tumorsample_vs_normalsample>.somatic_sv.vcf.gz.tbi`
  - `VCF` with Tabix index

#### TIDDIT

[TIDDIT](https://github.com/SciLifeLab/TIDDIT) identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions.

Germline calls are provided for all samples, to enable comparison of both, tumor and normal, for possible mixup.
Low quality calls are removed internally, to simplify processing of variant calls but they are saved by `Sarek`.
For further reading and documentation see the [TIDDIT manual](https://github.com/SciLifeLab/TIDDIT/blob/master/README.md).

For normal-only and tumor-only samples:

**Output directory: `results/variantcalling/<sample>/tiddit`**

- `<sample>.vcf.gz` and `<sample>.vcf.gz.tbi`
  - `VCF` with Tabix index
- `<sample>.ploidy.tab`
  - tab file describing the estimated ploidy and coverage across each contig

For tumor/normal paired samples:

**Output directory: `results/variantcalling/<tumorsample_vs_normal_sample>/tiddit`**

- `<tumorsample_vs_normal_sample>.normal.vcf.gz` and `<tumorsample_vs_normal_sample>.normal.vcf.gz.tbi`
  - `VCF` with Tabix index
- `<tumorsample_vs_normal_sample>.tumor.vcf.gz` and `<tumorsample_vs_normal_sample>.tumor.vcf.gz.tbi`
  - `VCF` with Tabix index
- `<tumorsample_vs_normal_sample>_sv_merge.vcf`
  - `VCF`
- `<tumorsample_vs_normal_sample>.ploidy.tab`
  - tab file describing the estimated ploidy and coverage across each contig

### Sample heterogeneity, ploidy and CNVs

#### ASCAT

[ASCAT](https://github.com/Crick-CancerGenomics/ascat) is a software for performing allele-specific copy number analysis of tumor samples and for estimating tumor ploidy and purity (normal contamination).
It infers tumor purity and ploidy and calculates whole-genome allele-specific copy number profiles.
`ASCAT` is written in `R` and available here: [github.com/Crick-CancerGenomics/ascat](https://github.com/Crick-CancerGenomics/ascat).
The `ASCAT` process gives several images as output, described in detail in this [book chapter](http://www.ncbi.nlm.nih.gov/pubmed/22130873).

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/ASCAT`**

- `[TUMORSAMPLE].aberrationreliability.png`
  - Image with information about aberration reliability
- `[TUMORSAMPLE].ASCATprofile.png`
  - Image with information about ASCAT profile
- `[TUMORSAMPLE].ASPCF.png`
  - Image with information about ASPCF
- `[TUMORSAMPLE].rawprofile.png`
  - Image with information about raw profile
- `[TUMORSAMPLE].sunrise.png`
  - Image with information about sunrise
- `[TUMORSAMPLE].tumour.png`
  - Image with information about tumor
- `[TUMORSAMPLE].cnvs.txt`
  - file with information about CNVS
- `[TUMORSAMPLE].LogR.PCFed.txt`
  - file with information about LogR
- `[TUMORSAMPLE].purityploidy.txt`
  - file with information about purity ploidy

The text file `[TUMORSAMPLE].cnvs.txt` countains predictions about copy number state for all the segments.
The output is a tab delimited text file with the following columns:

- _chr_: chromosome number
- _startpos_: start position of the segment
- _endpos_: end position of the segment
- _nMajor_: number of copies of one of the allels (for example the chromosome inherited from the father)
- _nMinor_: number of copies of the other allele (for example the chromosome inherited of the mother)

The file `[TUMORSAMPLE].cnvs.txt` contains all segments predicted by ASCAT, both those with normal copy number (nMinor = 1 and nMajor =1) and those corresponding to copy number aberrations.

For further reading and documentation see the [ASCAT manual](https://www.crick.ac.uk/research/labs/peter-van-loo/software).

#### Control-FREEC

[Control-FREEC](https://github.com/BoevaLab/FREEC) is a tool for detection of copy-number changes and allelic imbalances (including loss of heterozygoity (LOH)) using deep-sequencing data.
`Control-FREEC` automatically computes, normalizes, segments copy number and beta allele frequency profiles, then calls copy number alterations and LOH.
And also detects subclonal gains and losses and evaluate the most likely average ploidy of the sample.

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/ControlFREEC`**

- `[TUMORSAMPLE]_vs_[NORMALSAMPLE].config.txt`
  - Configuration file used to run Control-FREEC
- `[TUMORSAMPLE].pileup.gz_CNVs` and `[TUMORSAMPLE].pileup.gz_normal_CNVs`
  - file with coordinates of predicted copy number alterations
- `[TUMORSAMPLE].pileup.gz_ratio.txt` and `[TUMORSAMPLE].pileup.gz_normal_ratio.txt`
  - file with ratios and predicted copy number alterations for each window
- `[TUMORSAMPLE].pileup.gz_BAF.txt` and `[NORMALSAMPLE].pileup.gz_BAF.txt`
  - file with beta allele frequencies for each possibly heterozygous SNP position

For further reading and documentation see the [Control-FREEC manual](http://boevalab.com/FREEC/tutorial.html).

#### CNVKit

### MSI status

[Microsatellite instability](https://en.wikipedia.org/wiki/Microsatellite_instability) is a genetic condition associated to deficiencies in the mismatch repair (MMR) system which causes a tendency to accumulate a high number of mutations (SNVs and indels).
An altered distribution of microsatellite length is associated to a missed replication slippage which would be corrected under normal MMR conditions.

#### MSIsensorPro

[MSIsensor](https://github.com/ding-lab/msisensor) is a tool to detect the MSI status of a tumor scanning the length of the microsatellite regions.
It requires a normal sample for each tumour to differentiate the somatic and germline cases.

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMORSAMPLE]_vs_[NORMALSAMPLE]/MSIsensor`**

- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]_msisensor`
  - MSI score output, contains information about the number of somatic sites.
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]_msisensor_dis`
  - The normal and tumor length distribution for each microsatellite position.
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]_msisensor_germline`
  - Somatic sites detected.
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]_msisensor_somatic`
  - Germline sites detected.

For further reading see the [MSIsensor paper](https://www.ncbi.nlm.nih.gov/pubmed/24371154).

## Variant annotation

This directory contains results from the final annotation steps: two tools are used for annotation, [snpEff](http://snpeff.sourceforge.net/) and [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html).
Only a subset of the `VCF` files are annotated, and only variants that have a `PASS` filter.
Currently, `FreeBayes` results are not annotated as we are lacking a decent somatic filter.

### snpEff

[snpeff](http://snpeff.sourceforge.net/) is a genetic variant annotation and effect prediction toolbox.
It annotates and predicts the effects of variants on genes (such as amino acid changes) using multiple databases for annotations.
The generated `VCF` header contains the software version and the used command line.

For all samples:

**Output directory: `results/Annotation/[SAMPLE]/snpEff`**

- `VariantCaller_Sample_snpEff.ann.vcf.gz` and `VariantCaller_Sample_snpEff.ann.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#outputSummary)

### VEP

[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html), based on `Ensembl`, is a tool to determine the effects of all sorts of variants, including SNPs, indels, structural variants, CNVs.
The generated `VCF` header contains the software version, also the version numbers for additional databases like `Clinvar` or `dbSNP` used in the `VEP` line.
The format of the [consequence annotations](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html) is also in the `VCF` header describing the `INFO` field.
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

For all samples:

**Output directory: `results/Annotation/[SAMPLE]/VEP`**

- `VariantCaller_Sample_VEP.ann.vcf.gz` and `VariantCaller_Sample_VEP.ann.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [VEP manual](https://www.ensembl.org/info/docs/tools/vep/index.html)

## QC and reporting

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

For all samples:

**Output directory: `results/Reports/[SAMPLE]/fastqc`**

- `sample_R1_XXX_fastqc.html` and `sample_R2_XXX_fastqc.html`
  - `FastQC` report containing quality metrics for your untrimmed raw `FASTQ` files
- `sample_R1_XXX_fastqc.zip` and `sample_R2_XXX_fastqc.zip`
  - Zip archive containing the FastQC report, tab-delimited data file and plot images

> **NB:** The `FastQC` plots displayed in the `MultiQC` report shows _untrimmed_ reads.
> They may contain adapter sequence and potentially regions with low quality.

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

#### FastP

#### Mosdepth

[Qualimap bamqc](http://qualimap.bioinfo.cipf.es/) reports information for the evaluation of the quality of the provided alignment data.
In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced.

Plot will show:

- Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/bamQC`**

- `VariantCaller_[SAMPLE].bcf.tools.stats.out`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [Qualimap bamqc manual](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#id7)

#### GATK MarkDuplicates reports

More information in the [GATK MarkDuplicates section](#gatk-markduplicates)

Duplicates can arise during sample preparation _e.g._ library construction using PCR.
Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.
These duplication artifacts are referred to as optical duplicates.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/MarkDuplicates`**

- `[SAMPLE].bam.metrics`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [MarkDuplicates manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_markduplicates_MarkDuplicates.php).

#### samtools stats

[samtools stats](https://www.htslib.org/doc/samtools.html) collects statistics from `BAM` files and outputs in a text format.

Plots will show:

- Alignment metrics.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/SamToolsStats`**

- `[SAMPLE].bam.samtools.stats.out`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [`samtools` manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS)

#### bcftools stats

[bcftools](https://samtools.github.io/bcftools/) is a program for variant calling and manipulating `VCF` files.

Plot will show:

- Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/BCFToolsStats`**

- `VariantCaller_[SAMPLE].bcf.tools.stats.out`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [bcftools stats manual](https://samtools.github.io/bcftools/bcftools.html#stats)

#### VCFtools

[VCFtools](https://vcftools.github.io/) is a program package designed for working with `VCF` files.

Plots will show:

- the summary counts of each type of transition to transversion ratio for each `FILTER` category.
- the transition to transversion ratio as a function of alternative allele count (using only bi-allelic SNPs).
- the transition to transversion ratio as a function of SNP quality threshold (using only bi-allelic SNPs).

For all samples:

**Output directory: `results/Reports/[SAMPLE]/VCFTools`**

- `VariantCaller_[SAMPLE].FILTER.summary`
  - Raw statistics used by `MultiQC`
- `VariantCaller_[SAMPLE].TsTv.count`
  - Raw statistics used by `MultiQC`
- `VariantCaller_[SAMPLE].TsTv.qual`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [VCFtools manual](https://vcftools.github.io/man_latest.html#OUTPUT%20OPTIONS)

#### snpEff reports

[snpeff](http://snpeff.sourceforge.net/) is a genetic variant annotation and effect prediction toolbox.
It annotates and predicts the effects of variants on genes (such as amino acid changes) using multiple databases for annotations.

Plots will shows :

- locations of detected variants in the genome and the number of variants for each location.
- the putative impact of detected variants and the number of variants for each impact.
- the effect of variants at protein level and the number of variants for each effect type.
- the quantity as function of the variant quality score.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/snpEff`**

- `VariantCaller_Sample_snpEff.csv`
  - Raw statistics used by `MultiQC`
- `VariantCaller_Sample_snpEff.html`
  - Statistics to be visualised with a web browser
- `VariantCaller_Sample_snpEff.genes.txt`
  - TXT (tab separated) summary counts for variants affecting each transcript and gene

For further reading and documentation see the [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#outputSummary)

#### VEP reports

[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html), based on `Ensembl`, is a tools to determine the effects of all sorts of variants, including SNPs, indels, structural variants, CNVs.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/VEP`**

- `VariantCaller_Sample_VEP.summary.html`
  - Summary of the VEP run to be visualised with a web browser

For further reading and documentation see the [VEP manual](https://www.ensembl.org/info/docs/tools/vep/index.html)

### Reporting

#### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project.
Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the `MultiQC` output for future traceability.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

For more information about how to use `MultiQC` reports, see [https://multiqc.info](https://multiqc.info).

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

## Reference files
