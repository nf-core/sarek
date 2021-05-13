# nf-core/sarek: Output <!-- omit in toc -->

## Introduction <!-- omit in toc -->

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [Map to Reference](#map-to-reference)
    - [bwa](#bwa)
    - [BWA-mem2](#bwa-mem2)
  - [Mark Duplicates](#mark-duplicates)
    - [GATK MarkDuplicates](#gatk-markduplicates)
  - [Base (Quality Score) Recalibration](#base-quality-score-recalibration)
    - [GATK BaseRecalibrator](#gatk-baserecalibrator)
    - [GATK ApplyBQSR](#gatk-applybqsr)
  - [TSV files](#tsv-files)
  - [TSV files with `--skip_markduplicates`](#tsv-files-with---skip_markduplicates)
  - [TSV files with `--sentieon`](#tsv-files-with---sentieon)
- [Variant Calling](#variant-calling)
  - [SNVs and small indels](#snvs-and-small-indels)
    - [FreeBayes](#freebayes)
    - [GATK HaplotypeCaller](#gatk-haplotypecaller)
    - [GATK GenotypeGVCFs](#gatk-genotypegvcfs)
    - [GATK Mutect2](#gatk-mutect2)
    - [samtools mpileup](#samtools-mpileup)
    - [Strelka2](#strelka2)
    - [Sentieon DNAseq](#sentieon-dnaseq)
    - [Sentieon DNAscope](#sentieon-dnascope)
    - [Sentieon TNscope](#sentieon-tnscope)
  - [Structural Variants](#structural-variants)
    - [Manta](#manta)
    - [TIDDIT](#tiddit)
    - [Sentieon DNAscope SV](#sentieon-dnascope-sv)
  - [Sample heterogeneity, ploidy and CNVs](#sample-heterogeneity-ploidy-and-cnvs)
    - [ConvertAlleleCounts](#convertallelecounts)
    - [ASCAT](#ascat)
    - [Control-FREEC](#control-freec)
  - [MSI status](#msi-status)
    - [MSIsensor](#msisensor)
- [Variant annotation](#variant-annotation)
  - [snpEff](#snpeff)
  - [VEP](#vep)
- [QC and reporting](#qc-and-reporting)
  - [QC](#qc)
    - [FastQC](#fastqc)
    - [bamQC](#bamqc)
    - [GATK MarkDuplicates reports](#gatk-markduplicates-reports)
    - [samtools stats](#samtools-stats)
    - [bcftools stats](#bcftools-stats)
    - [VCFtools](#vcftools)
    - [snpEff reports](#snpeff-reports)
    - [VEP reports](#vep-reports)
  - [Reporting](#reporting)
    - [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

## Preprocessing

`Sarek` pre-processes raw `FASTQ` files or `unmapped BAM` files, based on [GATK best practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows).

### Map to Reference

#### bwa

[bwa](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome.

Such files are intermediate and not kept in the final files delivered to users.

#### BWA-mem2

[BWA-mem2](https://github.com/bwa-mem2/bwa-mem2) is a software package for mapping low-divergent sequences against a large reference genome.

Such files are intermediate and not kept in the final files delivered to users.

### Mark Duplicates

#### GATK MarkDuplicates

By default, `Sarek` will use [GATK MarkDuplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/360042912511-MarkDuplicatesSpark), `Spark` implementation of [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360042477492-MarkDuplicates-Picard), which locates and tags duplicate reads in a `BAM` or `SAM` file, where duplicate reads are defined as originating from a single fragment of DNA.

Specify `--no_gatk_spark` to use `GATK MarkDuplicates` instead.

This directory is the location for the `BAM` files delivered to users.
Besides the `duplicates-marked BAM` files, the recalibration tables (`*.recal.table`) are also stored, and can be used to create `recalibrated BAM` files.

For all samples:

**Output directory: `results/Preprocessing/[SAMPLE]/DuplicatesMarked`**

- `[SAMPLE].md.bam` and `[SAMPLE].md.bai`
  - `BAM` file and index

For further reading and documentation see the [data pre-processing for variant discovery from the GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).

### Base (Quality Score) Recalibration

#### GATK BaseRecalibrator

[GATK BaseRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360042477672-BaseRecalibrator) generates a recalibration table based on various co-variates.

For all samples:

**Output directory: `results/Preprocessing/[SAMPLE]/DuplicatesMarked`**

- `[SAMPLE].recal.table`
  - Recalibration table associated to the `duplicates-marked BAM` file.

#### GATK ApplyBQSR

[GATK ApplyBQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360042476852-ApplyBQSR) recalibrates the base qualities of the input reads based on the recalibration table produced by the [GATK BaseRecalibrator](#gatk-baserecalibrator) tool.

This directory is the location for the final `recalibrated BAM` files.
`Recalibrated BAM` files are usually 2-3 times larger than the `duplicates-marked BAM` files.
To re-generate `recalibrated BAM` file you have to apply the recalibration table delivered to the `DuplicatesMarked\` folder either using `Sarek` ( [`--step recalibrate`](usage.md#step-recalibrate) ) , or doing this recalibration yourself.

For all samples:

**Output directory: `results/Preprocessing/[SAMPLE]/Recalibrated`**

- `[SAMPLE].recal.bam` and `[SAMPLE].recal.bam.bai`
  - `BAM` file and index

For further reading and documentation see the [data pre-processing for variant discovery from the GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).

### TSV files

The `TSV` files are auto-generated and can be used by `Sarek` for further processing and/or variant calling.

For further reading and documentation see the [`--input`](usage.md#--input) section in the usage documentation.

For all samples:

**Output directory: `results/Preprocessing/TSV`**

- `duplicates_marked_no_table.tsv`, `duplicates_marked.tsv` and `recalibrated.tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps.
- `duplicates_marked_no_table_[SAMPLE].tsv`, `duplicates_marked_[SAMPLE].tsv` and `recalibrated_[SAMPLE].tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps for a specific sample.

### TSV files with `--skip_markduplicates`

> **WARNING** Only with [`--skip_markduplicates`](usage.md#--skip_markduplicates)

For all samples:

**Output directory: `results/Preprocessing/TSV`**

- `mapped.tsv`, `mapped_no_duplicates_marked.tsv` and `recalibrated.tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps.
- `mapped_[SAMPLE].tsv`, `mapped_no_duplicates_marked_[SAMPLE].tsv` and `recalibrated_[SAMPLE].tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps for a specific sample.

### TSV files with `--sentieon`

> **WARNING** Only with [`--sentieon`](usage.md#--sentieon)

For all samples:

**Output directory: `results/Preprocessing/TSV`**

- `sentieon_deduped.tsv` and `recalibrated_sentieon.tsv`
  - `TSV` files to start `Sarek` from `variantcalling` step.
- `sentieon_deduped_[SAMPLE].tsv` and `recalibrated_sentieon_[SAMPLE].tsv`
  - `TSV` files to start `Sarek` from `variantcalling` step for a specific sample.

## Variant Calling

All the results regarding Variant Calling are collected in this directory.
If some results from a variant caller do not appear here, please check out the [`--tools`](usage.md#--tools) section in the usage documentation.

`Recalibrated BAM` files can used as an input to start the Variant Calling.

### SNVs and small indels

#### FreeBayes

[FreeBayes](https://github.com/ekg/freebayes) is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs, indels, MNPs, and complex events smaller than the length of a short-read sequencing alignment.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/FreeBayes`**

- `FreeBayes_[SAMPLE].vcf.gz` and `FreeBayes_[SAMPLE].vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [FreeBayes manual](https://github.com/ekg/freebayes/blob/master/README.md#user-manual-and-guide).

#### GATK HaplotypeCaller

[GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller) calls germline SNPs and indels via local re-assembly of haplotypes.

Germline calls are provided for all samples, to enable comparison of both, tumor and normal, for possible mixup.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/HaploTypeCaller`**

- `HaplotypeCaller_[SAMPLE].vcf.gz` and `HaplotypeCaller_[SAMPLE].vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [HaplotypeCaller manual](https://gatk.broadinstitute.org/hc/en-us/articles/360042913231-HaplotypeCaller).

#### GATK GenotypeGVCFs

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

[samtools mpileup](https://www.htslib.org/doc/samtools.html) generates pileup of a `BAM` file.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/mpileup`**

- `[SAMPLE].pileup.gz`
  - The pileup format is a text-based format for summarizing the base calls of aligned reads to a reference sequence. Alignment records are grouped by sample (`SM`) identifiers in `@RG` header lines.

For further reading and documentation see the [samtools manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS).

#### Strelka2

[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/Strelka`**

- `Strelka_Sample_genome.vcf.gz` and `Strelka_Sample_genome.vcf.gz.tbi`
  - `VCF` with Tabix index
- `Strelka_Sample_variants.vcf.gz` and `Strelka_Sample_variants.vcf.gz.tbi`
  - `VCF` with Tabix index

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Strelka`**

- `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz` and `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz.tbi`
  - `VCF` with Tabix index
- `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz` and `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz.tbi`
  - `VCF` with Tabix index

Using [Strelka Best Practices](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic-configuration-example) with the `candidateSmallIndels` from `Manta`:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Strelka`**

- `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz` and `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz.tbi`
  - `VCF` with Tabix index
- `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz` and `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [Strelka2 user guide](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md).

#### Sentieon DNAseq

> **WARNING** Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon DNAseq](https://www.sentieon.com/products/#dnaseq) implements the same mathematics used in the Broad Institute's BWA-GATK HaplotypeCaller 3.3-4.1 Best Practices Workflow pipeline.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/SentieonDNAseq`**

- `DNAseq_Sample.vcf.gz` and `DNAseq_Sample.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [Sentieon DNAseq user guide](https://support.sentieon.com/manual/DNAseq_usage/dnaseq/).

#### Sentieon DNAscope

> **WARNING** Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon DNAscope](https://www.sentieon.com/products) calls SNPs and small indels.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/SentieonDNAscope`**

- `DNAscope_Sample.vcf.gz` and `DNAscope_Sample.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [Sentieon DNAscope user guide](https://support.sentieon.com/manual/DNAscope_usage/dnascope/).

#### Sentieon TNscope

> **WARNING** Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon TNscope](https://www.sentieon.com/products/#tnscope) calls SNPs and small indels on an Tumor/Normal pair.

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/SentieonTNscope`**

- `TNscope_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz` and `TNscope_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [Sentieon TNscope user guide](https://support.sentieon.com/manual/TNscope_usage/tnscope/).

### Structural Variants

#### Manta

[Manta](https://github.com/Illumina/manta) calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
`Manta` provides a candidate list for small indels that can be fed to `Strelka` following [Strelka Best Practices](https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic-configuration-example).

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/Manta`**

- `Manta_[SAMPLE].candidateSmallIndels.vcf.gz` and `Manta_[SAMPLE].candidateSmallIndels.vcf.gz.tbi`
  - `VCF` with Tabix index
- `Manta_[SAMPLE].candidateSV.vcf.gz` and `Manta_[SAMPLE].candidateSV.vcf.gz.tbi`
  - `VCF` with Tabix index

For Normal sample only:

- `Manta_[NORMALSAMPLE].diploidSV.vcf.gz` and `Manta_[NORMALSAMPLE].diploidSV.vcf.gz.tbi`
  - `VCF` with Tabix index

For a Tumor sample only:

- `Manta_[TUMORSAMPLE].tumorSV.vcf.gz` and `Manta_[TUMORSAMPLE].tumorSV.vcf.gz.tbi`
  - `VCF` with Tabix index

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Manta`**

- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSmallIndels.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSmallIndels.vcf.gz.tbi`
  - `VCF` with Tabix index
- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSV.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSV.vcf.gz.tbi`
  - `VCF` with Tabix index
- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].diploidSV.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].diploidSV.vcf.gz.tbi`
  - `VCF` with Tabix index
- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].somaticSV.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].somaticSV.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [Manta user guide](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md).

#### TIDDIT

[TIDDIT](https://github.com/SciLifeLab/TIDDIT) identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions.

Germline calls are provided for all samples, to enable comparison of both, tumor and normal, for possible mixup.
Low quality calls are removed internally, to simplify processing of variant calls but they are saved by `Sarek`.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/TIDDIT`**

- `TIDDIT_[SAMPLE].vcf.gz` and `TIDDIT_[SAMPLE].vcf.gz.tbi`
  - `VCF` with Tabix index
- `TIDDIT_[SAMPLE].signals.tab`
  - tab file describing coverage across the genome, binned per 50 bp
- `TIDDIT_[SAMPLE].ploidy.tab`
  - tab file describing the estimated ploidy and coverage across each contig
- `TIDDIT_[SAMPLE].old.vcf`
  - `VCF` including the low qualiy calls
- `TIDDIT_[SAMPLE].wig`
  - wiggle file containing coverage across the genome, binned per 50 bp
- `TIDDIT_[SAMPLE].gc.wig`
  - wiggle file containing fraction of gc content, binned per 50 bp

For further reading and documentation see the [TIDDIT manual](https://github.com/SciLifeLab/TIDDIT/blob/master/README.md).

#### Sentieon DNAscope SV

> **WARNING** Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon DNAscope](https://www.sentieon.com/products) can perform structural variant calling in addition to calling SNPs and small indels.

For all samples:

**Output directory: `results/VariantCalling/[SAMPLE]/SentieonDNAscope`**

- `DNAscope_SV_Sample.vcf.gz` and `DNAscope_SV_Sample.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [Sentieon DNAscope user guide](https://support.sentieon.com/manual/DNAscope_usage/dnascope/).

### Sample heterogeneity, ploidy and CNVs

#### ConvertAlleleCounts

Running ASCAT on NGS data requires that the `BAM` files are converted into BAF and LogR values.
This can be done using the software [AlleleCount](https://github.com/cancerit/alleleCount) followed by the provided [ConvertAlleleCounts](https://github.com/nf-core/sarek/blob/master/bin/convertAlleleCounts.r) R-script.

For a Tumor/Normal pair:

**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/ASCAT`**

- `[TUMORSAMPLE].BAF` and `[NORMALSAMPLE].BAF`
  - file with beta allele frequencies
- `[TUMORSAMPLE].LogR` and `[NORMALSAMPLE].LogR`
  - file with total copy number on a logarithmic scale

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

- *chr*: chromosome number
- *startpos*: start position of the segment
- *endpos*: end position of the segment
- *nMajor*: number of copies of one of the allels (for example the chromosome inherited from the father)
- *nMinor*: number of copies of the other allele (for example the chromosome inherited of the mother)

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

### MSI status

[Microsatellite instability](https://en.wikipedia.org/wiki/Microsatellite_instability) is a genetic condition associated to deficiencies in the mismatch repair (MMR) system which causes a tendency to accumulate a high number of mutations (SNVs and indels).
An altered distribution of microsatellite length is associated to a missed replication slippage which would be corrected under normal MMR conditions.

#### MSIsensor

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

- *Consequence*: impact of the variation, if there is any
- *Codons*: the codon change, i.e. cGt/cAt
- *Amino_acids*: change in amino acids, i.e. R/H if there is any
- *Gene*: ENSEMBL gene name
- *SYMBOL*: gene symbol
- *Feature*: actual transcript name
- *EXON*: affected exon
- *PolyPhen*: prediction based on [PolyPhen](http://genetics.bwh.harvard.edu/pph2/)
- *SIFT*: prediction by [SIFT](http://sift.bii.a-star.edu.sg/)
- *Protein_position*: Relative position of amino acid in protein
- *BIOTYPE*: Biotype of transcript or regulatory feature

For all samples:

**Output directory: `results/Annotation/[SAMPLE]/VEP`**

- `VariantCaller_Sample_VEP.ann.vcf.gz` and `VariantCaller_Sample_VEP.ann.vcf.gz.tbi`
  - `VCF` with Tabix index

For further reading and documentation see the [VEP manual](https://www.ensembl.org/info/docs/tools/vep/index.html)

## QC and reporting

### QC

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads.
It provides information about the quality score distribution across your reads, per base sequence content (`%A/T/G/C`), adapter contamination and overrepresented sequences.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/fastqc`**

- `sample_R1_XXX_fastqc.html` and `sample_R2_XXX_fastqc.html`
  - `FastQC` report containing quality metrics for your untrimmed raw `FASTQ` files
- `sample_R1_XXX_fastqc.zip` and `sample_R2_XXX_fastqc.zip`
  - Zip archive containing the FastQC report, tab-delimited data file and plot images

> **NB:** The `FastQC` plots displayed in the `MultiQC` report shows _untrimmed_ reads.
> They may contain adapter sequence and potentially regions with low quality.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

#### bamQC

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

**Output files:**

- `multiqc/`  
  - `multiqc_report.html`
    - Standalone HTML file that can be viewed in your web browser
  - `multiqc_data/`
    - Directory containing parsed statistics from the different tools used in the pipeline
  - `multiqc_plots/`
    - Directory containing static images from the report in various formats

For more information about how to use `MultiQC` reports, see [https://multiqc.info](https://multiqc.info).

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  - Documentation for interpretation of results in HTML format: `results_description.html`.
