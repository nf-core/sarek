# nf-core/sarek: Output <!-- omit in toc -->

This document describes the output produced by the pipeline.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Preprocessing](#preprocessing)
  - [Map to Reference](#map-to-reference)
    - [BWA mem](#bwa-mem)
  - [Mark Duplicates](#mark-duplicates)
    - [GATK MarkDuplicatesSpark](#gatk-markduplicatesspark)
  - [Base (Quality Score) Recalibration](#base-quality-score-recalibration)
    - [GATK BaseRecalibrator](#gatk-baserecalibrator)
    - [GATK ApplyBQSR](#gatk-applybqsr)
  - [TSV files](#tsv-files)
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
    - [MarkDuplicates reports](#markduplicates-reports)
    - [samtools stats](#samtools-stats)
    - [bcftools stats](#bcftools-stats)
    - [VCFtools](#vcftools)
    - [snpEff reports](#snpeff-reports)
    - [VEP reports](#vep-reports)
  - [Reporting](#reporting)
    - [MultiQC](#multiqc)

## Preprocessing

Sarek preprocesses raw FastQ files or unmapped BAM files, based on [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/).

BAM files with Recalibration tables can also be used as an input to start with the recalibration of said BAM files, for more information see [TSV files output information](#tsv-files)

### Map to Reference

#### BWA mem

[BWA mem](http://bio-bwa.sourceforge.net/) is a software package for mapping low-divergent sequences against a large reference genome.

Such files are intermediate and not kept in the final files delivered to users.

### Mark Duplicates

#### GATK MarkDuplicatesSpark

[GATK MarkDuplicatesSpark](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_spark_transforms_markduplicates_MarkDuplicatesSpark.php) is a Spark implementation of [Picard MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php) and locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.

If the pipeline is run with the option `--no_gatk_spark` then [GATK MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.4.0/picard_sam_markduplicates_MarkDuplicates.php) is used instead.

This directory is the location for the BAM files delivered to users.
Besides the duplicate marked BAM files, the recalibration tables (`*.recal.table`) are also stored, and can be used to create base recalibrated files.

For further reading and documentation see the [data pre-processing workflow from the GATK best practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165).

For all samples:
**Output directory: `results/Preprocessing/[SAMPLE]/DuplicateMarked`**

- `[SAMPLE].md.bam`, `[SAMPLE].md.bai` and `[SAMPLE].recal.table`
  - BAM file and index with Recalibration Table

### Base (Quality Score) Recalibration

#### GATK BaseRecalibrator

[GATK BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php) generates a recalibration table based on various covariates.

Such files are intermediate and not kept in the final files delivered to users.

#### GATK ApplyBQSR

[GATK ApplyBQSR](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_ApplyBQSR.php) recalibrates the base qualities of the input reads based on the recalibration table produced by the [`BaseRecalibrator`](#gatk-baserecalibrator) tool.

This directory is usually empty, it is the location for the final recalibrated BAM files.
Recalibrated BAM files are usually 2-3 times larger than the duplicate marked BAM files.
To re-generate recalibrated BAM file you have to apply the recalibration table delivered to the `DuplicateMarked` directory either within Sarek, or doing this recalibration step yourself.

For further reading and documentation see the [data pre-processing workflow from the GATK best practices](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165).

For all samples:
**Output directory: `results/Preprocessing/[SAMPLE]/Recalibrated`**

- `[SAMPLE].recal.bam` and `[SAMPLE].recal.bam.bai`
  - BAM file and index

### TSV files

The TSV files are auto-generated and can be used by Sarek for further processing and/or variant calling.

For further reading and documentation see the [input documentation](https://github.com/nf-core/sarek/blob/master/docs/input.md).

For all samples:
**Output directory: `results/Preprocessing/TSV`**

- `duplicateMarked.tsv` and `recalibrated.tsv`
  - TSV files to start Sarek from `recalibration` or `variantcalling` steps.
- `duplicateMarked_[SAMPLE].tsv` and `recalibrated_[SAMPLE].tsv`
  - TSV files to start Sarek from `recalibration` or `variantcalling` steps for a specific sample.

> `/!\` Only with [`--sentieon`](usage.md#--sentieon)

For all samples:
**Output directory: `results/Preprocessing/TSV`**

- `recalibrated_sentieon.tsv`
  - TSV files to start Sarek from `variantcalling` step.
- `recalibrated_sentieon_[SAMPLE].tsv`
  - TSV files to start Sarek from `variantcalling` step for a specific sample.

## Variant Calling

All the results regarding Variant Calling are collected in this directory.

Recalibrated BAM files can also be used as an input to start the Variant Calling, for more information see [TSV files output information](#tsv-files)

### SNVs and small indels

#### FreeBayes

[FreeBayes](https://github.com/ekg/freebayes) is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs, indels, MNPs, and complex events smaller than the length of a short-read sequencing alignment..

For further reading and documentation see the [FreeBayes manual](https://github.com/ekg/freebayes/blob/master/README.md#user-manual-and-guide).

For a Tumor/Normal pair only:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/FreeBayes`**

- `FreeBayes_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz` and `FreeBayes_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.tbi`
  - VCF with Tabix index

#### GATK HaplotypeCaller

[GATK HaplotypeCaller](https://github.com/broadinstitute/gatk) calls germline SNPs and indels via local re-assembly of haplotypes.

Germline calls are provided for all samples, to able comparison of both tumor and normal for possible mixup.

For further reading and documentation see the [HaplotypeCaller manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/HaploTypeCaller`**

- `HaplotypeCaller_[SAMPLE].vcf.gz` and `HaplotypeCaller_[SAMPLE].vcf.gz.tbi`
  - VCF with Tabix index

#### GATK GenotypeGVCFs

[GATK GenotypeGVCFs](https://github.com/broadinstitute/gatk) performs joint genotyping on one or more samples pre-called with HaplotypeCaller.

Germline calls are provided for all samples, to able comparison of both tumor and normal for possible mixup.

For further reading and documentation see the [GenotypeGVCFs manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/HaplotypeCallerGVCF`**

- `HaplotypeCaller_[SAMPLE].g.vcf.gz` and `HaplotypeCaller_[SAMPLE].g.vcf.gz.tbi`
  - VCF with Tabix index

#### GATK Mutect2

[GATK Mutect2](https://github.com/broadinstitute/gatk) calls somatic SNVs and indels via local assembly of haplotypes.

For further reading and documentation see the [Mutect2 manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php).
It is recommended to have panel of normals [PON](https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2) for this version of Mutect2 using at least 40 normal samples, and you can add your PON file to get filtered somatic calls.

For a Tumor/Normal pair only:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Mutect2`**

Files created:

- `Mutect2_unfiltered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz` and `Mutect2_unfiltered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.tbi`
  - unfiltered (raw) Mutect2 calls VCF with Tabix index
- `Mutect2_filtered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz` and `Mutect2_filtered_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.tbi`
  - filtered Mutect2 calls VCF with Tabix index: these entries has a PASS filter, you can get these when supplying a panel of normals using the `--pon` option
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.stats`
  - a stats file generated during calling raw variants (needed for filtering)
- `[TUMORSAMPLE]_contamination.table`
  - a text file exported when panel-of-normals provided about sample contamination

#### samtools mpileup

[samtools mpileup](https://www.htslib.org/doc/samtools.html) generate pileup for a BAM file.

For further reading and documentation see the [samtools manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/mpileup`**

- `[SAMPLE].pileup.gz`
  - The pileup format is a text-based format for summarizing the base calls of aligned reads to a reference sequence. Alignment records are grouped by sample (SM) identifiers in @RG header lines.

#### Strelka2

[Strelka2](https://github.com/Illumina/strelka) is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs.

For further reading and documentation see the [Strelka2 user guide](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/Strelka`**

- `Strelka_Sample_genome.vcf.gz` and `Strelka_Sample_genome.vcf.gz.tbi`
  - VCF with Tabix index
- `Strelka_Sample_variants.vcf.gz` and `Strelka_Sample_variants.vcf.gz.tbi`
  - VCF with Tabix index

For a Tumor/Normal pair:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Strelka`**

- `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz` and `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz.tbi`
  - VCF with Tabix index
- `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz` and `Strelka_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz.tbi`
  - VCF with Tabix index

Using [Strelka Best Practices](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic-configuration-example) with the `candidateSmallIndels` from `Manta`:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Strelka`**

- `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz` and `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_indels.vcf.gz.tbi`
  - VCF with Tabix index
- `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz` and `StrelkaBP_[TUMORSAMPLE]_vs_[NORMALSAMPLE]_somatic_snvs.vcf.gz.tbi`
  - VCF with Tabix index

#### Sentieon DNAseq

> `/!\` Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon DNAseq](https://www.sentieon.com/products/#dnaseq) implements the same mathematics used in the Broad Institute's BWA-GATK HaplotypeCaller 3.3-4.1 Best Practices Workflow pipeline.

For further reading and documentation see the [Sentieon DNAseq user guide](https://support.sentieon.com/manual/DNAseq_usage/dnaseq/).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/SentieonDNAseq`**

- `DNAseq_Sample.vcf.gz` and `DNAseq_Sample.vcf.gz.tbi`
  - VCF with Tabix index

#### Sentieon DNAscope

> `/!\` Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon DNAscope](https://www.sentieon.com/products) calls SNPs and small indels.

For further reading and documentation see the [Sentieon DNAscope user guide](https://support.sentieon.com/manual/DNAscope_usage/dnascope/).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/SentieonDNAscope`**

- `DNAscope_Sample.vcf.gz` and `DNAscope_Sample.vcf.gz.tbi`
  - VCF with Tabix index

#### Sentieon TNscope

> `/!\` Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon TNscope](https://www.sentieon.com/products/#tnscope) calls SNPs and small indels on an Tumor/Normal pair.

For further reading and documentation see the [Sentieon TNscope user guide](https://support.sentieon.com/manual/TNscope_usage/tnscope/).

For a Tumor/Normal pair:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/SentieonTNscope`**

- `TNscope_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz` and `TNscope_[TUMORSAMPLE]_vs_[NORMALSAMPLE].vcf.gz.tbi`
  - VCF with Tabix index

### Structural Variants

#### Manta

[Manta](https://github.com/Illumina/manta) calls structural variants (SVs) and indels from mapped paired-end sequencing reads.
It is optimized for analysis of germline variation in small sets of individuals and somatic variation in tumor/normal sample pairs.
`Manta` provides a candidate list for small indels also that can be fed to `Strelka` following [Strelka Best Practices](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#somatic-configuration-example).

For further reading and documentation see the [Manta user guide](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/Manta`**

- `Manta_[SAMPLE].candidateSmallIndels.vcf.gz` and `Manta_[SAMPLE].candidateSmallIndels.vcf.gz.tbi`
  - VCF with Tabix index
- `Manta_[SAMPLE].candidateSV.vcf.gz` and `Manta_[SAMPLE].candidateSV.vcf.gz.tbi`
  - VCF with Tabix index

For Normal sample only:

- `Manta_[NORMALSAMPLE].diploidSV.vcf.gz` and `Manta_[NORMALSAMPLE].diploidSV.vcf.gz.tbi`
  - VCF with Tabix index

For a Tumor sample only:

- `Manta_[TUMORSAMPLE].tumorSV.vcf.gz` and `Manta_[TUMORSAMPLE].tumorSV.vcf.gz.tbi`
  - VCF with Tabix index

For a Tumor/Normal pair only:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/Manta`**

- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSmallIndels.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSmallIndels.vcf.gz.tbi`
  - VCF with Tabix index
- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSV.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].candidateSV.vcf.gz.tbi`
  - VCF with Tabix index
- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].diploidSV.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].diploidSV.vcf.gz.tbi`
  - VCF with Tabix index
- `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].somaticSV.vcf.gz` and `Manta_[TUMORSAMPLE]_vs_[NORMALSAMPLE].somaticSV.vcf.gz.tbi`
  - VCF with Tabix index

#### TIDDIT

[TIDDIT](https://github.com/SciLifeLab/TIDDIT) identifies intra and inter-chromosomal translocations, deletions, tandem-duplications and inversions.

Germline calls are provided for all samples, to able comparison of both tumor and normal for possible mixup.
Low quality calls are removed internally, to simplify processing of variant calls but they are saved by Sarek.

For further reading and documentation see the [TIDDIT manual](https://github.com/SciLifeLab/TIDDIT/blob/master/README.md).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/TIDDIT`**

- `TIDDIT_[SAMPLE].vcf.gz` and `TIDDIT_[SAMPLE].vcf.gz.tbi`
  - VCF with Tabix index
- `TIDDIT_[SAMPLE].signals.tab`
  - tab file describing coverage across the genome, binned per 50 bp
- `TIDDIT_[SAMPLE].ploidy.tab`
  - tab file describing the estimated ploidy and coverage across each contig
- `TIDDIT_[SAMPLE].old.vcf`
  - VCF including the low qualiy calls
- `TIDDIT_[SAMPLE].wig`
  - wiggle file containing coverage across the genome, binned per 50 bp
- `TIDDIT_[SAMPLE].gc.wig`
  - wiggle file containing fraction of gc content, binned per 50 bp

#### Sentieon DNAscope SV

> `/!\` Only with [`--sentieon`](usage.md#--sentieon)

[Sentieon DNAscope](https://www.sentieon.com/products) can perform structural variant calling in addition to calling SNPs and small indels.

For further reading and documentation see the [Sentieon DNAscope user guide](https://support.sentieon.com/manual/DNAscope_usage/dnascope/).

For all samples:
**Output directory: `results/VariantCalling/[SAMPLE]/SentieonDNAscope`**

- `DNAscope_SV_Sample.vcf.gz` and `DNAscope_SV_Sample.vcf.gz.tbi`
  - VCF with Tabix index

### Sample heterogeneity, ploidy and CNVs

#### ConvertAlleleCounts

[ConvertAlleleCounts](https://github.com/nf-core/sarek/blob/master/bin/convertAlleleCounts.r) is a R-script for converting output from AlleleCount to BAF and LogR values.

For a Tumor/Normal pair only:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/ASCAT`**

- `[TUMORSAMPLE].BAF` and `[NORMALSAMPLE].BAF`
  - file with beta allele frequencies
- `[TUMORSAMPLE].LogR` and `[NORMALSAMPLE].LogR`
  - file with total copy number on a logarithmic scale

#### ASCAT

[ASCAT](https://github.com/Crick-CancerGenomics/ascat) is a method to derive copy number profiles of tumor cells, accounting for normal cell admixture and tumor aneuploidy.
ASCAT infers tumor purity and ploidy and calculates whole-genome allele-specific copy number profiles.

For further reading and documentation see [the Sarek documentation about ASCAT](https://github.com/nf-core/sarek/blob/master/docs/ascat.md) or the [ASCAT manual](https://www.crick.ac.uk/research/labs/peter-van-loo/software).

For a Tumor/Normal pair only:
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

#### Control-FREEC

[Control-FREEC](https://github.com/BoevaLab/FREEC) is a tool for detection of copy-number changes and allelic imbalances (including LOH) using deep-sequencing data.
Control-FREEC automatically computes, normalizes, segments copy number and beta allele frequency profiles, then calls copy number alterations and LOH.
And also detects subclonal gains and losses and evaluate the likeliest average ploidy of the sample.

For further reading and documentation see the [Control-FREEC manual](http://boevalab.com/FREEC/tutorial.html).

For a Tumor/Normal pair only:
**Output directory: `results/VariantCalling/[TUMOR_vs_NORMAL]/ControlFREEC`**

- `[TUMORSAMPLE]_vs_[NORMALSAMPLE].config.txt`
  - Configuration file used to run Control-FREEC
- `[TUMORSAMPLE].pileup.gz_CNVs` and `[TUMORSAMPLE].pileup.gz_normal_CNVs`
  - file with coordinates of predicted copy number alterations
- `[TUMORSAMPLE].pileup.gz_ratio.txt` and `[TUMORSAMPLE].pileup.gz_normal_ratio.txt`
  - file with ratios and predicted copy number alterations for each window
- `[TUMORSAMPLE].pileup.gz_BAF.txt` and `[NORMALSAMPLE].pileup.gz_BAF.txt`
  - file with beta allele frequencies for each possibly heterozygous SNP position

### MSI status

[Microsatellite instability](https://en.wikipedia.org/wiki/Microsatellite_instability)
is a genetic condition associated to deficienceies in the
mismatch repair (MMR) system which causes a tendency to accumulate a high
number of mutations (SNVs and indels).

#### MSIsensor

[MSIsensor](https://github.com/ding-lab/msisensor) is a tool to detect the MSI
status of a tumor scaning the length of the microsatellite regions. An altered
distribution of  microsatellite length is associated to a missed replication
slippage which would be corrected under normal mismatch repair (MMR) conditions. It requires
a normal sample for each tumour to differentiate the somatic and germline
cases.

For further reading see the [MSIsensor paper](https://www.ncbi.nlm.nih.gov/pubmed/24371154).

For a Tumor/Normal pair only:
**Output directory: `results/VariantCalling/[TUMORSAMPLE]_vs_[NORMALSAMPLE]/MSIsensor`**

- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]`_msisensor
  - MSI score output, contains information about the number of somatic sites.
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]`_msisensor_dis
  - The normal and tumor length distribution for each microsatellite position.
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]`_msisensor_germline
  - somatic sites detected
- `[TUMORSAMPLE]_vs_[NORMALSAMPLE]`_msisensor_somatic
  - germ line sites detected

## Variant annotation

This directory contains results from the final annotation steps: two software are used for annotation, [snpEff](http://snpeff.sourceforge.net/) and [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html).
Only a subset of the VCF files are annotated, and only variants that have a PASS filter.
FreeBayes results are not annotated in the moment yet as we are lacking a decent somatic filter.
For HaplotypeCaller the germline variations are annotated for both the tumor and the normal sample.

### snpEff

[snpeff](http://snpeff.sourceforge.net/) is a genetic variant annotation and effect prediction toolbox.
It annotates and predicts the effects of variants on genes (such as amino acid changes) using multiple databases for annotations.
The generated VCF header contains the software version and the used command line.

For further reading and documentation see the [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#outputSummary)

For all samples:
**Output directory: `results/Annotation/[SAMPLE]/snpEff`**

- `VariantCaller_Sample_snpEff.ann.vcf.gz` and `VariantCaller_Sample_snpEff.ann.vcf.gz.tbi`
  - VCF with Tabix index

### VEP

[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html), based on Ensembl, is a tools to determine the effects of all sorts of variants, including SNPs, indels, structural variants, CNVs.
The generated VCF header contains the software version, also the version numbers for additional databases like Clinvar or dbSNP used in the "VEP" line.
The format of the [consequence annotations](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html) is also in the VCF header describing the INFO field.
In the moment it contains:

- Consequence: impact of the variation, if there is any
- Codons: the codon change, i.e. cGt/cAt
- Amino_acids: change in amino acids, i.e. R/H if there is any
- Gene: ENSEMBL gene name
- SYMBOL: gene symbol
- Feature: actual transcript name
- EXON: affected exon
- PolyPhen: prediction based on [PolyPhen](http://genetics.bwh.harvard.edu/pph2/)
- SIFT: prediction by [SIFT](http://sift.bii.a-star.edu.sg/)
- Protein_position: Relative position of amino acid in protein
- BIOTYPE: Biotype of transcript or regulatory feature

For further reading and documentation see the [VEP manual](https://www.ensembl.org/info/docs/tools/vep/index.html)

For all samples:
**Output directory: `results/Annotation/[SAMPLE]/VEP`**

- `VariantCaller_Sample_VEP.ann.vcf.gz` and `VariantCaller_Sample_VEP.ann.vcf.gz.tbi`
  - VCF with Tabix index

## QC and reporting

### QC

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads.
It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C).
You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

For all samples:
**Output directory: `results/Reports/[SAMPLE]/fastqc`**

- `sample_R1_XXX_fastqc.html` and `sample_R2_XXX_fastqc.html`
  - FastQC report, containing quality metrics for each pair of the raw fastq files
- `sample_R1_XXX_fastqc.zip` and `sample_R2_XXX_fastqc.zip`
  - zip file containing the FastQC reports, tab-delimited data files and plot images

#### bamQC

[Qualimap bamqc](http://qualimap.bioinfo.cipf.es/) reports information for the evaluation of the quality of the provided alignment data. In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced.

Plot will show:

- Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

For all samples:
**Output directory: `results/Reports/[SAMPLE]/bamQC`**

- `VariantCaller_[SAMPLE].bcf.tools.stats.out`
  - RAW statistics used by MultiQC

For more information about how to use Qualimap bamqc reports, see [Qualimap bamqc manual](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#id7)

#### MarkDuplicates reports

[[GATK MarkDuplicatesSpark](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_spark_transforms_markduplicates_MarkDuplicatesSpark.php), Spark implementation of [Picard MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/picard_sam_markduplicates_MarkDuplicates.php) locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.

If the pipeline is run with the option `--no_gatk_spark` then [GATK MarkDuplicates](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.4.0/picard_sam_markduplicates_MarkDuplicates.php) is used instead.

Collecting duplicate metrics slows down performance.
To disable them use `--skipQC MarkDuplicates`.

Duplicates can arise during sample preparation _e.g._ library construction using PCR.
Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.
These duplication artifacts are referred to as optical duplicates.

For all samples:
**Output directory: `results/Reports/[SAMPLE]/MarkDuplicates`**

- `[SAMPLE].bam.metrics`
  - RAW statistics used by MultiQC

For further reading and documentation see the [MarkDuplicates manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_markduplicates_MarkDuplicates.php).

#### samtools stats

[samtools stats](https://www.htslib.org/doc/samtools.html) collects statistics from BAM files and outputs in a text format.
Plots will show:

- Alignment metrics.

For all samples:
**Output directory: `results/Reports/[SAMPLE]/SamToolsStats`**

- `[SAMPLE].bam.samtools.stats.out`
  - RAW statistics used by MultiQC

For further reading and documentation see the [samtools manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS)

#### bcftools stats

[bcftools](https://samtools.github.io/bcftools/) is a program for variant calling and manipulating files in the Variant Call Format.
Plot will show:

- Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

For all samples:
**Output directory: `results/Reports/[SAMPLE]/BCFToolsStats`**

- `VariantCaller_[SAMPLE].bcf.tools.stats.out`
  - RAW statistics used by MultiQC

For further reading and documentation see the [bcftools stats manual](https://samtools.github.io/bcftools/bcftools.html#stats)

#### VCFtools

[VCFtools](https://vcftools.github.io/) is a program package designed for working with VCF files.
Plots will show:

- the summary counts of each type of transition to transversion ratio for each FILTER category.
- the transition to transversion ratio as a function of alternative allele count (using only bi-allelic SNPs).
- the transition to transversion ratio as a function of SNP quality threshold (using only bi-allelic SNPs).

For all samples:
**Output directory: `results/Reports/[SAMPLE]/VCFTools`**

- `VariantCaller_[SAMPLE].FILTER.summary`
  - RAW statistics used by MultiQC
- `VariantCaller_[SAMPLE].TsTv.count`
  - RAW statistics used by MultiQC
- `VariantCaller_[SAMPLE].TsTv.qual`
  - RAW statistics used by MultiQC

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
  - RAW statistics used by MultiQC
- `VariantCaller_Sample_snpEff.html`
  - Statistics to be visualised with a web browser
- `VariantCaller_Sample_snpEff.genes.txt`
  - TXT (tab separated) summary counts for variants affecting each transcript and gene

For further reading and documentation see the [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#outputSummary)

#### VEP reports

[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html), based on Ensembl, is a tools to determine the effects of all sorts of variants, including SNPs, indels, structural variants, CNVs.

For all samples:
**Output directory: `results/Reports/[SAMPLE]/VEP`**

- `VariantCaller_Sample_VEP.summary.html`
  - Summary of the VEP run to be visualised with a web browser

For further reading and documentation see the [VEP manual](https://www.ensembl.org/info/docs/tools/vep/index.html)

### Reporting

#### MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project.
Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

For the whole Sarek run:
**Output directory: `results/Reports/MultiQC`**

- `multiqc_report.html`
  - MultiQC report - a standalone HTML file that can be viewed in your web browser
- `multiqc_data/`
  - Directory containing parsed statistics from the different tools used in the pipeline

For further reading and documentation see the [MultiQC website](http://multiqc.info)
