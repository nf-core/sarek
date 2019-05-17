# nf-core/sarek: Output

This document describes the output produced by the pipeline.
Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

1. **Preprocessing** _(based on [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/))_
    * Map reads to Reference
        * [BWA](#BWA)
    * Mark Duplicates
        * [GATK MarkDuplicates](#MarkDuplicates)
    * Base (Quality Score) Recalibration
        * [GATK BaseRecalibrator](#BaseRecalibrator)
        * [GATK GatherBQSRReports](#GatherBQSRReports)
        * [GATK ApplyBQSR](#ApplyBQSR)
2. **Variant calling**
    * SNVs and small indels
        * [Freebayes](#Freebayes)
        * [GATK HaplotypeCaller](#HaplotypeCaller)
        * [GATK GenotypeGVCFs](#GenotypeGVCFs)
        * [GATK MuTect2](#MuTect2)
        * [Strelka2](#Strelka2)
    * Structural variants
        * [Manta](#Manta)
    * Sample heterogeneity, ploidy and CNVs
        * [alleleCounter](#alleleCounter)
        * [ConvertAlleleCounts](#ConvertAlleleCounts)
        * [ASCAT](#ASCAT)
        * [samtools mpileup](#mpileup)
        * [Control-FREEC](#Control-FREEC)
3. **Annotation**
    * Variant annotation
        * [snpEff](#snpEff)
        * [VEP (Variant Effect Predictor)](#VEP)
4. **QC and Reporting**
    * QC
        * [FastQC](#FastQC)
        * [Qualimap bamqc](#bamQC)
        * [GATK MarkDuplicates](#MarkDuplicates-reports)
        * [samtools stats](#Samtools-stats)
        * [bcftools stats](#BCFtools)
        * [VCFtools](#VCFtools)
        * [snpeff](#snpEff-reports)
    * Reporting
        * [MultiQC](#MultiQC)

## Preprocessing
Sarek preprocessing raw FastQ files or unmapped BAM files is based on [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/)).
BAM files with Recalibration tables can also be used as an input to start with the recalibration of said BAM files.

### BWA
[BWA](http://bio-bwa.sourceforge.net/)

### MarkDuplicates
[GATK MarkDuplicates](https://github.com/broadinstitute/gatk)

### BaseRecalibrator
[GATK BaseRecalibrator](https://github.com/broadinstitute/gatk)

### GatherBQSRReports
[GATK GatherBQSRReports](https://github.com/broadinstitute/gatk)

### ApplyBQSR
[GATK ApplyBQSR](https://github.com/broadinstitute/gatk)

## Variant Calling

### Freebayes
[Freebayes](https://github.com/ekg/freebayes)

### HaplotypeCaller
[GATK HaplotypeCaller](https://github.com/broadinstitute/gatk)

### GenotypeGVCFs
[GATK GenotypeGVCFs](https://github.com/broadinstitute/gatk)

### MuTect2
[MuTect2](https://github.com/broadinstitute/gatk)

### Strelka2
[Strelka2](https://github.com/Illumina/strelka)

### Manta
[Manta](https://github.com/Illumina/manta)

### alleleCounter

### ConvertAlleleCounts

### ASCAT

### samtools mpileup

### Control-FREEC

## Annotation

### snpEff
[snpeff](http://snpeff.sourceforge.net/)

### VEP
[VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html)

## QC and reports

### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads.
It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C).
You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `results/Reports/[SAMPLE]/fastqc`**

* `sample_R1_XXX_fastqc.html` and `sample_R2_XXX_fastqc.html`
  * FastQC report, containing quality metrics for each pair of the raw fastq files
* `sample_R1_XXX_fastqc.zip` and `sample_R2_XXX_fastqc.zip`
  * zip file containing the FastQC reports, tab-delimited data files and plot images

### bamQC
[Qualimap bamqc](http://qualimap.bioinfo.cipf.es/) reports information for the evaluation of the quality of the provided alignment data. In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced.

Plot will show:
* Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

**Output directory: `results/Reports/[SAMPLE]/bamQC`**
* `VariantCaller_Sample.bcf.tools.stats.out`
  * RAW statistics used by MultiQC

For more information about how to use Qualimap bamqc reports, see [Qualimap bamqc manual](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#id7)

### MarkDuplicates -reports
[GATK MarkDuplicates](https://github.com/broadinstitute/gatk) locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.
Duplicates can arise during sample preparation e.g.
library construction using PCR.
Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.
These duplication artifacts are referred to as optical duplicates.

Plot will show:
* Duplication metrics

**Output directory: `results/Reports/[SAMPLE]/MarkDuplicates`**
* `Sample.bam.metrics`
  * RAW statistics used by MultiQC

For more information about how to use MarkDuplicates reports, see [MarkDuplicates manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_markduplicates_MarkDuplicates.php)

### samtools stats
[samtools stats](https://www.htslib.org/doc/samtools.html) collects statistics from BAM files and outputs in a text format.
Plots will show:
* Alignment metrics.

**Output directory: `results/Reports/[SAMPLE]/SamToolsStats`**
* `Sample.bam.samtools.stats.out`
  * RAW statistics used by MultiQC

For more information about how to use samtools stats reports, see [samtools stats manual](http://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS)

### bcftools stats
[bcftools stats](https://samtools.github.io/bcftools/) is a program for variant calling and manipulating files in the Variant Call Format.
Plot will show:
* Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

**Output directory: `results/Reports/[SAMPLE]/BCFToolsStats`**
* `VariantCaller_Sample.bcf.tools.stats.out`
  * RAW statistics used by MultiQC

For more information about how to use bcftools stats reports, see [bcftools stats manual](https://samtools.github.io/bcftools/bcftools.html#stats)

### VCFtools
[VCFtools](https://vcftools.github.io/) is a program package designed for working with VCF files.
Plots will show:
* the summary counts of each type of transition to transversion ratio for each FILTER category.
* the transition to transversion ratio as a function of alternative allele count (using only bi-allelic SNPs).
* the transition to transversion ratio as a function of SNP quality threshold (using only bi-allelic SNPs).

**Output directory: `results/Reports/[SAMPLE]/VCFTools`**
* `VariantCaller_Sample.FILTER.summary`
  * RAW statistics used by MultiQC
* `VariantCaller_Sample.TsTv.count`
  * RAW statistics used by MultiQC
* `VariantCaller_Sample.TsTv.qual`
  * RAW statistics used by MultiQC

For more information about how to use VCFtools reports, see [VCFtools manual](https://vcftools.github.io/man_latest.html#OUTPUT%20OPTIONS)

### snpEff reports
[snpeff](http://snpeff.sourceforge.net/) is a genetic variant annotation and effect prediction toolbox.
It annotates and predicts the effects of variants on genes (such as amino acid changes).
Plots will shows :
* locations of detected variants in the genome and the number of variants for each location.
* the putative impact of detected variants and the number of variants for each impact.
* the effect of variants at protein level and the number of variants for each effect type.
* the quantity as function of the variant quality score.

**Output directory: `results/Reports/[SAMPLE]/snpEff`**
* `VariantCaller_Sample_snpEff.csv`
  * RAW statistics used by MultiQC
* `VariantCaller_Sample_snpEff.html`
  * Statistics to be visualised with a web browser
* `VariantCaller_Sample_snpEff.txt`
  * TXT (tab separated) summary counts for variants affecting each transcript and gene

For more information about how to use snpEff reports, see [snpEff manual](http://snpeff.sourceforge.net/SnpEff_manual.html#outputSummary)

### MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project.
Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/Reports/MultiQC`**

* `multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
