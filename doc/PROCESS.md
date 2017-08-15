# Workflow processes

Several processes are run within the workflow. We divide them for the moment into 3 main steps:

## Preprocessing [Stable]:

- MapReads - Map reads [Stable]
- MergeBams - Merge BAMs if multilane samples [Stable]
- MarkDuplicates - Mark Duplicates [Stable]
- RealignerTargetCreator - Create realignment target intervals [Stable]
- IndelRealigner - Realign Bams as T/N pair [Stable]
- CreateRecalibrationTable - Create Recalibration Table [Stable]
- RecalibrateBam - Recalibreate Bam [Stable]

## Variant Calling:

- RunHaplotypecaller - Run HaplotypeCaller for GermLine Variant Calling (Parrallelized processes) [Stable]
- RunMutect1 - Run MuTect1 for Variant Calling (Parrallelized processes) [Stable]
- RunMutect2 - Run MuTect2 for Variant Calling (Parrallelized processes) [Stable]
- RunFreeBayes - Run FreeBayes for Variant Calling (Parrallelized processes) [Working]
- ConcatVCF - Merge results from HaplotypeCaller, MuTect1 and MuTect2 parrallelized processes [Stable]
- RunStrelka - Run Strelka for Variant Calling [Stable]
- RunManta - Run Manta for Structural Variant Calling [Need more tests]
- RunAlleleCount - Run AlleleCount to prepare for Ascat [Stable]
- RunConvertAlleleCounts - Run convertAlleleCounts to prepare for Ascat [Stable]
- RunAscat - Run Ascat for CNV [Stable]

## Report and QC:

- RunFastQC - Run FastQC for QC on fastq files [Stable]
- MarkDuplicates - Mark Duplicates [Stable]
- RunSamtoolsStats - Run Samtools stats on BAM files [Stable]
- RunBamQC - Run Qualimap bamQC on BAM files [Stable]
- GenerateMultiQCconfig - Generate a config file for MultiQC [Stable]
- RunMultiQC - Run MultiQC for report and QC [Stable]
- RunBcftoolsStats - Run BCFTools stats on vcf files [Working]

## Annotation:

- RunSnpeff - Run snpEff on vcf files [Working]
- RunVEP - Run VEP on vcf files [Working]

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link] [![](images/NGI-final-small.png "NGI")][ngi-link]
[![](doc/images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
