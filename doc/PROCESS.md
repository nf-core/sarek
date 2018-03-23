# Workflow processes

Several processes are run within the workflow.
We divide them for the moment into 5 main steps:

## Preprocessing:

- MapReads - Map reads with BWA
- MergeBams - Merge BAMs if multilane samples
- MarkDuplicates - Mark Duplicates with Picard
- RealignerTargetCreator - Create realignment target intervals
- IndelRealigner - Realign BAMs as T/N pair
- CreateRecalibrationTable - Create Recalibration Table with BaseRecalibrator
- RecalibrateBam - Recalibrate Bam with PrintReads

## Germline Variant Calling:

- CreateIntervalBeds - Create and sort intervals into bed files
- RunHaplotypecaller - Run HaplotypeCaller for GermLine Variant Calling (Parallelized processes)
- RunGenotypeGVCFs - Run HaplotypeCaller for GermLine Variant Calling (Parallelized processes)
- ConcatVCF - Merge results from HaplotypeCaller
- RunSingleStrelka - Run Strelka for Germline Variant Calling
- RunSingleManta - Run Manta for Single Structural Variant Calling

## Somatic Variant Calling:

- CreateIntervalBeds - Create and sort intervals into bed files
- RunMutect1 - Run MuTect1 for Variant Calling (Parallelized processes)
- RunMutect2 - Run MuTect2 for Variant Calling (Parallelized processes)
- RunFreeBayes - Run FreeBayes for Variant Calling (Parallelized processes)
- ConcatVCF - Merge results from Freebayes, MuTect1 and MuTect2
- RunStrelka - Run Strelka for Variant Calling
- RunManta - Run Manta for Structural Variant Calling
- RunSingleManta - Run Manta for Single Structural Variant Calling
- RunAlleleCount - Run AlleleCount to prepare for ASCAT
- RunConvertAlleleCounts - Run convertAlleleCounts to prepare for ASCAT
- RunAscat - Run ASCAT for CNV

## Report and QC:

- RunFastQC - Run FastQC for QC on fastq files
- RunSamtoolsStats - Run Samtools stats on recalibrated BAM files
- RunBamQC - Run qualimap BamQC on recalibrated BAM files
- RunBcftoolsStats - Run BCFTools stats on vcf files

## Annotation:

- RunSnpeff - Run snpEff for annotation of vcf files
- RunVEP - Run VEP for annotation of vcf files

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
