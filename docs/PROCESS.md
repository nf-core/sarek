# Workflow processes

Several processes are run within the workflow.
We divide them for the moment into 5 main steps:

## Preprocessing:

- MapReads - Map reads with BWA
- MergeBams - Merge BAMs if multilane samples
- MarkDuplicates - Mark Duplicates with GATK4
- CreateRecalibrationTable - Create Recalibration Table with BaseRecalibrator
- RecalibrateBam - Recalibrate Bam with PrintReads

## Germline Variant Calling:

- CreateIntervalBeds - Create and sort intervals into bed files
- RunHaplotypecaller - Run HaplotypeCaller for Germline Variant Calling (Parallelized processes)
- RunGenotypeGVCFs - Run HaplotypeCaller for Germline Variant Calling (Parallelized processes)
- ConcatVCF - Merge results from paralellized callers
- RunSingleStrelka - Run Strelka for Germline Variant Calling
- RunSingleManta - Run Manta for Single Structural Variant Calling

## Somatic Variant Calling:

- CreateIntervalBeds - Create and sort intervals into bed files
- RunMutect2 - Run MuTect2 for Variant Calling (Parallelized processes)
- RunFreeBayes - Run FreeBayes for Variant Calling (Parallelized processes)
- ConcatVCF - Merge results from paralellized variant callers
- RunStrelka - Run Strelka for Variant Calling
- RunStrelkaBP - Run Strelka Best Practices for Variant Calling
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
- RunVcftools - Run VCFTools on vcf files
- GetVersionAlleleCount - Get version of tools
- GetVersionASCAT - Get version of tools
- GetVersionSnpeff - Get version of tools
- GetVersionVEP - Get version of tools
- GetVersionAll - Get version of tools
- RunMultiQC - Run MultiQC on reports

## Annotation:

- RunSnpeff - Run snpEff for annotation of vcf files
- RunVEP - Run VEP for annotation of vcf files
- CompressVCF - Compress and index vcf files using tabix
