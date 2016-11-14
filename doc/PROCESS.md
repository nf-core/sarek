# Workflow processes
Several processes are run within the workflow. We divide them for the moment into 2 main steps:

## Preprocessing:
 - MapReads - Map reads
 - MergeBams - Merge BAMs if multilane samples
 - MarkDuplicates - Mark Duplicates
 - CreateIntervals - Create Intervals
 - RealignBams - Realign Bams as T/N pair
 - CreateRecalibrationTable - Create Recalibration Table
 - RecalibrateBam - Recalibreate Bam

## Variant Calling:
 - RunHaplotypecaller - Run HaplotypeCaller for GermLine Variant Calling (Parrallelized processes)
 - RunMutect1 - Run MuTect1 for Variant Calling (Parrallelized processes)
 - RunMutect2 - Run MuTect2 for Variant Calling (Parrallelized processes)
 - RunVardict - Run VarDict for Variant Calling (Parrallelized processes)
 - ConcatVCF - Merge results from HaplotypeCaller, MuTect1, MuTect2 and VarDict parrallelized processes
 - RunStrelka - Run Strelka for Variant Calling
 - RunManta - Run Manta for Structural Variant Calling
 - alleleCount - Run Ascat for CNV
 - convertAlleleCounts - Run Ascat for CNV
 - runASCAT - Run Ascat for CNV