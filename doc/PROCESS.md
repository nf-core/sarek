# Workflow processes
Several processes are run within the workflow. We divide them for the moment into 2 main steps:

## Preprocessing:
- Mapping - Map reads with BWA
- MergeBam - Merge BAMs if multilane samples
- RenameSingleBam - Rename BAM if non-multilane sample
- MarkDuplicates - using Picard
- CreateIntervals - using GATK
- Realign - using GATK
- CreateRecalibrationTable - using GATK
- RecalibrateBam - using GATK

## Variant Calling:
- RunMutect1 - run MuTect1 on multiple intervals
- concatFiles - merge MuTect1 results
- RunMutect2 - run MuTect1 on multiple intervals
- concatFiles - merge MuTect2 results
- VarDict - run VarDict on multiple intervals
- VarDictCollatedVCF - merge Vardict results
- RunStrelka - run Strelka
- Manta - run Manta
- alleleCount - preprocess runASCAT
- convertAlleleCounts - preprocess runASCAT
- runASCAT - run Ascat
- RunHaplotypeCaller - run VarDict on multiple intervals
- concatFiles - merge HaplotypeCaller results