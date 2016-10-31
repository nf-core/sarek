[![Stories in Ready](https://badge.waffle.io/SciLifeLab/CAW.png?label=ready&title=Ready)](https://waffle.io/SciLifeLab/CAW)

# CAW
Nextflow Cancer Analysis Workflow Prototype developed at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

## Version
0.8.3

## Authors
- Sebastian DiLorenzo (@Sebastian-D)
- Jesper Eisfeldt (@J35P312)
- Maxime Garcia (@MaxUlysse)
- Szilveszter Juhos (@szilvajuhos)
- Max Käller (@gulfshores)
- Malin Larsson (@malinlarsson)
- Björn Nystedt (@bjornnystedt)
- Pall Olason (@pallolason)
- Pelin Sahlén (@pelinakan)

## Installation and first execution
See the [Nextflow documentation from SciLifeLab](https://github.com/SciLifeLab/NGI-NextflowDocs/README.md)

See the [workflow installation documentation](doc/UPPMAX.md)

## Usage
I would recommand to run Nextflow within a screen session (cf [help on screen](https://www.howtoforge.com/linux_screen)).
```bash
nextflow run SciLifeLab/CAW --sample <file.tsv> [--steps STEP[,STEP]]
```
All variables and parameters are specified in the config (cf [configuration options](#config)) and the sample files.

### Steps
To configure which processes will be runned or skipped in the workflow. Different steps to be separated by commas.
Possible values are:
- preprocessing (default, will start workflow with FASTQ files)
- recalibrate (will start workflow with non-recalibrated BAM files)
- skipPreprocessing (will start workflow with recalibrated BAM files)
- MuTect1 (use MuTect1 for VC)
- MuTect2 (use MuTect2 for VC)
- VarDict (use VarDict for VC)
- Strelka (use Strelka for VC)
- HaplotypeCaller (use HaplotypeCaller for normal bams VC)
- Manta (use Manta for SV)
- ascat (use ascat for CNV)

## Verbose
To have more information about files being processed, you can use the verbose option
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv --steps preprocessing --verbose
```

## Nextflow parameters
### config
More informations on [Nextflow documentation](https://www.nextflow.io/docs/latest/basic.html#configuration-options)
```bash
-c <file.config>
```
If no config file is specified, Nextflow will look for one in Nextflow intallation `$NXF_HOME/config` of for one in the current directory `nextflow.config`.

The config file provided as an example is a [config file](https://raw.githubusercontent.com/SciLifeLab/CAW/master/config/milou.config) specific to Swedish UPPMAX milou cluster, but can be easily modified to suit any clusters.

You can use this file as an example to make your own config file. And you can even if needed make several config files (for example if you want to have a config file for each UPPMAX project identifier).

### clean
Use `nextflow clean -f` to remove everything contained in the `work` directory. Do not worry, non-recalibrated bam, indexes and recalibration tables as well as recalibrated bams and index are stored respectively in the `Preprocessing/NonRecalibrated` and `Preprocessing/Recalibrated` directories. And variant calling files are stored in the `VariantCalling` directory.
```bash
nextflow clean -f
```

### resume
Use `-resume` to restart the workflow where it last failed.
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv --steps preprocessing -resume
```

### info
Use `info` to get information about the workflow.
```bash
nextflow info SciLifeLab/CAW
```

### pull
Use `pull` to update the workflow.
```bash
nextflow pull SciLifeLab/CAW
```

## Nextflow processes
Several processes are run within Nextflow
We divide them for the moment into 2 main steps

### Preprocessing:
- Mapping - Map reads with BWA
- MergeBam - Merge BAMs if multilane samples
- RenameSingleBam - Rename BAM if non-multilane sample
- MarkDuplicates - using Picard
- CreateIntervals - using GATK
- Realign - using GATK
- CreateRecalibrationTable - using GATK
- RecalibrateBam - using GATK

### Variant Calling:
- RunMutect2 - using MuTect2 shipped in GATK v3.6
- concatFiles - merge MuTect2 results
- VarDict - run VarDict on multiple intervals
- VarDictCollatedVCF - merge Vardict results
- RunStrelka - using Strelka 1.0.15
- Manta - run Manta 1.0.0

## TSV file for sample
It's a Tab Separated Value file, based on: `subject status sample lane fastq1 fastq2`, `subject status sample bam bai recal` or `subject status sample bam bai`
Quite straight-forward:
- `subject` is the ID of the Patient
- `status` is the status of the Patient, (0 for Normal and 1 for Tumor)
- `sample` is the Sample ID (It is possible to have more than one tumor sample for each patient)
- `lane` is used when the sample is multiplexed on several lanes
- `fastq1` is the path to the first pair of the fastq file
- `fastq2` is the path to the second pair of the fastq file
- `bam` is the realigned bam file
- `bai` is the index
- `recal` is the recalibration table

## Example of TSV files
See the [workflow TSV example documentation](doc/EXAMPLE.md)

## Tools and dependencies
- nextflow 0.17.3
- bwa 0.7.8
- samtools 1.3
- picard 1.118
- GATK 3.6
- R 3.2.3
- gcc 4.9.2
- perl 5.18.4
- strelka 1.0.15
- manta 1.0.0

## Use cases
See the [workflow use cases documentation](doc/USE_CASES.md)
