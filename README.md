[![Stories in Ready](https://badge.waffle.io/SciLifeLab/CAW.png?label=ready&title=Ready)](https://waffle.io/SciLifeLab/CAW)

# CAW
Nextflow Cancer Analysis Workflow Prototype developed at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

## Version
0.8.5

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
To use this pipeline, you need to have a working version of Nextflow installed.
- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [workflow first execution documentation](doc/UPPMAX.md)

## Usage
I would recommand to run Nextflow within a screen session (cf [help on screen](https://www.howtoforge.com/linux_screen)). The typical command line is:
```bash
nextflow run SciLifeLab/CAW --sample <file.tsv>
```
All variables and parameters are specified in the config (cf [configuration documentation](#config)).

### Project
To specify your UPPMAX project number ID. It can also be specified in your `config` file.
```bash
nextflow run SciLifeLab/CAW --sample <file.tsv> --project <UPPMAX_Project>
```

### Steps
To configure which processes will be runned or skipped in the workflow. Different steps to be separated by commas. Possible values are:
- preprocessing (default, will start workflow with FASTQ files)
- realign (will start workflow with BAM files (with T/N BAMs that were not realigned together))
- skipPreprocessing (will skip entire preprocessing (Only with T/N BAMs that were realigned together))
- MuTect1 (use MuTect1 for VC)
- MuTect2 (use MuTect2 for VC)
- VarDict (use VarDict for VC)
- Strelka (use Strelka for VC)
- HaplotypeCaller (use HaplotypeCaller for normal bams VC)
- Manta (use Manta for SV)
- ascat (use ascat for CNV)

### Verbose
To have more information about files being processed, you can use the verbose option
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv --verbose
```

## Nextflow options
See the [options documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/OPTIONS.md)

## config and profiles
More informations on the [SciLifeLab Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md). The default profile is `standard`
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv -profile standard
```
A standard profile is defined in [`nextflow.config`](https://raw.githubusercontent.com/SciLifeLab/CAW/master/nextflow.config). You can use the [`config/milou.config`](https://raw.githubusercontent.com/SciLifeLab/CAW/master/config/milou.config) file as a base to make a new `config` file that you can specify directly (or add as a profile):
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv -c config/milou.config
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
It's a Tab Separated Value file, based on: `subject status sample lane fastq1 fastq2` or `subject status sample bam bai`
Quite straight-forward:
- `subject` is the ID of the Patient
- `status` is the status of the Patient, (0 for Normal and 1 for Tumor)
- `sample` is the Sample ID (It is possible to have more than one tumor sample for each patient)
- `lane` is used when the sample is multiplexed on several lanes
- `fastq1` is the path to the first pair of the fastq file
- `fastq2` is the path to the second pair of the fastq file
- `bam` is the bam file
- `bai` is the index

## Example of TSV files
See the [workflow TSV example documentation](doc/EXAMPLE.md)

## Tools and dependencies
- nextflow 0.22.1
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
