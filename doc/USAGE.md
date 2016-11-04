# Usage
I would recommand to run Nextflow within a screen session (cf [help on screen](https://www.howtoforge.com/linux_screen)). The typical command line is:
```bash
nextflow run SciLifeLab/CAW --sample <file.tsv>
```
All variables and parameters are specified in the config (cf [configuration documentation](#profiles)).
All samples are specified in the TSV files (cf [TSV documentation](TSV.md)).

## Steps
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

## Project
To specify your UPPMAX project number ID. It can also be specified in your `config` file (cf [configuration documentation](#Profiles)).
```bash
nextflow run SciLifeLab/CAW --sample <file.tsv> --project <UPPMAX_Project>
```

## Verbose
To have more information about files being processed, you can use the verbose option:
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv --verbose
```

# Nextflow options
See the [options documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/OPTIONS.md)

## Profiles
More informations on the [SciLifeLab Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md). The default profile is `standard`. If you want you can use your own profile:
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv -profile myprofile
```
A standard profile is defined in [`nextflow.config`](https://raw.githubusercontent.com/SciLifeLab/CAW/master/nextflow.config). You can use the [`config/milou.config`](https://raw.githubusercontent.com/SciLifeLab/CAW/master/config/milou.config) file as a base to make a new `config` file that you can specify directly (or add as a profile):
```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv -c config/milou.config
```
