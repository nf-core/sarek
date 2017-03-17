# Usage

I would recommand to run Nextflow within a [screen](https://www.gnu.org/software/screen/) or [tmux](https://tmux.github.io/) session. It is recommanded to run only one instance of CAW for one patient in the same directory. The typical reduced command line is:

```bash
nextflow run SciLifeLab/CAW --sample <file.tsv> --step <step> --tools <tool>
```

All parameters, options and variables can be specified with configuration files and profile (cf [configuration documentation](#profiles)).

## Options

### --callName `Name`

Specify a name for MultiQC report (optionnal)

### --contactMail `email`

Specify an email for MultiQC report (optionnal)

### --help

Display help

### --project `ProjectID`

Specify a project number ID on a UPPMAX cluster. (optionnal if not on such a cluster)

### --sample `file.tsv`

Use the given TSV file as sample (cf [TSV documentation](TSV.md)).

### --step `step`

Choose from wich step the workflow will start. Choose only one step. Possible values are:

- preprocessing (default, will start workflow with FASTQ files)
- realign (will start workflow with BAM files (with T/N BAMs that were not realigned together))
- recalibrate (will start workflow with BAM files and Recalibration Tables (Only with T/N BAMs that were realigned together))
- skipPreprocessing (will skip entire preprocessing (Only with T/N BAMs that were realigned together))

### --test

Test run CAW on a smaller dataset, that way you don't have to specify `--sample data/tsv/tiny.tsv --intervals repeats/tiny.list`

### --tools `tool1[,tool2,tool3...]`

Choose which tools will be used in the workflow. Different tools to be separated by commas. Possible values are:

- Ascat (use ascat for CNV)
- HaplotypeCaller (use HaplotypeCaller for VC)
- Manta (use Manta for SV)
- MultiQC (Make a QC report)
- MuTect1 (use MuTect1 for VC)
- MuTect2 (use MuTect2 for VC)
- Strelka (use Strelka for VC)
- VarDict (use VarDict for VC)
- snpEff (use snpEff for Annotation)

### --verbose

Display more information about files being processed.

### --version

Display version number and information.

## Parameters
Simpler to specify in the config file.

### --runTime `time`
### --singleCPUMem `memory`

## References [(cf [References documentation](REFERENCES.md))]
Could be usefull if you wish to change one reference for testing.

### --acLoci `file`

## COSMIC files
- --cosmic `file`
- --cosmicIndex `file`

### Files from the GATK Bundle
- --dbsnp `file`
- --dbsnpIndex `file`
- --kgIndels `file`
- --kgIndex `file`
- --genome `file`
- --genomeDict `file`
- --genomeIndex `file`
- --millsIndels `file`
- --millsIndex `file`

### BWA indexes
- --genomeAmb `file`
- --genomeAnn `file`
- --genomeBwt `file`
- --genomePac `file`
- --genomeSa `file`

## --intervals `file`

## --snpeffDb `db`
	Which database to use for snpEff

## --vardictHome `path`
	Path to Vardict

# Nextflow options

See the [options documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/OPTIONS.md)

## Profiles

More informations on the [SciLifeLab Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md). The default profile is `standard`. If you want you can use your own profile:

```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv -profile myprofile
```

A standard profile is defined in [`nextflow.config`](../nextflow.config). You can use the [`config/milou.config`](../config/milou.config) file as a base to make a new `config` file that you can specify directly (or add as a profile):

```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv -c config/milou.config
```

## Update to latest version

To update workflow to the latest version use:

```bash
nextflow pull SciLifeLab/CAW
```

## Run the latest version

If there is a feature or bugfix you want to use in a resumed or re-analyzed run, you have to update the workflow to the latest version. By default it is not updated automatically, so use something like:

```bash
nextflow run -latest SciLifeLab/CAW --sample mysample.tsv -resume
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link] [![](images/NGI-final-small.png "NGI")][ngi-link]

[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: http://www.scilifelab.se/
