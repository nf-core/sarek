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

### --noReports

Disable QC tools and MultiQC to generate a HTML report.

### --project `ProjectID`

Specify a project number ID on a UPPMAX cluster. (optionnal if not on such a cluster)

### --sample `file.tsv`

Use the given TSV file as sample (cf [TSV documentation](TSV.md)).

### --step `step`

Choose from wich step the workflow will start. Choose only one step. Possible values are:

- mapping (default, will start workflow with FASTQ files)
- realign (will start workflow with BAM files (with T/N BAMs that were not realigned together))
- recalibrate (will start workflow with BAM files and Recalibration Tables (Only with T/N BAMs that were realigned together))
- variantcalling (will skip entire preprocessing (Only with T/N BAMs that were realigned together))
- annotate (will annotate Variant Calling output. By default it will try to annotate all available vcfs. Use with ```--annotateTools``` or ```--annotateVCF``` to specify what to annotate)

`--step` option is case insensitive to avoid easy introduction of errors when choosing a step. So you can write `--step variantCalling` or `--step variantcalling` without worrying about case sensitivity.

### --test

Test run CAW on a smaller dataset, that way you don't have to specify `--sample data/tsv/tiny.tsv`

### --tools `tool1[,tool2,tool3...]`

Choose which tools will be used in the workflow. Different tools to be separated by commas. Possible values are:

- ascat (use ascat for CNV)
- haplotypecaller (use HaplotypeCaller for VC)
- manta (use Manta for SV)
- mutect1 (use MuTect1 for VC)
- mutect2 (use MuTect2 for VC)
- strelka (use Strelka for VC)
- snpeff (use snpEff for Annotation)
- vep (use VEP for Annotation)

`--tools` option is case insensitive to avoid easy introduction of errors when choosing tools. So you can write `--tools mutect2,snpEff` or `--tools MuTect2,snpeff` without worrying about case sensitivity.

### --annotateTools `tool1[,tool2,tool3...]`

Choose which tools to annotate. Different tools to be separated by commas. Possible values are:
- haplotypecaller (Annotate HaplotypeCaller output)
- manta (Annotate Manta output)
- mutect1 (Annotate MuTect1 output)
- mutect2 (Annotate MuTect2 output)
- strelka (Annotate Strelka output)

### --annotateVCF `file1[,file2,file3...]`

Choose which vcf to annotate. Different vcf to be separated by commas.

### --verbose

Display more information about files being processed.

### --version

Display version number and information.

## Containers

### --containerPath `Path to the singularity containers (default=containers/)`

### --repository `Docker-hub repository (default=maxulysse)`

### --tag `tag of the containers to use (default=current version)`

## References

If needed, you can specify each reference file by command line.

### --acLoci `acLoci file`

### --bwaIndex `bwaIndex file`

### --cosmic `cosmic file`

### --cosmicIndex `cosmicIndex file`

### --dbsnp `dbsnp file`

### --dbsnpIndex `dbsnpIndex file`

### --genomeDict `genomeDict file`

### --genomeFile `genomeFile file`

### --genomeIndex `genomeIndex file`

### --intervals `intervals file`

### --knownIndels `knownIndels file`

### --knownIndelsIndex `knownIndelsIndex file`

### --snpeffDb `snpeffDb file`

## Parameters

Simpler to specify in the configuration files, but it's still possible to specify every thing in the command line.

### --runTime `time`

### --singleCPUMem `memory`

### --totalMemory `memory`

# Nextflow options

See the [options documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/OPTIONS.md)

## Profiles

More informations on the [SciLifeLab Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md). The default profile is `standard`. You can use your own profile:

```bash
nextflow run SciLifeLab/CAW --sample mysample.tsv -profile myprofile
```

A standard profile is defined in [`nextflow.config`](../nextflow.config). You can use the files in the [`configuration/`](../configuration) directory as a base to make a new `.config` file that you can specify directly (or add as a profile):

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

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
