# Usage

I would recommand to run Nextflow within a [screen](https://www.gnu.org/software/screen/) or [tmux](https://tmux.github.io/) session.
It is recommended to run only one instance of Sarek for one patient in the same directory.
Sarek uses several scripts, a wrapper is currently being made to simplify the command lines.
Currently the typical reduced command lines are:

```bash
nextflow run SciLifeLab/Sarek/main.nf --sample <file.tsv> --step <step>
nextflow run SciLifeLab/Sarek/germlineVC.nf --sample <file.tsv> --tools <tool>
nextflow run SciLifeLab/Sarek/somaticVC.nf --sample <file.tsv> --tools <tool>
nextflow run SciLifeLab/Sarek/annotate.nf --tools <tool> (--annotateTools <tools>||--annotateVCF <vcfs>)
nextflow run SciLifeLab/Sarek/runMultiQC.nf
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

Disable all QC tools and MultiQC to generate a HTML report.

### --onlyQC

Run only QC tools and MultiQC to generate a HTML report.

### --outDir

Choose an output directory

### --project `ProjectID`

Specify a project number ID on a UPPMAX cluster. (optionnal if not on such a cluster)

### --sample `file.tsv`

Use the given TSV file as sample (cf [TSV documentation](TSV.md)).

### --step `step`

Choose from wich step the workflow will start. Choose only one step. Possible values are:

- mapping (default, will start workflow with FASTQ files)
- realign (will start workflow with BAM files (with T/N BAMs that were not realigned together))
- recalibrate (will start workflow with BAM files and Recalibration Tables (Only with T/N BAMs that were realigned together))

`--step` option is case insensitive to avoid easy introduction of errors when choosing a step.
### --test

Test run Sarek on a smaller dataset, that way you don't have to specify `--sample data/tsv/tiny.tsv`

### --tools `tool1[,tool2,tool3...]`

Choose which tools will be used in the workflow. Different tools to be separated by commas. Possible values are:

- haplotypecaller (use `HaplotypeCaller` for VC) (germlineVC)
- manta (use `Manta` for SV) (germlineVC,somaticVC)
- strelka (use `Strelka` for VC) (germlineVC,somaticVC)
- ascat (use `ASCAT` for CNV) (somaticVC)
- mutect1 (use `MuTect1` for VC) (somaticVC)
- mutect2 (use `MuTect2` for VC) (somaticVC)
- snpeff (use `snpEff` for Annotation) (annotate)
- vep (use `VEP` for Annotation) (annotate)

`--tools` option is case insensitive to avoid easy introduction of errors when choosing tools. So you can write `--tools mutect2,ascat` or `--tools MuTect2,ASCAT` without worrying about case sensitivity.

### --annotateTools `tool1[,tool2,tool3...]`

Choose which tools to annotate. Different tools to be separated by commas. Possible values are:
- haplotypecaller (Annotate `HaplotypeCaller` output)
- manta (Annotate `Manta` output)
- mutect1 (Annotate `MuTect1` output)
- mutect2 (Annotate `MuTect2` output)
- strelka (Annotate `Strelka` output)

### --annotateVCF `file1[,file2,file3...]`

Choose vcf to annotate. Different vcfs to be separated by commas.

### --verbose

Display more information about files being processed.

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
nextflow run SciLifeLab/Sarek --sample mysample.tsv -profile myprofile
```

A standard profile is defined in [`nextflow.config`](../nextflow.config). You can use the files in the [`configuration/`](../configuration) directory as a base to make a new `.config` file that you can specify directly (or add as a profile):

```bash
nextflow run SciLifeLab/Sarek --sample mysample.tsv -c config/milou.config
```

## Update to latest version

To update workflow to the latest version use:

```bash
nextflow pull SciLifeLab/Sarek
```

## Run the latest version

If there is a feature or bugfix you want to use in a resumed or re-analyzed run, you have to update the workflow to the latest version. By default it is not updated automatically, so use something like:

```bash
nextflow run -latest SciLifeLab/Sarek/main.nf ... -resume
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
