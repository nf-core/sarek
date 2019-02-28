# Parameters

A list of all possible parameter that can be used for the different scripts included in Sarek.

## Common for all scripts

### --help

Display help

### --noReports

Disable all QC tools and MultiQC.

### --outDir

Choose an output directory

### --project `ProjectID`

Specify a project number ID on a UPPMAX cluster.
(optional if not on such a cluster)

## --publishDirMode

Specify wich mode `publishDir` directive need to follow, in Sarek link is the default mode.

From [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#publishdir):

| Mode         | Description |
|--------------|-------------|
| copy         | Copies the output files into the published directory |
| copyNoFollow | Copies the output files into the published directory without following symlinks ie. copies the links themselves |
| link         | Creates a hard link in the published directory for each process output file (default) |
| move         | Moves the output files into the published directory. Note: this is only supposed to be used for a terminating process i.e. a process whose output is not consumed by any other downstream process |
| rellink      | Creates a relative symbolic link in the published directory for each process output file |
| symlink      | Creates an absolute symbolic link in the published directory for each process output file |

### --sample `file.tsv`

Use the given TSV file as sample (cf [TSV documentation](TSV.md)).
Is not used for `annotate.nf` and `runMultiQC.nf`.

### --tools `tool1[,tool2,tool3...]`

Choose which tools will be used in the workflow.
Different tools to be separated by commas.
Possible values are:

- haplotypecaller (use `HaplotypeCaller` for VC) (germlineVC.nf)
- manta (use `Manta` for SV) (germlineVC.nf,somaticVC.nf)
- strelka (use `Strelka` for VC) (germlineVC.nf,somaticVC.nf)
- ascat (use `ASCAT` for CNV) (somaticVC.nf)
- mutect2 (use `MuTect2` for VC) (somaticVC.nf)
- snpeff (use `snpEff` for Annotation) (annotate.nf)
- vep (use `VEP` for Annotation) (annotate.nf)

`--tools` option is case insensitive to avoid easy introduction of errors when choosing tools.
So you can write `--tools mutect2,ascat` or `--tools MuTect2,ASCAT` without worrying about case sensitivity.

### --awsqueue `BatchQueueName`

Only required if you use the awsbatch profile. This parameter specifies the queue for which jobs are submitted in AWS Batch.

### --awsqueue_tiny `BatchQueueName`

Only used if you use the awsbatch profile. This parameter specifies a queue used for certain small jobs that might still require a significant amount of disk storage.

### --localReportDir `Directory`

Only used if you use the awsbatch profile. This parameter specifies an output directory for nextflow reports, such as Sarek_timeline.html, which currently is not fully supported to store on s3.

### --verbose

Display more information about files being processed.

## Preprocessing script (`main.nf`)
### --step `step`

Choose from wich step the workflow will start.
Choose only one step.
Possible values are:

- mapping (default, will start workflow with FASTQ files)
- recalibrate (will start workflow with BAM files and Recalibration Tables

`--step` option is case insensitive to avoid easy introduction of errors when choosing a step.

### --test

Test run Sarek on a smaller dataset, that way you don't have to specify `--sample Sarek-data/testdata/tsv/tiny.tsv`

### --onlyQC

Run only QC tools and MultiQC to generate a HTML report.

## Annotate script (`annotate.nf`)

### --annotateTools `tool1[,tool2,tool3...]`

Choose which tools to annotate.
Different tools to be separated by commas.
Possible values are:
- haplotypecaller (Annotate `HaplotypeCaller` output)
- manta (Annotate `Manta` output)
- mutect2 (Annotate `MuTect2` output)
- strelka (Annotate `Strelka` output)

### --annotateVCF `file1[,file2,file3...]`

Choose vcf to annotate.
Different vcfs to be separated by commas.

### --annotate_cache

Use available cache on local machine.
Works only when used with `--snpEff_cache` and/or `--vep_cache` to specify the cache directories.
See [annotation documentation](ANNOTATION.md) for more information.

### --snpEff_cache `/Path/To/snpEffCache`

Specify the location of snpEff cache.
Works only when used with `--annotate_cache`.
See [annotation documentation](ANNOTATION.md) for more information.

### --vep_cache `/Path/To/VEPcache`

Specify the location of VEP cache.
Works only when used with `--annotate_cache`.
See [annotation documentation](ANNOTATION.md) for more information.

## MultiQC script (`runMultiQC.nf`)
### --callName `Name`

Specify a name for MultiQC report (optional)

### --contactMail `email`

Specify an email for MultiQC report (optional)


## References

For most use cases, the reference information is already in the configuration file [`conf/genomes.config`](https://github.com/SciLifeLab/Sarek/blob/master/conf/genomes.config).
However, if needed, you can specify any reference file at the command line.

### --acLoci `acLoci file`

### --bwaIndex `bwaIndex file`

### --dbsnp `dbsnp file`

### --dbsnpIndex `dbsnpIndex file`

### --genomeDict `genomeDict file`

### --genomeFile `genomeFile file`

### --genomeIndex `genomeIndex file`

### --intervals `intervals file`

### --knownIndels `knownIndels file`

### --knownIndelsIndex `knownIndelsIndex file`

### --snpeffDb `snpeffDb file`

## Hardware Parameters

For most use cases, the reference information is already in the appropriate [configuration files](https://github.com/SciLifeLab/Sarek/blob/master/conf/).
However, it is still possible to specify these parameters at the command line as well.

### --runTime `time`

### --singleCPUMem `memory`

### --totalMemory `memory`
