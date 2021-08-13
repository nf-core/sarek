# nf-core/sarek: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/sarek/usage](https://nf-co.re/sarek/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

Sarek is a workflow designed to detect variants on whole genome or targeted sequencing data.
Initially designed for Human, and Mouse, it can work on any species with a reference genome.
Sarek can also handle tumour / normal pairs and could include additional relapses.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```console
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column         | Description                                                                                                                                                                            |
|----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sample`       | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`      | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`      | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/sarek --input samplesheet.csv --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/sarek
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `2.6.1`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.6.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile.
Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of `Docker` or `Singularity` containers for full pipeline reproducibility, however when this is not possible, `Conda` is also supported.

The pipeline also dynamically loads configurations from [github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time.
For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`.
This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `conda`
    * A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/software/star/align/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB. The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: STAR_ALIGN {
        memory = 100.GB
    }
}
```

> **NB:** We specify just the process name i.e. `STAR_ALIGN` in the config file and not the full task name string that is printed to screen in the error message or on the terminal whilst the pipeline is running i.e. `RNASEQ:ALIGN_STAR:STAR_ALIGN`. You may get a warning suggesting that the process selector isn't recognised but you can ignore that if the process name has been specified correctly. This is something that needs to be fixed upstream in core Nextflow.

### Tool-specific options

For the ultimate flexibility, we have implemented and are using Nextflow DSL2 modules in a way where it is possible for both developers and users to change tool-specific command-line arguments (e.g. providing an additional command-line argument to the `STAR_ALIGN` process) as well as publishing options (e.g. saving files produced by the `STAR_ALIGN` process that aren't saved by default by the pipeline). In the majority of instances, as a user you won't have to change the default options set by the pipeline developer(s), however, there may be edge cases where creating a simple custom config file can improve the behaviour of the pipeline if for example it is failing due to a weird error that requires setting a tool-specific parameter to deal with smaller / larger genomes.

The command-line arguments passed to STAR in the `STAR_ALIGN` module are a combination of:

* Mandatory arguments or those that need to be evaluated within the scope of the module, as supplied in the [`script`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L49-L55) section of the module file.

* An [`options.args`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L56) string of non-mandatory parameters that is set to be empty by default in the module but can be overwritten when including the module in the sub-workflow / workflow context via the `addParams` Nextflow option.

The nf-core/rnaseq pipeline has a sub-workflow (see [terminology](https://github.com/nf-core/modules#terminology)) specifically to align reads with STAR and to sort, index and generate some basic stats on the resulting BAM files using SAMtools. At the top of this file we import the `STAR_ALIGN` module via the Nextflow [`include`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/subworkflows/nf-core/align_star.nf#L10) keyword and by default the options passed to the module via the `addParams` option are set as an empty Groovy map [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/subworkflows/nf-core/align_star.nf#L5); this in turn means `options.args` will be set to empty by default in the module file too. This is an intentional design choice and allows us to implement well-written sub-workflows composed of a chain of tools that by default run with the bare minimum parameter set for any given tool in order to make it much easier to share across pipelines and to provide the flexibility for users and developers to customise any non-mandatory arguments.

When including the sub-workflow above in the main pipeline workflow we use the same `include` statement, however, we now have the ability to overwrite options for each of the tools in the sub-workflow including the [`align_options`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/workflows/rnaseq.nf#L225) variable that will be used specifically to overwrite the optional arguments passed to the `STAR_ALIGN` module. In this case, the options to be provided to `STAR_ALIGN` have been assigned sensible defaults by the developer(s) in the pipeline's [`modules.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L70-L74) and can be accessed and customised in the [workflow context](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/workflows/rnaseq.nf#L201-L204) too before eventually passing them to the sub-workflow as a Groovy map called `star_align_options`. These options will then be propagated from `workflow -> sub-workflow -> module`.

As mentioned at the beginning of this section it may also be necessary for users to overwrite the options passed to modules to be able to customise specific aspects of the way in which a particular tool is executed by the pipeline. Given that all of the default module options are stored in the pipeline's `modules.config` as a [`params` variable](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L24-L25) it is also possible to overwrite any of these options via a custom config file.

Say for example we want to append an additional, non-mandatory parameter (i.e. `--outFilterMismatchNmax 16`) to the arguments passed to the `STAR_ALIGN` module. Firstly, we need to copy across the default `args` specified in the [`modules.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L71) and create a custom config file that is a composite of the default `args` as well as the additional options you would like to provide. This is very important because Nextflow will overwrite the default value of `args` that you provide via the custom config.

As you will see in the example below, we have:

* appended `--outFilterMismatchNmax 16` to the default `args` used by the module.
* changed the default `publish_dir` value to where the files will eventually be published in the main results directory.
* appended `'bam':''` to the default value of `publish_files` so that the BAM files generated by the process will also be saved in the top-level results directory for the module. Note: `'out':'log'` means any file/directory ending in `out` will now be saved in a separate directory called `my_star_directory/log/`.

```nextflow
params {
    modules {
        'star_align' {
            args          = "--quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --readFilesCommand zcat --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outFilterMismatchNmax 16"
            publish_dir   = "my_star_directory"
            publish_files = ['out':'log', 'tab':'log', 'bam':'']
        }
    }
}
```

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

    * For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Singularity:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```

## Troubleshooting

### TSV file

> **NB** Delimiter is the tab (`\t`) character, and no header is required

There are different kinds of `TSV` files that can be used as input, depending on the input files available (`FASTQ`, `unmapped BAM`, `recalibrated BAM`...).
The `TSV` file should correspond to the correct step.
For all possible `TSV` files, described in the next sections, here is an explanation of what the columns refer to:

`Sarek` auto-generates `TSV` files for all and for each individual samples, depending of the options specified.

* `subject` designates the subject, it should be the ID of the subject, and it must be unique for each subject, but one subject can have multiple samples (e.g.
normal and tumor)
* `sex` are the sex chromosomes of the subject, (ie `XX`, `XY`...) and will only be used for Copy-Number Variation in a tumor/pair.
* `status` is the status of the measured sample, (`0` for Normal or `1` for Tumor)
* `sample` designates the sample, it should be the ID of the sample (it is possible to have more than one tumor sample for each subject, i.e. a tumor and a relapse), it must be unique, but samples can have multiple lanes (which will later be merged)
* `lane` is used when the sample is multiplexed on several lanes, it must be unique for each lane in the same sample (but does not need to be the original lane name), and must contain at least one character
* `fastq1` is the path to the first pair of the `FASTQ` file
* `fastq2` is the path to the second pair of the `FASTQ` file
* `bam` is the path to the `BAM` file
* `bai` is the path to the `BAM` index file
* `recaltable` is the path to the recalibration table
* `mpileup` is the path to the mpileup file

It is recommended to use the absolute path of the files, but relative path should also work.

If necessary, a tumor sample can be associated to a normal sample as a pair, if specified with the same `subject`and a different `sample`.
An additional tumor sample (such as a relapse for example), can be added if specified with the same `subject` and a different `sample`.

`Sarek` will output results in a different directory for each sample.
If multiple samples are specified in the `TSV` file, `Sarek` will consider all files to be from different samples.
Multiple `TSV` files can be specified if the path is enclosed in quotes.

Output from Variant Calling and/or Annotation will be in a specific directory for each sample (or normal/tumor pair if applicable).

#### --input &lt;FASTQ&gt; --step mapping

The `TSV` file to start with the mapping step (`--step mapping`) with paired-end `FASTQs` should contain the columns:

`subject sex status sample lane fastq1 fastq2`

In this example (`example_fastq.tsv`), there are 3 read groups.

| | | | | | | |
|-|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID|1|/samples/normal1_1.fastq.gz|/samples/normal1_2.fastq.gz|
|SUBJECT_ID|XX|0|SAMPLE_ID|2|/samples/normal2_1.fastq.gz|/samples/normal2_2.fastq.gz|
|SUBJECT_ID|XX|0|SAMPLE_ID|3|/samples/normal3_1.fastq.gz|/samples/normal3_2.fastq.gz|

```bash
--input example_fastq.tsv
```

Or, for a normal/tumor pair:

In this example (`example_pair_fastq.tsv`), there are 3 read groups for the normal sample and 2 for the tumor sample.

| | | | | | | |
|-|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|1|/samples/normal1_1.fastq.gz|/samples/normal1_2.fastq.gz|
|SUBJECT_ID|XX|0|SAMPLE_ID1|2|/samples/normal2_1.fastq.gz|/samples/normal2_2.fastq.gz|
|SUBJECT_ID|XX|0|SAMPLE_ID1|3|/samples/normal3_1.fastq.gz|/samples/normal3_2.fastq.gz|
|SUBJECT_ID|XX|1|SAMPLE_ID2|1|/samples/tumor1_1.fastq.gz|/samples/tumor1_2.fastq.gz|
|SUBJECT_ID|XX|1|SAMPLE_ID2|2|/samples/tumor2_1.fastq.gz|/samples/tumor2_2.fastq.gz|

```bash
--input example_pair_fastq.tsv
```

#### --input &lt;uBAM&gt; --step mapping

The `TSV` file to start with the mapping step (`--step mapping`) with `unmapped BAM` files should contain the columns:

`subject sex status sample lane bam`

In this example (`example_ubam.tsv`), there are 3 read groups.

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID|1|/samples/normal_1.bam|
|SUBJECT_ID|XX|0|SAMPLE_ID|2|/samples/normal_2.bam|
|SUBJECT_ID|XX|0|SAMPLE_ID|3|/samples/normal_3.bam|

```bash
--input example_ubam.tsv
```

Or, for a normal/tumor pair:

In this example (`example_pair_ubam.tsv`), there are 3 read groups for the normal sample and 2 for the tumor sample.

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|1|/samples/normal_1.bam|
|SUBJECT_ID|XX|0|SAMPLE_ID1|2|/samples/normal_2.bam|
|SUBJECT_ID|XX|0|SAMPLE_ID1|3|/samples/normal_3.bam|
|SUBJECT_ID|XX|1|SAMPLE_ID2|1|/samples/tumor_1.bam|
|SUBJECT_ID|XX|1|SAMPLE_ID2|2|/samples/tumor_2.bam|

```bash
--input example_pair_ubam.tsv
```

#### --input &lt;TSV&gt; --step prepare_recalibration

To start from the preparation of the recalibration step (`--step prepare_recalibration`), a `TSV` file needs to be given as input containing the paths to the `non-recalibrated BAM` files.
The `Sarek`-generated `TSV` file is stored under `results/Preprocessing/TSV/duplicates_marked_no_table.tsv` and will automatically be used as an input when specifying the parameter `--step prepare_recalibration`.

The `TSV` contains the following columns:

`subject sex status sample bam bai`

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID|/samples/normal.md.bam|/samples/normal.md.bai|

Or, for a normal/tumor pair:

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|/samples/normal.md.bam|/samples/normal.md.bai|
|SUBJECT_ID|XX|1|SAMPLE_ID2|/samples/tumor.md.bam|/samples/tumor.md.bai|

#### --input &lt;TSV&gt; --step prepare_recalibration --skip_markduplicates

The `Sarek`-generated `TSV` file is stored under `results/Preprocessing/TSV/mapped.tsv` and will automatically be used as an input when specifying the parameter `--step prepare_recalibration --skip_markduplicates`.
The `TSV` file contains the same columns, but the content is slightly different:

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID|/samples/normal.bam|/samples/normal.bai|

Or, for a normal/tumor pair:

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|/samples/normal.bam|/samples/normal.bai|
|SUBJECT_ID|XX|1|SAMPLE_ID2|/samples/tumor.bam|/samples/tumor.bai|

#### --input &lt;TSV&gt; --step recalibrate

To start from the recalibrate step (`--step recalibrate`), a `TSV` file needs to be given as input containing the paths to the `non-recalibrated BAM` file and the associated recalibration table.
The `Sarek`-generated `TSV` file is stored under `results/Preprocessing/TSV/duplicates_marked.tsv` and will automatically be used as an input when specifying the parameter `--step recalibrate`.

The `TSV` contains the following columns:

`subject sex status sample bam bai recaltable`

| | | | | | | |
|-|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID|/samples/normal.md.bam|/samples/normal.md.bai|/samples/normal.recal.table|

Or, for a normal/tumor pair:

| | | | | | | |
|-|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|/samples/normal.md.bam|/samples/normal.md.bai|/samples/normal.recal.table|
|SUBJECT_ID|XX|1|SAMPLE_ID2|/samples/tumor.md.bam|/samples/tumor.md.bai|/samples/tumor.recal.table|

#### --input &lt;TSV&gt; --step recalibrate --skip_markduplicates

The `Sarek`-generated `TSV` file is stored under `results/Preprocessing/TSV/mapped_no_duplicates_marked.tsv` and will automatically be used as an input when specifying the parameter `--step recalibrate --skip_markduplicates`.
The `TSV` file contains the same columns, but the content is slightly different:

| | | | | | | |
|-|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID|/samples/normal.bam|/samples/normal.bai|/samples/normal.recal.table|

Or, for a normal/tumor pair:

| | | | | | | |
|-|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|/samples/normal.bam|/samples/normal.bai|/samples/normal.recal.table|
|SUBJECT_ID|XX|1|SAMPLE_ID2|/samples/tumor.bam|/samples/tumor.bai|/samples/tumor.recal.table|

#### --input &lt;TSV&gt; --step variant_calling

To start from the variant calling step (`--step variant_calling`), a `TSV` file needs to be given as input containing the paths to the `recalibrated BAM` file and the associated index.
The `Sarek`-generated `TSV` file is stored under `results/Preprocessing/TSV/recalibrated.tsv` and will automatically be used as an input when specifying the parameter `--step variant_calling`.

The `TSV` file should contain the columns:

`subject sex status sample bam bai`

Here is an example for two samples from the same subject:

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID|/samples/normal.recal.bam|/samples/normal.recal.bai|

Or, for a normal/tumor pair:

| | | | | | |
|-|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|/samples/normal.recal.bam|/samples/normal.recal.bai|
|SUBJECT_ID|XX|1|SAMPLE_ID2|/samples/tumor.recal.bam|/samples/tumor.recal.bai|

#### --input &lt;TSV&gt; --step Control-FREEC

To start from the Control-FREEC step (`--step Control-FREEC`), a `TSV` file needs to be given as input containing the paths to the mpileup files.
The `Sarek`-generated `TSV` file is stored under `results/VariantCalling/TSV/control-freec_mpileup.tsv` and will automatically be used as an input when specifying the parameter `--step Control-FREEC`.

The `TSV` file should contain the columns:

`subject sex status sample mpileup`

Here is an example for one normal/tumor pair from one subjects:

| | | | | |
|-|-|-|-|-|
|SUBJECT_ID|XX|0|SAMPLE_ID1|/samples/normal.pileup|
|SUBJECT_ID|XX|1|SAMPLE_ID2|/samples/tumor.pileup|

### --input &lt;sample/&gt; --step mapping

Use this to specify the location to a directory with `FASTQ` files for the `mapping` step of a single germline sample only.
For example:

```bash
--input </path/to/directory>
```

> **NB** All of the found `FASTQ` files are considered to belong to the same sample.

The input folder, containing the `FASTQ` files for one subject (ID) should be organized into one sub-folder for every sample.
The given directory is searched recursively for `FASTQ` files that are named `*_R1_*.fastq.gz`, and a matching pair with the same name except `_R2_` instead of `_R1_` is expected to exist alongside.
All `FASTQ` files for that sample should be collected here.

```text
ID
+--sample1
+------sample1_<lib>_<flowcell-index>_lane1_R1_1000.fastq.gz
+------sample1_<lib>_<flowcell-index>_lane1_R2_1000.fastq.gz
+------sample1_<lib>_<flowcell-index>_lane2_R1_1000.fastq.gz
+------sample1_<lib>_<flowcell-index>_lane2_R2_1000.fastq.gz
+--sample2
+------sample2_<lib>_<flowcell-index>_lane1_R1_1000.fastq.gz
+------sample2_<lib>_<flowcell-index>_lane1_R2_1000.fastq.gz
+--sample3
+------sample3_<lib>_<flowcell-index>_lane1_R1_1000.fastq.gz
+------sample3_<lib>_<flowcell-index>_lane1_R2_1000.fastq.gz
+------sample3_<lib>_<flowcell-index>_lane2_R1_1000.fastq.gz
+------sample3_<lib>_<flowcell-index>_lane2_R2_1000.fastq.gz
```

`FASTQ` filename structure:

* `<sample>_<lib>_<flowcell-index>_<lane>_R1_<XXX>.fastq.gz` and
* `<sample>_<lib>_<flowcell-index>_<lane>_R2_<XXX>.fastq.gz`

Where:

* `sample` = sample id
* `lib` = identifier of library preparation
* `flowcell-index` = identifier of flow cell for the sequencing run
* `lane` = identifier of the lane of the sequencing run

Read group information will be parsed from `FASTQ` file names according to this:

* `RGID` = "sample_lib_flowcell_index_lane"
* `RGPL` = "Illumina"
* `PU` = sample
* `RGLB` = lib

Each `FASTQ` file pair gets its own read group (`@RG`) in the resulting `BAM` file in the following way.

* The sample name (`SM`) is derived from the the last component of the path given to `--input`.
That is, you should make sure that that directory has a meaningful name! For example, with `--input=/my/fastqs/sample123`, the sample name will be `sample123`.
* The read group id is set to *flowcell.samplename.lane*.
The flowcell id and lane number are auto-detected from the name of the first read in the `FASTQ` file.

### --input &lt;VCF&gt; --step annotate

Input files for Sarek can be specified using the path to a `VCF` file given to the `--input` command only with the annotation step (`--step annotate`).
As `Sarek` will use `bgzip` and `tabix` to compress and index `VCF` files annotated, it expects `VCF` files to be sorted.
Multiple `VCF` files can be specified, using a [glob path](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob), if enclosed in quotes.
For example:

```bash
--step annotate --input "results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,Strelka,TIDDIT}/*.vcf.gz"
```

### Sentieon

Sentieon is a commercial solution to process genomics data with high computing efficiency, fast turnaround time, exceptional accuracy, and 100% consistency.

Please refer to the [nf-core/configs](https://github.com/nf-core/configs#adding-a-new-pipeline-specific-config) repository on how to make a pipeline-specific configuration file based on the [munin-sarek specific configuration file](https://github.com/nf-core/configs/blob/master/conf/pipeline/sarek/munin.config).

Or ask us on the [nf-core Slack](http://nf-co.re/join/slack) on the following channels: [#sarek](https://nfcore.slack.com/channels/sarek) or [#configs](https://nfcore.slack.com/channels/configs).

#### Alignment

> Sentieon BWA matches BWA-MEM with > 2X speedup.

This tool is enabled by default within `Sarek` if both `--sentieon` and `--step mapping` are specified.

#### Germline SNV/INDEL Variant Calling - DNAseq

> Precision FDA award-winning software.
> Matches GATK 3.3-4.1, and without down-sampling.
> Results up to 10x faster and 100% consistent every time.

This tool is enabled within `Sarek` if both `--sentieon` and `--tools DNAseq` are specified.

#### Germline SNV/INDEL Variant Calling - DNAscope

> Improved accuracy and genome characterization.
> Machine learning enhanced filtering producing top variant calling accuracy.

This tool is enabled within `Sarek` if both `--sentieon` and `--tools DNAscope` are specified.

#### Somatic SNV/INDEL Variant Calling - TNscope

> Winner of ICGC-TCGA DREAM challenge.
> Improved accuracy, machine learning enhanced filtering.
> Supports molecular barcodes and unique molecular identifiers.

This tool is enabled within `Sarek` if both `--sentieon` and `--tools TNscope` are specified.

#### Structural Variant Calling

> Germline and somatic SV calling, including translocations, inversions, duplications and large INDELs

This tool is enabled within `Sarek` if both `--sentieon` and `--tools DNAscope` are specified.

### Containers

With `Nextflow DSL2`, each process use its own `Conda` environment or container from `biocontainers`.

For annotation, cache has to be downloaded, or specifically designed containers are available with cache.

`sareksnpeff`, our `snpeff` container is designed using [Conda](https://conda.io/).

[![sareksnpeff-docker status](https://img.shields.io/docker/automated/nfcore/sareksnpeff.svg)](https://hub.docker.com/r/nfcore/sareksnpeff)

Based on [nfcore/base:1.12.1](https://hub.docker.com/r/nfcore/base/tags), it contains:

* **[snpEff](http://snpeff.sourceforge.net/)** 4.3.1t
* Cache for `GRCh37`, `GRCh38`, `GRCm38`, `CanFam3.1` or `WBcel235`

`sarekvep`, our `vep` container is designed using [Conda](https://conda.io/).

[![sarekvep-docker status](https://img.shields.io/docker/automated/nfcore/sarekvep.svg)](https://hub.docker.com/r/nfcore/sarekvep)

Based on [nfcore/base:1.12.1](https://hub.docker.com/r/nfcore/base/tags), it contains:

* **[GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)** 1.0
* **[VEP](https://github.com/Ensembl/ensembl-vep)** 99.2
* Cache for `GRCh37`, `GRCh38`, `GRCm38`, `CanFam3.1` or `WBcel235`

### Using downloaded cache

Both `snpEff` and `VEP` enable usage of cache.
If cache is available on the machine where `Sarek` is run, it is possible to run annotation using cache.
You need to specify the cache directory using `--snpeff_cache` and `--vep_cache` in the command lines or within configuration files.
The cache will only be used when `--annotation_cache` and cache directories are specified (either in command lines or in a configuration file).

Example:

```bash
nextflow run nf-core/sarek --tools snpEff --step annotate --sample <file.vcf.gz> --snpeff_cache </path/to/snpEff/cache> --annotation_cache
nextflow run nf-core/sarek --tools VEP --step annotate --sample <file.vcf.gz> --vep_cache </path/to/VEP/cache> --annotation_cache
```

### Spark related issues

If you have problems running processes that make use of Spark such as ```MarkDuplicates```.
You are probably experiencing issues with the limit of open files in your system.
You can check your current limit by typing the following:

```bash
ulimit -n
```

The default limit size is usually 1024 which is quite low to run Spark jobs.
In order to increase the size limit permanently you can:

Edit the file ```/etc/security/limits.conf``` and add the lines:

```bash
*     soft   nofile  65535
*     hard   nofile  65535
```

Edit the file ```/etc/sysctl.conf``` and add the line:

```bash
fs.file-max = 65535
```

Edit the file ```/etc/sysconfig/docker``` and add the new limits to OPTIONS like this:

```bash
OPTIONS=”—default-ulimit nofile=65535:65535"
```

Re-start your session.

Note that the way to increase the open file limit in your system may be slightly different or require additional steps.

### Download cache

A `Nextflow` helper script has been designed to help downloading `snpEff` and `VEP` caches.
Such files are meant to be shared between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --snpeff_cache </path/to/snpEff/cache> --snpeff_db <snpEff DB version> --genome <GENOME>
nextflow run download_cache.nf --vep_cache </path/to/VEP/cache> --species <species> --vep_cache_version <VEP cache version> --genome <GENOME>
```

### Using VEP CADD plugin

To enable the use of the `VEP` `CADD` plugin:

* Download the `CADD` files
* Specify them (either on the command line, like in the example or in a configuration file)
* use the `--cadd_cache` flag

Example:

```bash
nextflow run nf-core/sarek --step annotate --tools VEP --sample <file.vcf.gz> --cadd_cache \
    --cadd_indels </path/to/CADD/cache/InDels.tsv.gz> \
    --cadd_indels_tbi </path/to/CADD/cache/InDels.tsv.gz.tbi> \
    --cadd_wg_snvs </path/to/CADD/cache/whole_genome_SNVs.tsv.gz> \
    --cadd_wg_snvs_tbi </path/to/CADD/cache/whole_genome_SNVs.tsv.gz.tbi>
```

### Downloading CADD files

An helper script has been designed to help downloading `CADD` files.
Such files are meant to be share between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --cadd_cache </path/to/CADD/cache> --cadd_version <CADD version> --genome <GENOME>
```
