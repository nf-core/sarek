# nf-core/sarek: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/sarek/usage](https://nf-co.re/sarek/usage)

## Introduction

Sarek is a workflow designed to detect germline and somatic variants on whole genome, whole exome, or targeted sequencing data.

Initially designed for human and mouse, it can work on any species if a reference genome is available.
Sarek is designed to handle single samples, such as single-normal or single-tumor samples, and tumor-normal pairs including additional relapses.

## Running the pipeline

### Quickstart

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/sarek --input samplesheet.csv  --outdir <OUTDIR> --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Input: Samplesheet configurations

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` to specify its location. It has to be a comma-separated file with at least 3 columns, and a header row as shown in the examples below.

It is recommended to use the absolute path of the files, but a relative path should also work.

If necessary, a tumor sample can be associated to a normal sample as a pair, if specified with the same `patient` ID, a different `sample`, and the respective `status`.
An additional tumor sample (such as a relapse for example), can be added if specified with the same `patient` ID, a different `sample`, and the `status` value `1`.

`Sarek` will output results in a different directory for _each sample_.
If multiple samples IDs are specified in the `CSV` file, `Sarek` will consider all files to be from different samples.

Output from Variant Calling and/or Annotation will be in a specific directory for each sample and tool configuration (or normal/tumor pair if applicable).

Multiple `CSV` files can be specified if the path is enclosed in quotes.

```console
--input '[path to samplesheet file(s)]'
```

| Column    | Description                                                                                                                                                                                                                                                                                                     |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `patient` | **Custom patient ID**; designates the patient/subject; must be unique for each patient, but one patient can have multiple samples (e.g. normal and tumor).                                                                                                                                                      |
| `sex`  | **Sex chromosomes of the patient**; i.e. XX, XY..., only used for Copy-Number Variation analysis in a tumor/pair<br /> _Optional, Default: `NA`_                                                                                                                                                                |
| `status`  | **Normal/tumor status of sample**; can be `0` (normal) or `1` (tumor).<br /> _Optional, Default: `0`_                                                                                                                                                                                                           |
| `sample`  | **Custom sample ID** for each tumor and normal sample; more than one tumor sample for each subject is possible, i.e. a tumor and a relapse; samples can have multiple lanes for which the _same_ ID must be used to merge them later (see also `lane`). Sample IDs must be unique for unique biological samples |
| `lane`    | Lane ID, used when the `sample` is multiplexed on several lanes. Must be unique for each lane in the same sample (but does not need to be the original lane name), and must contain at least one character <br /> _Required for `--step_mapping`_                                                               |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                                                                                                                      |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                                                                                                                      |
| `bam`     | Full path to (u)BAM file                                                                                                                                                                                                                                                                                        |
| `bai`     | Full path to BAM index file                                                                                                                                                                                                                                                                                     |
| `cram`    | Full path to CRAM file                                                                                                                                                                                                                                                                                          |
| `crai`    | Full path to CRAM index file                                                                                                                                                                                                                                                                                    |
| `table`   | Full path to recalibration table file                                                                                                                                                                                                                                                                           |
| `vcf`     | Full path to vcf file                                                                                                                                                                                                                                                                                           |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

#### Start with mapping (`--step mapping` [default])

This step can be started either from `fastq` files or `(u)bam`s. The `CSV` must contain at least the columns `patient`, `sample`, `lane`, and either `fastq_1/fastq_2` or `bam`.

##### Examples

Minimal config file:

```console
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_1.fastq.gz,test_2.fastq.gz
```

```console
patient,sample,lane,bam
patient1,test_sample,lane_1,test.bam
```

In this example, there are 3 read groups:

```console
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,test_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,test_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
```

```console
patient,sample,lane,bam
patient1,test_sample,1,test_L001.bam
patient1,test_sample,2,test_L002.bam
patient1,test_sample,3,test_L003.bam
```

##### Full samplesheet

In this example, all possible columns are used. There are 3 read groups for the normal sample, 2 for the tumor sample, 1 for the relapse, including the `sex` and `status` information per patient:

```console
patient,sex,status,sample,lane,fastq_1,fastq_2
patient1,XX,0,normal_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,XX,0,normal_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,XX,0,normal_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
patient1,XX,1,tumor_sample,lane_1,test2_L001_1.fastq.gz,test2_L001_2.fastq.gz
patient1,XX,1,tumor_sample,lane_2,test2_L002_1.fastq.gz,test2_L002_2.fastq.gz
patient1,XX,1,relapse_sample,lane_1,test3_L001_1.fastq.gz,test3_L001_2.fastq.gz
```

```console
patient,sex,status,sample,lane,bam
patient1,XX,0,normal_sample,lane_1,test_L001.bam
patient1,XX,0,normal_sample,lane_2,test_L002.bam
patient1,XX,0,normal_sample,lane_3,test_L003.bam
patient1,XX,1,tumor_sample,lane_1,test2_L001.bam
patient1,XX,1,tumor_sample,lane_2,test2_L002.bam
patient1,XX,1,relapse_sample,lane_1,test3_L001.bam
```

#### Start with duplicate marking (`--step markduplicates`)

##### Duplicate Marking

For starting from duplicate marking, the `CSV` file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`

Example:

```console
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```console
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The `Sarek`-generated `CSV` file is stored under `results/csv/mapped.csv` if in a previous run `--save_bam_mapped` was set and will automatically be used as an input when specifying the parameter `--step prepare_recalibration`. Otherwise this file will need to be manually generated.

##### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```console
patient,sex,status,sample,bam,bai
patient1,XX,0,test_sample,test_mapped.bam,test_mapped.bam.bai
patient1,XX,1,tumor_sample,test2_mapped.bam,test2_mapped.bam.bai
patient1,XX,1,relapse_sample,test3_mapped.bam,test3_mapped.bam.bai
```

```console
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_mapped.cram,test_mapped.cram.crai
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai
```

##### Prepare Recalibration

For starting directly from preparing recalibration, the `CSV` file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```console
patient,sample,bam,bai
patient1,test_sample,test_md.bam,test_md.bam.bai
```

```console
patient,sample,cram,crai
patient1,test_sample,test_md.cram,test_md.cram.crai
```

The `Sarek`-generated `CSV` file is stored under `results/csv/markduplicates_no_table.csv` and will automatically be used as an input when specifying the parameter `--step prepare_recalibration`.

##### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```console
patient,sex,status,sample,bam,bai
patient1,XX,0,test_sample,test_md.bam,test_md.bam.bai
patient1,XX,1,tumor_sample,test2_md.bam,test2_md.bam.bai
patient1,XX,1,relapse_sample,test3_md.bam,test3_md.bam.bai
```

```console
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_md.cram,test_md.cram.crai
patient1,XX,1,tumor_sample,test2_md.cram,test2_md.cram.crai
patient1,XX,1,relapse_sample,test3_md.cram,test3_md.cram.crai
```

#### Start with base quality score recalibration (`--step recalibrate`)

For starting from base quality score recalibration the `CSV` file must contain at least the columns `patient`, `sample`, `bam`, `bai`, `table` or `patient`, `sample`, `cram`, `crai`, `table` containing the paths to _non-recalibrated CRAM/BAM_ files and the associated recalibration table.

Example:

```console
patient,sample,bam,bai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

```console
patient,sample,cram,crai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

The `Sarek`-generated `CSV` file is stored under `results/csv/markduplicates.csv` and will automatically be used as an input when specifying the parameter `--step recalibrate`.

##### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```console
patient,sex,status,sample,cram,crai,table
patient1,XX,0,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai,test2.table
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai,test3.table
```

#### Start with variant calling (`--step variant_calling`)

For starting from the variant calling step, the `CSV` file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```console
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```console
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The `Sarek`-generated `CSV` file is stored under `results/csv/recalibrated.csv` and will automatically be used as an input when specifying the parameter `--step variant_calling`.

##### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```console
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_mapped.cram,test_mapped.cram.crai
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai
```

#### Start with annotation (`--step annotate`)

For starting from the annotation step, the `CSV` file must contain at least the columns `patient`, `sample`, `vcf`.

As `Sarek` will use `bgzip` and `tabix` to compress and index the annotated `VCF` files, it expects the input `VCF` files to be sorted.

Example:

```console
patient,sample,vcf
patient1,test_sample,test,vcf
```

The `Sarek`-generated `CSV` file is stored under `results/csv/variantcalled.csv` and will automatically be used as an input when specifying the parameter `--step annotation`.

##### Full samplesheet

In this example, all possible columns are used including the `variantcaller` information per sample:

```console
patient,sample,variantcaller,vcf
test,sample3,strelka,sample3.variants.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.diploid_sv.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.somatic_sv.vcf.gz
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/sarek
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `3.0.0`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 3.0.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

### Automatic Resubmission

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

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

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

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

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

## Tutorials

### Which tool for which data type

| Tool            | WGS | WES |  Panel | Tumor |  Normal | Somatic |
| :-------------- | :-: | :-: | :----: | :---: | :-----: | :-----: |
| Strelka2        |  x  |  x  |   x    |   x   |    x    |    x    |
| Freebayes       |  x  |  x  |   x    |   x   |    x    |    x    |
| mutect          |  x  |  x  |   x    |   x   |    -    |    x    |
| Haplotypecaller |  x  |  x  |   x    |   -   |    x    |    -    |
| Deepvariant     |  x  |  x  |   x    |   -   |    x    |    -    |
| cnvkit          |  x  |  x  |   -    |   x   |    x    |    x    |
| Msisensor       |  x  |  x  |   x    |   x   |    -    |    x    |
| controlfreec    |  x  |  x  |   x    |   x   |    -    |    x    |
| ascat           |  x  |  x  |   -    |   -   |    -    |    x    |
| manta           |  x  |  x  |   x    |   x   |    x    |    x    |
| tiddit          |  x  |  x  |   x    |   x   |    x    |    x    |

<!---
Strelka2: WGS, WES, Panel; Tumor, normal, somatic
Freebayes: WGS, WES, Panel; Tumor,normal, somatic
mutect: WGS, WES, Panel; Tumor, somatic
Haplotypecaller: WGS, WES, Panel; normal
Deepvariant: WGS, WES, Panel; normal
cnvkit: WGS, WES, Panel; Tumor, normal, somatic
Msisensor: WGS, WES, Panel; Tumor, normal, somatic
controlfreec: WGS, WES, Panel; Tumor, normal, somatic
ascat: WGS, WES, Panel; Tumor, normal, somatic
manta: WGS, WES, Panel; Tumor, normal, somatic
tiddit: WGS, WES, Panel; Tumor, normal, somatic
--->

### How to run sarek when not all reference files are in igenomes

### How to deal with a (custom) annotation cache

Examples for both snpeff and VEP

#### Containers

With `Nextflow DSL2`, each process use its own `Conda` environment or container from `biocontainers`.

For annotation, cache has to be downloaded, or specifically designed containers are available with cache.

`sareksnpeff`, our `snpeff` container is designed using [Conda](https://conda.io/).

[![sareksnpeff-docker status](https://img.shields.io/docker/automated/nfcore/sareksnpeff.svg)](https://hub.docker.com/r/nfcore/sareksnpeff)

Based on [nfcore/base:1.12.1](https://hub.docker.com/r/nfcore/base/tags), it contains:

- **[snpEff](http://snpeff.sourceforge.net/)** 4.3.1t
- Cache for `GRCh37`, `GRCh38`, `GRCm38`, `CanFam3.1` or `WBcel235`

`sarekvep`, our `vep` container is designed using [Conda](https://conda.io/).

[![sarekvep-docker status](https://img.shields.io/docker/automated/nfcore/sarekvep.svg)](https://hub.docker.com/r/nfcore/sarekvep)

Based on [nfcore/base:1.12.1](https://hub.docker.com/r/nfcore/base/tags), it contains:

- **[GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)** 1.0
- **[VEP](https://github.com/Ensembl/ensembl-vep)** 99.2
- Cache for `GRCh37`, `GRCh38`, `GRCm38`, `CanFam3.1` or `WBcel235`

#### Using downloaded cache

Both `snpEff` and `VEP` enable usage of cache.
If cache is available on the machine where `Sarek` is run, it is possible to run annotation using cache.
You need to specify the cache directory using `--snpeff_cache` and `--vep_cache` in the command lines or within configuration files.
The cache will only be used when `--annotation_cache` and cache directories are specified (either in command lines or in a configuration file).

Example:

```bash
nextflow run nf-core/sarek --tools snpEff --step annotate --sample <file.vcf.gz> --snpeff_cache </path/to/snpEff/cache> --annotation_cache
nextflow run nf-core/sarek --tools VEP --step annotate --sample <file.vcf.gz> --vep_cache </path/to/VEP/cache> --annotation_cache
```

#### Download cache

A `Nextflow` helper script has been designed to help downloading `snpEff` and `VEP` caches.
Such files are meant to be shared between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --snpeff_cache </path/to/snpEff/cache> --snpeff_db <snpEff DB version> --genome <GENOME>
nextflow run download_cache.nf --vep_cache </path/to/VEP/cache> --species <species> --vep_cache_version <VEP cache version> --genome <GENOME>
```

#### Using VEP CADD plugin

To enable the use of the `VEP` `CADD` plugin:

- Download the `CADD` files
- Specify them (either on the command line, like in the example or in a configuration file)
- use the `--cadd_cache` flag

Example:

```bash
nextflow run nf-core/sarek --step annotate --tools VEP --sample <file.vcf.gz> --cadd_cache \
    --cadd_indels </path/to/CADD/cache/InDels.tsv.gz> \
    --cadd_indels_tbi </path/to/CADD/cache/InDels.tsv.gz.tbi> \
    --cadd_wg_snvs </path/to/CADD/cache/whole_genome_SNVs.tsv.gz> \
    --cadd_wg_snvs_tbi </path/to/CADD/cache/whole_genome_SNVs.tsv.gz.tbi>
```

#### Downloading CADD files

An helper script has been designed to help downloading `CADD` files.
Such files are meant to be share between multiple users, so this script is mainly meant for people administrating servers, clusters and advanced users.

```bash
nextflow run download_cache.nf --cadd_cache </path/to/CADD/cache> --cadd_version <CADD version> --genome <GENOME>
```

### How to set sarek up to use sentieon

#### Sentieon

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

## Troubleshooting & FAQ

### Spark related issues

If you have problems running processes that make use of Spark such as `MarkDuplicates`.
You are probably experiencing issues with the limit of open files in your system.
You can check your current limit by typing the following:

```bash
ulimit -n
```

The default limit size is usually 1024 which is quite low to run Spark jobs.
In order to increase the size limit permanently you can:

Edit the file `/etc/security/limits.conf` and add the lines:

```bash
*     soft   nofile  65535
*     hard   nofile  65535
```

Edit the file `/etc/sysctl.conf` and add the line:

```bash
fs.file-max = 65535
```

Edit the file `/etc/sysconfig/docker` and add the new limits to OPTIONS like this:

```bash
OPTIONS=”—default-ulimit nofile=65535:65535"
```

Re-start your session.

Note that the way to increase the open file limit in your system may be slightly different or require additional steps.
