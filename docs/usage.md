# nf-core/sarek: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/sarek/usage](https://nf-co.re/sarek/usage)

# Introduction

Sarek is a workflow designed to detect germline and somatic variants on whole genome, whole exome, or targeted sequencing data.

Initially designed for human and mouse, it can work on any species if a reference genome is available.
Sarek is designed to handle single samples, such as single-normal or single-tumor samples, and tumor-normal pairs including additional relapses.

# Running the pipeline

## Quickstart

The typical command for running the pipeline is as follows:

```console
nextflow run nf-core/sarek --input samplesheet.csv  --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow.log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Input: Sample sheet configurations

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` to specify its location. It has to be a comma-separated file with at least 3 columns, and a header row as shown in the examples below.

It is recommended to use the absolute path of the files, but a relative path should also work.

If necessary, a tumor sample can be associated to a normal sample as a pair, if specified with the same `patient` ID, a different `sample`, and the respective `status`.
An additional tumor sample (such as a relapse for example), can be added if specified with the same `patient` ID, a different `sample`, and the `status` value `1`.

Sarek will output results in a different directory for _each sample_.
If multiple samples IDs are specified in the CSV file, Sarek will consider all files to be from different samples.

Output from Variant Calling and/or Annotation will be in a specific directory for each sample and tool configuration (or normal/tumor pair if applicable).

Multiple CSV files can be specified if the path is enclosed in quotes.

```console
--input '[path to sample sheet file(s)]'
```

### Overview: Samplesheet Columns

| Column    | Description                                                                                                                                                                                                                                                                                                                       |
| --------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `patient` | **Custom patient ID**; designates the patient/subject; must be unique for each patient, but one patient can have multiple samples (e.g. normal and tumor). <br /> _Required_                                                                                                                                                      |
| `sex`     | **Sex chromosomes of the patient**; i.e. XX, XY..., only used for Copy-Number Variation analysis in a tumor/pair<br /> _Optional, Default: `NA`_                                                                                                                                                                                  |
| `status`  | **Normal/tumor status of sample**; can be `0` (normal) or `1` (tumor).<br /> _Optional, Default: `0`_                                                                                                                                                                                                                             |
| `sample`  | **Custom sample ID** for each tumor and normal sample; more than one tumor sample for each subject is possible, i.e. a tumor and a relapse; samples can have multiple lanes for which the _same_ ID must be used to merge them later (see also `lane`). Sample IDs must be unique for unique biological samples <br /> _Required_ |
| `lane`    | Lane ID, used when the `sample` is multiplexed on several lanes. Must be unique for each lane in the same sample (but does not need to be the original lane name), and must contain at least one character <br /> _Required for `--step mapping`_                                                                                 |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`.                                                                                                                                                                                                        |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`.                                                                                                                                                                                                        |
| `bam`     | Full path to (u)BAM file                                                                                                                                                                                                                                                                                                          |
| `bai`     | Full path to BAM index file                                                                                                                                                                                                                                                                                                       |
| `cram`    | Full path to CRAM file                                                                                                                                                                                                                                                                                                            |
| `crai`    | Full path to CRAM index file                                                                                                                                                                                                                                                                                                      |
| `table`   | Full path to recalibration table file                                                                                                                                                                                                                                                                                             |
| `vcf`     | Full path to vcf file                                                                                                                                                                                                                                                                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Start with mapping (`--step mapping` [default])

This step can be started either from FastQ files or (u)BAMs. The CSV must contain at least the columns `patient`, `sample`, `lane`, and either `fastq_1/fastq_2` or `bam`.

#### Examples

Minimal config file:

```console
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_1.fastq.gz,test_2.fastq.gz
```

```console
patient,sample,lane,bam
patient1,test_sample,lane_1,test.bam
```

In this example, the sample is multiplexed over three lanes:

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

#### Full samplesheet

In this example, all possible columns are used. There are three lanes for the normal sample, two for the tumor sample, and one for the relapse sample, including the `sex` and `status` information per patient:

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

### Start with duplicate marking (`--step markduplicates`)

#### Duplicate Marking

For starting from duplicate marking, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`

> **NB:** When using [GATK4 MarkduplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/5358833264411-MarkDuplicatesSpark) reads should be name-sorted for efficient execution

Example:

```console
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```console
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/mapped.csv` if in a previous run `--save_bam_mapped` was set and will automatically be used as an input when specifying the parameter `--step markduplicates`. Otherwise this file will need to be manually generated.

#### Full samplesheet

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

### Start with preparing the recalibration tables (`--step prepare_recalibration`)

For starting directly from preparing the recalibration tables, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```console
patient,sample,bam,bai
patient1,test_sample,test_md.bam,test_md.bam.bai
```

```console
patient,sample,cram,crai
patient1,test_sample,test_md.cram,test_md.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/markduplicates_no_table.csv` and will automatically be used as an input when specifying the parameter `--step prepare_recalibration`.

#### Full samplesheet

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

### Start with base quality score recalibration (`--step recalibrate`)

For starting from base quality score recalibration the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai`, `table` or `patient`, `sample`, `cram`, `crai`, `table` containing the paths to _non-recalibrated CRAM/BAM_ files and the associated recalibration table.

Example:

```console
patient,sample,bam,bai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

```console
patient,sample,cram,crai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

The Sarek-generated CSV file is stored under `results/csv/markduplicates.csv` and will automatically be used as an input when specifying the parameter `--step recalibrate`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```console
patient,sex,status,sample,cram,crai,table
patient1,XX,0,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai,test2.table
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai,test3.table
```

### Start with variant calling (`--step variant_calling`)

For starting from the variant calling step, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```console
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```console
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/recalibrated.csv` and will automatically be used as an input when specifying the parameter `--step variant_calling`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```console
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_mapped.cram,test_mapped.cram.crai
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai
```

### Start with annotation (`--step annotate`)

For starting from the annotation step, the CSV file must contain at least the columns `patient`, `sample`, `vcf`.

As Sarek will use [bgzip](http://www.htslib.org/doc/bgzip.html) and [tabix](http://www.htslib.org/doc/tabix.html) to compress and index the annotated VCF files, it expects the input VCF files to be sorted and compressed.

Example:

```console
patient,sample,vcf
patient1,test_sample,test.vcf.gz
```

The Sarek-generated CSV file is stored under `results/csv/variantcalled.csv` and will automatically be used as an input when specifying the parameter `--step annotation`.

#### Full samplesheet

In this example, all possible columns are used including the `variantcaller` information per sample:

```console
patient,sample,variantcaller,vcf
test,sample3,strelka,sample3.variants.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.diploid_sv.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.somatic_sv.vcf.gz
```

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull nf-core/sarek
```

## Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `3.0.1`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 3.0.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

# Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

## `-profile`

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

## `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

## `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

# Custom configuration

## Resource requests

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

## Updating containers

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

## nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

# Troubleshooting & FAQ

## How to test the pipeline

When using default parameters only, sarek runs preprocessing and exits after base quality score recalibration. This is reflected in the default test profile:

```console
nextflow run nf-core/sarek -r 3.0.1 -profile test,<container/institute>
```

Expected run output:

```
[91/018ca5] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:BWAMEM1_INDEX (genome.fasta)                       [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:BWAMEM2_INDEX                                      -
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:DRAGMAP_HASHTABLE                                  -
[45/7ad672] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY (genome.fasta)      [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:MSISENSORPRO_SCAN                                  -
[79/7139ec] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:SAMTOOLS_FAIDX (genome.fasta)                      [100%] 1 of 1 ✔
[44/913bf9] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_DBSNP (dbsnp_146.hg38.vcf)                   [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_GERMLINE_RESOURCE                            -
[dc/348c16] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_KNOWN_INDELS (mills_and_1000G.indels.vcf)    [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_PON                                          -
[9f/53d6ad] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:CREATE_INTERVALS_BED (genome.interval_list)     [100%] 1 of 1 ✔
[57/a9312f] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:GATK4_INTERVALLISTTOBED (genome)                [100%] 1 of 1 ✔
[7e/b02b16] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_SPLIT (chr22_1-40001) [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                    -
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_UNMAP                -
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_MAP                  -
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_UNMAP                  -
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:SAMTOOLS_MERGE_UNMAP                     -
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:COLLATE_FASTQ_UNMAP                      -
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:COLLATE_FASTQ_MAP                        -
[-        ] process > NFCORE_SAREK:SAREK:ALIGNMENT_TO_FASTQ_INPUT:CAT_FASTQ                                -
[37/2d4ea9] process > NFCORE_SAREK:SAREK:RUN_FASTQC:FASTQC (test-test_L1)                                  [100%] 1 of 1 ✔
[a1/a64d09] process > NFCORE_SAREK:SAREK:GATK4_MAPPING:BWAMEM1_MEM (test)                                  [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:GATK4_MAPPING:BWAMEM2_MEM                                         -
[-        ] process > NFCORE_SAREK:SAREK:GATK4_MAPPING:DRAGMAP_ALIGN                                       -
[d3/488df3] process > NFCORE_SAREK:SAREK:MARKDUPLICATES:GATK4_MARKDUPLICATES (test)                        [100%] 1 of 1 ✔
[f1/0b56c6] process > NFCORE_SAREK:SAREK:MARKDUPLICATES:BAM_TO_CRAM:SAMTOOLS_BAMTOCRAM (test)              [100%] 1 of 1 ✔
[ae/e92179] process > NFCORE_SAREK:SAREK:MARKDUPLICATES:BAM_TO_CRAM:SAMTOOLS_STATS_CRAM (test)             [100%] 1 of 1 ✔
[8f/d06f35] process > NFCORE_SAREK:SAREK:MARKDUPLICATES:BAM_TO_CRAM:MOSDEPTH (test)                        [100%] 1 of 1 ✔
[38/af6ec2] process > NFCORE_SAREK:SAREK:PREPARE_RECALIBRATION:BASERECALIBRATOR (test)                     [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_RECALIBRATION:GATHERBQSRREPORTS                           -
[8b/f3ca07] process > NFCORE_SAREK:SAREK:RECALIBRATE:APPLYBQSR (test)                                      [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:RECALIBRATE:MERGE_INDEX_CRAM:MERGE_CRAM                           -
[a7/16bb3f] process > NFCORE_SAREK:SAREK:RECALIBRATE:MERGE_INDEX_CRAM:INDEX_CRAM (test)                    [100%] 1 of 1 ✔
[4d/309cb9] process > NFCORE_SAREK:SAREK:CRAM_QC:SAMTOOLS_STATS (test)                                     [100%] 1 of 1 ✔
[44/06eaf2] process > NFCORE_SAREK:SAREK:CRAM_QC:MOSDEPTH (test)                                           [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:SAMTOOLS_CRAMTOBAM_RECAL                                          [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:CUSTOM_DUMPSOFTWAREVERSIONS                                       [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:MULTIQC                                                           [100%] 1 of 1 ✔
```

The pipeline comes with a number of possible paths and tools that can be used. The easiest and fastest test to see that the preprocessing + variantcalling (in this case Strelka2) works, is to run:

```console
nextflow run nf-core/sarek -r 3.0.1 -profile test,<container/institute> --tools strelka
```

Due to the small test data size, unfortunately not everything can be tested from top-to-bottom, but often is done by utilizing the pipeline's `--step` parameter. Annotation has to tested separatly from the remaining workflow, since we use references for `C.elegans`, while the remaining tests are run on downsampled human data.

```console
nextflow run nf-core/sarek -r 3.0.1 -profile test,<container/institute> --tools snpeff --step annotation
```

If you are interested in any of the other tests that are run on every code change or would like to run them yourself, you can take a look at `tests/<filename>.yml`. For each entry the respective nextflow command run and the expected output is specified.

Some of the currently, available test profiles:

| Test profile    | Run command                                                                     |
| :-------------- | :------------------------------------------------------------------------------ |
| annotation      | `nextflow run main.nf -profile test,annotation,docker --tools snpeff.vep,merge` |
| no_intervals    | `nextflow run main.nf -profile test,no_intervals,docker`                        |
| targeted        | `nextflow run main.nf -profile test,targeted,docker`                            |
| tools_germline  | `nextflow run main.nf -profile test,tools_germline,docker --tools strelka`      |
| tools_tumoronly | `nextflow run main.nf -profile test,tools_tumoronly,docker --tools strelka`     |
| tools_somatic   | `nextflow run main.nf -profile test,tools_somatic,docker --tools strelka`       |
| trimming        | `nextflow run main.nf -profile test,trim_fastq,docker`                          |
| umi             | `nextflow run main.nf -profile test,umi,docker`                                 |
| use_gatk_spark  | `nextflow run main.nf -profile test,use_gatk_spark,docker`                      |

## How can the different steps be used

Sarek can be started at different points in the analysis by setting the parameter `--step`. Once started at a certain point, the pipeline runs through all the following steps without additional intervention. For example when starting from `--step mapping` (set by default) and `--tools strelka,vep`, the input reads will be aligned, duplicate marked, recalibrated, variant called with Strelka, and finally VEP will annotate the called variants.

## Which variant calling tool is implemented for which data type?

This list is by no means exhaustive and it will depend on the specific analysis you would like to run. This is a suggestion based on the individual docs of the tools specifically for human genomes and a garden-variety sequencing run as well as what has been added to the pipeline.

| Tool                                                                                                    | WGS | WES |  Panel |  Normal | Tumor | Somatic |
| :------------------------------------------------------------------------------------------------------ | :-: | :-: | :----: | :-----: | :---: | :-----: |
| [DeepVariant](https://github.com/google/deepvariant)                                                    |  x  |  x  |   x    |    x    |   -   |    -    |
| [FreeBayes](https://github.com/ekg/freebayes)                                                           |  x  |  x  |   x    |    x    |   x   |    x    |
| [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller) |  x  |  x  |   x    |    x    |   -   |    -    |
| [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2)                 |  x  |  x  |   x    |    -    |   x   |    x    |
| [mpileup](https://www.htslib.org/doc/samtools-mpileup.html)                                             |  x  |  x  |   x    |    x    |   x   |    -    |
| [Strelka2](https://github.com/Illumina/strelka)                                                         |  x  |  x  |   x    |    x    |   x   |    x    |
| [Manta](https://github.com/Illumina/manta)                                                              |  x  |  x  |   x    |    x    |   x   |    x    |
| [TIDDIT](https://github.com/SciLifeLab/TIDDIT)                                                          |  x  |  x  |   x    |    x    |   x   |    x    |
| [ASCAT](https://github.com/VanLoo-lab/ascat)                                                            |  x  |  x  |   -    |    -    |   -   |    x    |
| [CNVKit](https://cnvkit.readthedocs.io/en/stable/)                                                      |  x  |  x  |   -    |    x    |   x   |    x    |
| [Control-FREEC](https://github.com/BoevaLab/FREEC)                                                      |  x  |  x  |   x    |    -    |   x   |    x    |
| [MSIsensorPro](https://github.com/xjtu-omics/msisensor-pro)                                             |  x  |  x  |   x    |    -    |   -   |    x    |

## How to run ASCAT with whole-exome sequencing data?

While the ASCAT implementation in sarek is capable of running with whole-exome sequencing data, the needed references are currently not provided with the igenomes.config. According to the [developers](https://github.com/VanLoo-lab/ascat/issues/97) of ASCAT, loci and allele files (one file per chromosome) can be downloaded directly from the [Battenberg repository](https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52).

The GC correction file needs to be derived, so one has to concatenate all chromosomes into a single file and modify the header so it fits [this example](https://github.com/VanLoo-lab/ascat/tree/master/LogRcorrection#gc-correction-file-creation).

The RT correction file is missing for hg38 but can be derived using [ASCAT scripts](https://github.com/VanLoo-lab/ascat/tree/master/LogRcorrection#replication-timing-correction-file-creation) for hg19. For hg38, one needs to lift-over hg38 to hg19, run the script on hg19 positions and set coordinates back to hg38.

Please note that:

Row names (for GC and RT correction files) should be `${chr}_${position}` (there is no SNP/probe ID for HTS data).
ASCAT developers strongly recommend using a BED file for WES/TS data. This prevents considering SNPs covered by off-targeted reads that would add noise to log/BAF tracks.

## What are the bwa/bwa-mem2 parameters?

For mapping, sarek follows the parameter suggestions provided in this [paper](https://www.nature.com/articles/s41467-018-06159-4):

`-K 100000000` : for deterministic pipeline results, for more info see [here](https://github.com/CCDG/Pipeline-Standardization/issues/2)

`-Y`: force soft-clipping rather than default hard-clipping of supplementary alignments

In addition, currently the mismatch penalty for reads with tumor status in the sample sheet are mapped with a mismatch penalty of `-B 3`.

## How to create a panel-of-normals for Mutect2

For a detailed tutorial on how to create a panel-of-normals, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132).

## Spark related issues

If you have problems running processes that make use of Spark such as `MarkDuplicates`.
You are probably experiencing issues with the limit of open files in your system.
You can check your current limit by typing the following:

```console
ulimit -n
```

The default limit size is usually 1024 which is quite low to run Spark jobs.
In order to increase the size limit permanently you can:

Edit the file `/etc/security/limits.conf` and add the lines:

```console
*     soft   nofile  65535
*     hard   nofile  65535
```

Edit the file `/etc/sysctl.conf` and add the line:

```console
fs.file-max = 65535
```

Edit the file `/etc/sysconfig/docker` and add the new limits to OPTIONS like this:

```console
OPTIONS=”—default-ulimit nofile=65535:65535"
```

Re-start your session.

Note that the way to increase the open file limit in your system may be slightly different or require additional steps.

### Cannot delete work folder when using docker + Spark

Currently, when running spark-based tools in combination with docker, it is required to set `docker.userEmulation = false`. This can unfortunately causes permission issues when `work/` is being written with root permissions. In case this happens, you might need to configure docker to run without `userEmulation` (see [here](https://github.com/Midnighter/nf-core-adr/blob/main/docs/adr/0008-refrain-from-using-docker-useremulation-in-nextflow.md)).

## How to handle UMIs

Sarek can process UMI-reads, using [fgbio](http://fulcrumgenomics.github.io/fgbio/tools/latest/) tools.

In order to use reads containing UMI tags as your initial input, you need to include `--umi_read_structure [structure]` in your parameters.

This will enable pre-processing of the reads and UMI consensus reads calling, which will then be used to continue the workflow from the mapping steps. For post-UMI processing depending on the experimental setup, duplicate marking and base quality recalibration can be skipped with [`--skip_tools`].

### UMI Read Structure

This parameter is a string, which follows a [convention](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) to describe the structure of the umi.

As an example: if your reads contain a UMI only on the forward read, the string can only represent one structure (i.e. "2M11S+T"); should your reads contain a UMI on both reas, the string will contain two structures separated by a blank space (i.e. "2M11S+T 2M11S+T"); should your reads contain a UMI only on the reverse read, your structure must represent the template only for the forward read and template plus UMI for the reverse read (i.e. +T 12M11S+T). Please do refer to FGBIO documentation for more details, as providing the correct structure is essential and specific to the UMI kit used.

### Limitations and future updates

Recent updates to Samtools have been introduced, which can speed-up performance of fgbio tools used in this workflow.
The current workflow does not handle duplex UMIs (i.e. where opposite strands of a duplex molecule have been tagged with a different UMI), and best practices have been proposed to process this type of data.
Both changes will be implemented in a future release.

## How to run sarek when no(t all) reference files are in igenomes

For common genomes, such as GRCh38 and GRCh37, the pipeline is shipped with (almost) all necessary reference files. However, sometimes it is necessary to use custom references for some or all files:

### No igenomes reference files are used

If none of your required genome files are in igenomes, `--igenomes_ignore` must be set to ignore any igenomes input and `--genome null`. The `fasta` file is the only required input file and must be provided to run the pipeline. All other possible reference file can be provided in addition. For details, see the paramter documentation.

Minimal example for custom genomes:

```console
nextflow run nf-core/sarek --genome null --igenomes_ignore --fasta <custom.fasta>
```

### Overwrite specific reference files

If you don't want to use some of the provided reference genomes, they can be overwritten by either providing a new file or setting the respective file parameter to `false`, if it should be ignored:

Example for using a custom known indels file:

```console
nextflow run nf-core/sarek --known_indels <my_known_indels.vcf.gz> --genome GRCh38.GATK
```

Example for not using known indels, but all other provided reference file:

```console
nextflow run nf-core/sarek --known_indels false --genome GRCh38.GATK
```

### Where do the used reference genomes originate from

For GATK.GRCh38 the links for each reference file and the corresponding processes that use them is listed below. For GATK.GRCh37 the files originate from the same sources:

| File                  | Tools                                                                                                                                                                                                                                                                                                                                                                                                                                                | Origin                                                                                                                | Docs                                                                                 |
| :-------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------- |
| ascat_alleles         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/uouszfktzgoqfy7/G1000_alleles_hg38.zip                                                      | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci            | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/80cq0qgao8l1inj/G1000_loci_hg38.zip                                                         | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci_gc         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/80cq0qgao8l1inj/G1000_loci_hg38.zip                                                         | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci_rt         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                | https://www.dropbox.com/s/xlp99uneqh6nh6p/RT_G1000_hg38.zip                                                           | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| bwa                   | bwa-mem                                                                                                                                                                                                                                                                                                                                                                                                                                              | bwa index -p bwa/${fasta.baseName} $fasta                                                                             |                                                                                      |
| bwamem2               | bwa-mem2                                                                                                                                                                                                                                                                                                                                                                                                                                             | bwa-mem2 index -p bwamem2/${fasta} $fasta                                                                             |                                                                                      |
| dragmap               | DragMap                                                                                                                                                                                                                                                                                                                                                                                                                                              | dragen-os --build-hash-table true --ht-reference $fasta --output-directory dragmap                                    |                                                                                      |
| dbsnp                 | Baserecalibrator, ControlFREEC, GenotypeGVCF, HaplotypeCaller                                                                                                                                                                                                                                                                                                                                                                                        | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| dbsnp_tbi             | Baserecalibrator, ControlFREEC, GenotypeGVCF, HaplotypeCaller                                                                                                                                                                                                                                                                                                                                                                                        | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| dict                  | Baserecalibrator(Spark), CNNScoreVariant, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, MarkDulpicates(Spark), MergeVCFs, Mutect2, Variantrecalibrator                                                                                                                                                                                               | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| fasta                 | ApplyBQSR(Spark), ApplyVQSR, ASCAT, Baserecalibrator(Spark), BWA, BWAMem2, CNNScoreVariant, CNVKit, ControlFREEC, DragMap, DEEPVariant, EnsemblVEP, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, FreeBayes, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, interval building, Manta, MarkDuplicates(Spark),MergeVCFs,MSISensorPro, Mutect2, Samtools, snpEff, Strelka, Tiddit, Variantrecalibrator | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| fasta_fai             | ApplyBQSR(Spark), ApplyVQSR, ASCAT, Baserecalibrator(Spark), BWA, BWAMem2, CNNScoreVariant, CNVKit, ControlFREEC, DragMap, DEEPVariant, EnsemblVEP, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, FreeBayes, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, interval building, Manta, MarkDuplicates(Spark),MergeVCFs,MSISensorPro, Mutect2, Samtools, snpEff, Strelka, Tiddit, Variantrecalibrator | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| germline_resource     | GetPileupsummaries,Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                           | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| germline_resource_tbi | GetPileupsummaries,Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                           | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| intervals             | ApplyBQSR(Spark), ASCAT, Baserecalibrator(Spark), BCFTools, CNNScoreVariants, ControlFREEC, Deepvariant, FilterVariantTranches, FreeBayes, GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, Strelka, mpileup, MSISensorPro, Mutect2, VCFTools                                                                                                                                                                                                      | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_indels          | BaseRecalibrator(Spark), FilterVariantTranches                                                                                                                                                                                                                                                                                                                                                                                                       | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_indels_tbi      | BaseRecalibrator(Spark), FilterVariantTranches                                                                                                                                                                                                                                                                                                                                                                                                       | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_snps            | BaseRecalibrator(Spark), FilterVariantTranches, VariantRecalibrator                                                                                                                                                                                                                                                                                                                                                                                  | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_snps_tbi        | BaseRecalibrator(Spark), FilterVariantTranches, VariantRecalibrator                                                                                                                                                                                                                                                                                                                                                                                  | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |
| mappability           | ControlFREEC                                                                                                                                                                                                                                                                                                                                                                                                                                         | http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip                                                                | http://boevalab.inf.ethz.ch/FREEC/tutorial.html                                      |
| pon                   | Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                              | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- |
| pon_tbi               | Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                              | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- |

## How to customise SnpEff and VEP annotation

Sarek uses nf-core provided containers for both snpEff and VEP for several reference genomes ('CanFam3', 'GRCh37', 'GRCh38', 'GRCm38' and 'WBcel235').

### Using downloaded cache

Both `snpEff` and `VEP` enable usage of cache, if no pre-build container is available.
The cache needs to be made available on the machine where Sarek is run.
You need to specify the cache directory using `--snpeff_cache` and `--vep_cache` in the command lines or within configuration files.

Example:

```console
nextflow run nf-core/sarek --tools snpEff --step annotate --sample <file.vcf.gz> --snpeff_cache </path/to/snpEff/cache>
nextflow run nf-core/sarek --tools VEP --step annotate --sample <file.vcf.gz> --vep_cache </path/to/VEP/cache>
```

Similarly, when wanting to use a different cache than the one specified in the iGenomes config file, one can use `--snpeff_db`, `--snpeff_genome`, `--snpeff_version`, `--vep_cache_version`, `--vep_genome`, `--vep_species` and `--vep_version` to overwrite these default value related to the databases, genomes, versions and caches' versions used by these tools.

### Using VEP plugins

#### dbnsfp

Enable with `--vep_dbnsfp`. The following parameters are mandatory:

- `--dbnsfp`, to specify the path to the dbNSFP processed file.
- `--dbnsfp_tbi`, to specify the path to the dbNSFP tabix indexed file.

The following parameters are optionnal:

- `--dbnsfp_consequence`, to filter/limit outputs to a specific effect of the variant.
  - The set of consequence terms is defined by the Sequence Ontology and an overview of those used in VEP can be found [here](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html).
  - If one wants to filter using several consequences, then separate those by using '&' (i.e. `--dbnsfp_consequence '3_prime_UTR_variant&intron_variant'`.",
- `--dbnsfp_fields`, to retrieve individual values from the dbNSFP file.
  - The values correspond to the name of the columns in the dbNSFP file and are separated by comma.
  - The column names might differ between the different dbNSFP versions. Please check the Readme.txt file, which is provided with the dbNSFP file, to obtain the correct column names. The Readme file contains also a short description of the provided values and the version of the tools used to generate them.

For more details, see [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp).

#### LOFTEE

Enable with `--vep_loftee`.

For more details, see [here](https://github.com/konradjk/loftee).

#### SpliceAi

Enable with `--vep_spliceai`. The following parameters are mandatory:

- `--spliceai_snv`, to specify the path to SpliceAI raw scores snv file.
- `--spliceai_snv_tbi`, to specify the path to SpliceAI raw scores snv tabix indexed file.
- `--spliceai_indel`, to specify the path to SpliceAI raw scores indel file.
- `--spliceai_indel_tbi`, to specify the path to SpliceAI raw scores indel tabix indexed file.

For more details, see [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceai).

#### SpliceRegions

Enable with `--vep_spliceregion`.

For more details, see [here](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#spliceregion) and [here](https://www.ensembl.info/2018/10/26/cool-stuff-the-vep-can-do-splice-site-variant-annotation/)."

## Requested resources for the tools

Resource requests are difficult to generalize and are often dependent on input data size. Currently, the number of cpus and memory requested by default were adapted from tests on 5 ICGC paired whole-genome sequencing samples with approximately 40X and 80X depth.
For targeted data analysis, this is overshooting by a lot. In this case resources for each process can be limited by either setting `--max_memory` and `-max_cpus` or tailoring the request by process name as described [here](#resource-requests). If you are using sarek for a certain data type regulary, and would like to make these requests available to others on your system, an institution-specific, pipeline-specific config file can be added [here](https://github.com/nf-core/configs/tree/master/conf/pipeline/sarek).

## MultiQC related issues

### Plots for SnpEff are missing

When plots are missing, it is possible that the fasta and the custom SnpEff database are not matching https://pcingola.github.io/SnpEff/se_faq/#error_chromosome_not_found-details.
The SnpEff completes without throwing an error causing nextflow to complete successfully. An indication for the error are these lines in the `.command` files:

```text
ERRORS: Some errors were detected
Error type      Number of errors
ERROR_CHROMOSOME_NOT_FOUND      17522411
```

## How to set up sarek to use sentieon

Sarek 3.0.1 is currently not supporting sentieon. It is planned for the upcoming release 3.1. In the meantime, please revert to the last release 2.7.2.
