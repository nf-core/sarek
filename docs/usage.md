# nf-core/sarek: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/sarek/usage](https://nf-co.re/sarek/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

# Introduction

Sarek is a workflow designed to detect germline and somatic variants on whole genome, whole exome, or targeted sequencing data.

Initially designed for human and mouse, it can work on any species if a reference genome is available.
Sarek is designed to handle single samples, such as single-normal or single-tumor samples, and tumor-normal pairs including additional relapses.

# Running the pipeline

## Quickstart

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/sarek -r <VERSION> -profile <PROFILE> --input ./samplesheet.csv --outdir ./results --genome GATK.GRCh38 --tools <TOOLS>
```

`-r <VERSION>` is optional but strongly recommended for reproducibility and should match the latest version.

`-profile <PROFILE>` is mandatory and should reflect either your own institutional profile or any pipeline profile specified in the [profile section](##-profile).

This documentation imply that any `nextflow run nf-core/sarek` command is run with the appropriate `-r` and `-profile` commands.

This will launch the pipeline and perform variant calling with the tools specified in `--tools`, see the [parameter section](https://nf-co.re/sarek/latest/parameters#tools) for details on variant calling tools.
In the above example the pipeline runs with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/sarek -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GATK.GRCh38'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Input: Sample sheet configurations

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` to specify its location. It has to be a comma-separated file with at least 3 columns, and a header row as shown in the examples below.

It is recommended to use the absolute path of the files, but a relative path should also work.

If necessary, a tumor sample can be associated to a normal sample as a pair, if specified with the same `patient` ID, a different `sample`, and the respective `status`.
An additional tumor sample (such as a relapse for example), can be added if specified with the same `patient` ID, a different `sample`, and the `status` value `1`.

Sarek will output results in a different directory for _each sample_.
If multiple samples IDs are specified in the CSV file, Sarek will consider all files to be from different samples.

Output from Variant Calling and/or Annotation will be in a specific directory for each sample and tool configuration (or normal/tumor pair if applicable).

### Overview: Samplesheet Columns

| Column     | Description                                                                                                                                                                                                                                                                                                                       |
| ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `patient`  | **Custom patient ID**; designates the patient/subject; must be unique for each patient, but one patient can have multiple samples (e.g. normal and tumor). <br /> _Required_                                                                                                                                                      |
| `sex`      | **Sex chromosomes of the patient**; i.e. XX, XY..., only used for Copy-Number Variation analysis in a tumor/pair<br /> _Optional, Default: `NA`_                                                                                                                                                                                  |
| `status`   | **Normal/tumor status of sample**; can be `0` (normal) or `1` (tumor).<br /> _Optional, Default: `0`_                                                                                                                                                                                                                             |
| `sample`   | **Custom sample ID** for each tumor and normal sample; more than one tumor sample for each subject is possible, i.e. a tumor and a relapse; samples can have multiple lanes for which the _same_ ID must be used to merge them later (see also `lane`). Sample IDs must be unique for unique biological samples <br /> _Required_ |
| `lane`     | Lane ID, used when the `sample` is multiplexed on several lanes. Must be unique for each lane in the same sample (but does not need to be the original lane name), and must contain at least one character <br /> _Required for `--step mapping`_                                                                                 |
| `fastq_1`  | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`.                                                                                                                                                                                                        |
| `fastq_2`  | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`.                                                                                                                                                                                                        |
| `spring_1` | Full path to spring-compressed, gzipped FastQ file for read 1 or for reads 1 and 2. The Fastq file has to be first gzipped, then spring-compressed, and it must have the extension `.fastq.gz.spring` or `.fq.gz.spring`.                                                                                                         |
| `spring_2` | Full path to spring-compressed, gzipped FastQ file for read 2. The Fastq file has to be first gzipped, then spring-compressed, and it must have the extension `.fastq.gz.spring` or `.fq.gz.spring`.                                                                                                                              |
| `bam`      | Full path to (u)BAM file                                                                                                                                                                                                                                                                                                          |
| `bai`      | Full path to BAM index file                                                                                                                                                                                                                                                                                                       |
| `cram`     | Full path to CRAM file                                                                                                                                                                                                                                                                                                            |
| `crai`     | Full path to CRAM index file                                                                                                                                                                                                                                                                                                      |
| `table`    | Full path to recalibration table file                                                                                                                                                                                                                                                                                             |
| `vcf`      | Full path to vcf file                                                                                                                                                                                                                                                                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Start with mapping (`--step mapping` [default])

This step can be started either from FastQ files (gzip-compressed or gzip+spring-compressed) or (u)BAMs. The CSV must contain at least the columns `patient`, `sample`, `lane`, and `fastq_1/fastq_2`, `spring_1`, `spring_1/spring_2` or `bam`.

#### Examples

Minimal config file:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_1.fastq.gz,test_2.fastq.gz
```

```bash
patient,sample,lane,spring_1
patient1,test_sample,lane_1,test_R1_and_R2.fastq.gz.spring
```

```bash
patient,sample,lane,spring_1,spring_2
patient1,test_sample,lane_1,test_R1.fastq.gz.spring,test_R2.fastq.gz.spring
```

```bash
patient,sample,lane,bam
patient1,test_sample,lane_1,test.bam
```

In this example, the sample is multiplexed over three lanes:

```bash
patient,sample,lane,fastq_1,fastq_2
patient1,test_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,test_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,test_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
```

```bash
patient,sample,lane,bam
patient1,test_sample,1,test_L001.bam
patient1,test_sample,2,test_L002.bam
patient1,test_sample,3,test_L003.bam
```

#### Full samplesheet

In this example, all possible columns are used. There are three lanes for the normal sample, two for the tumor sample, and one for the relapse sample, including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,lane,fastq_1,fastq_2
patient1,XX,0,normal_sample,lane_1,test_L001_1.fastq.gz,test_L001_2.fastq.gz
patient1,XX,0,normal_sample,lane_2,test_L002_1.fastq.gz,test_L002_2.fastq.gz
patient1,XX,0,normal_sample,lane_3,test_L003_1.fastq.gz,test_L003_2.fastq.gz
patient1,XX,1,tumor_sample,lane_1,test2_L001_1.fastq.gz,test2_L001_2.fastq.gz
patient1,XX,1,tumor_sample,lane_2,test2_L002_1.fastq.gz,test2_L002_2.fastq.gz
patient1,XX,1,relapse_sample,lane_1,test3_L001_1.fastq.gz,test3_L001_2.fastq.gz
```

```bash
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

```bash
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```bash
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/mapped.csv` if in a previous run `--save_mapped` was set and will automatically be used as an input when specifying the parameter `--step markduplicates`. Otherwise this file will need to be manually generated.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,bam,bai
patient1,XX,0,test_sample,test_mapped.bam,test_mapped.bam.bai
patient1,XX,1,tumor_sample,test2_mapped.bam,test2_mapped.bam.bai
patient1,XX,1,relapse_sample,test3_mapped.bam,test3_mapped.bam.bai
```

```bash
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_mapped.cram,test_mapped.cram.crai
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai
```

### Start with preparing the recalibration tables (`--step prepare_recalibration`)

For starting directly from preparing the recalibration tables, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```bash
patient,sample,bam,bai
patient1,test_sample,test_md.bam,test_md.bam.bai
```

```bash
patient,sample,cram,crai
patient1,test_sample,test_md.cram,test_md.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/markduplicates_no_table.csv` and will automatically be used as an input when specifying the parameter `--step prepare_recalibration`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,bam,bai
patient1,XX,0,test_sample,test_md.bam,test_md.bam.bai
patient1,XX,1,tumor_sample,test2_md.bam,test2_md.bam.bai
patient1,XX,1,relapse_sample,test3_md.bam,test3_md.bam.bai
```

```bash
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_md.cram,test_md.cram.crai
patient1,XX,1,tumor_sample,test2_md.cram,test2_md.cram.crai
patient1,XX,1,relapse_sample,test3_md.cram,test3_md.cram.crai
```

### Start with base quality score recalibration (`--step recalibrate`)

For starting from base quality score recalibration the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai`, `table` or `patient`, `sample`, `cram`, `crai`, `table` containing the paths to _non-recalibrated CRAM/BAM_ files and the associated recalibration table.

Example:

```bash
patient,sample,bam,bai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

```bash
patient,sample,cram,crai,table
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
```

The Sarek-generated CSV file is stored under `results/csv/markduplicates.csv` and will automatically be used as an input when specifying the parameter `--step recalibrate`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,cram,crai,table
patient1,XX,0,test_sample,test_mapped.cram,test_mapped.cram.crai,test.table
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai,test2.table
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai,test3.table
```

### Start with variant calling (`--step variant_calling`)

For starting from the variant calling step, the CSV file must contain at least the columns `patient`, `sample`, `bam`, `bai` or `patient`, `sample`, `cram`, `crai`.

Example:

```bash
patient,sample,bam,bai
patient1,test_sample,test_mapped.bam,test_mapped.bam.bai
```

```bash
patient,sample,cram,crai
patient1,test_sample,test_mapped.cram,test_mapped.cram.crai
```

The Sarek-generated CSV file is stored under `results/csv/recalibrated.csv` and will automatically be used as an input when specifying the parameter `--step variant_calling`.

#### Full samplesheet

In this example, all possible columns are used including the `sex` and `status` information per patient:

```bash
patient,sex,status,sample,cram,crai
patient1,XX,0,normal_sample,test_mapped.cram,test_mapped.cram.crai
patient1,XX,1,tumor_sample,test2_mapped.cram,test2_mapped.cram.crai
patient1,XX,1,relapse_sample,test3_mapped.cram,test3_mapped.cram.crai
```

### Start with annotation (`--step annotate`)

For starting from the annotation step, the CSV file must contain at least the columns `patient`, `sample`, `vcf`.

As Sarek will use [bgzip](http://www.htslib.org/doc/bgzip.html) and [tabix](http://www.htslib.org/doc/tabix.html) to compress and index the annotated VCF files, it expects the input VCF files to be sorted and compressed.

Example:

```bash
patient,sample,vcf
patient1,test_sample,test.vcf.gz
```

The Sarek-generated CSV file is stored under `results/csv/variantcalled.csv` and will automatically be used as an input when specifying the parameter `--step annotation`.

#### Full samplesheet

In this example, all possible columns are used including the `variantcaller` information per sample:

```bash
patient,sample,variantcaller,vcf
test,sample3,strelka,sample3.variants.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.diploid_sv.vcf.gz
test,sample4_vs_sample3,manta,sample4_vs_sample3.somatic_sv.vcf.gz
```

## Updating the pipeline

When you launch a pipeline from the command-line with `nextflow run nf-core/sarek -params-file params.yaml`, Nextflow will automatically pull the pipeline code from GitHub and store it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/sarek
```

## Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `3.3.2`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 3.3.2`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

# Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

## `-profile`

Use this parameter to choose a configuration profile.
Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time.
For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
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
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

## `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

## `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
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

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

## nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

# Troubleshooting & FAQ

## How to test the pipeline

When using default parameters only, sarek runs preprocessing and `Strelka`.
This is reflected in the default test profile:

```bash
nextflow run nf-core/sarek -profile test,<container/institute> --outdir results
```

Expected run output:

```bash
[85/6b7739] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:BWAMEM1_INDEX (genome.fasta)                                                [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:BWAMEM2_INDEX                                                               -
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:DRAGMAP_HASHTABLE                                                           -
[22/cf54a8] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:GATK4_CREATESEQUENCEDICTIONARY (genome.fasta)                               [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:MSISENSORPRO_SCAN                                                           -
[28/dad25a] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:SAMTOOLS_FAIDX (genome.fasta)                                               [100%] 1 of 1 ✔
[23/3fe964] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_DBSNP (dbsnp_146.hg38.vcf)                                            [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_GERMLINE_RESOURCE                                                     -
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_KNOWN_SNPS                                                            -
[14/26e286] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_KNOWN_INDELS (mills_and_1000G.indels.vcf)                             [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:PREPARE_GENOME:TABIX_PON                                                                   -
[76/04d107] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:CREATE_INTERVALS_BED (genome.interval_list)                              [100%] 1 of 1 ✔
[d4/f97174] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:GATK4_INTERVALLISTTOBED (genome)                                         [100%] 1 of 1 ✔
[70/82ba3c] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_SPLIT (chr22_1-40001)                          [100%] 1 of 1 ✔
[d4/c2d0c4] process > NFCORE_SAREK:SAREK:PREPARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_COMBINED (genome)                              [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                                                  -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_UNMAP                                              -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_UNMAP_MAP                                                -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_UNMAP                                                -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_MERGE_UNMAP                                                   -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_UNMAP                                                    -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:COLLATE_FASTQ_MAP                                                      -
[-        ] process > NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:CAT_FASTQ                                                              -
[c4/f59e5a] process > NFCORE_SAREK:SAREK:FASTQC (test-test_L1)                                                                      [100%] 1 of 1 ✔
[0b/c5a999] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP:BWAMEM1_MEM (test)                                         [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP:BWAMEM2_MEM                                                -
[-        ] process > NFCORE_SAREK:SAREK:FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP:DRAGMAP_ALIGN                                              -
[c7/664cd1] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (test)                                             [100%] 1 of 1 ✔
[13/bc73b6] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:INDEX_MARKDUPLICATES (test)                                             [100%] 1 of 1 ✔
[2a/99608e] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS (test)                         [100%] 1 of 1 ✔
[f2/0420ca] process > NFCORE_SAREK:SAREK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH (test)                               [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:CRAM_TO_BAM                                                                                -
[eb/46945a] process > NFCORE_SAREK:SAREK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (test)                                         [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:BAM_BASERECALIBRATOR:GATK4_GATHERBQSRREPORTS                                               -
[ec/2377d4] process > NFCORE_SAREK:SAREK:BAM_APPLYBQSR:GATK4_APPLYBQSR (test)                                                       [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:MERGE_CRAM                                         -
[88/3af664] process > NFCORE_SAREK:SAREK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:INDEX_CRAM (test)                                  [100%] 1 of 1 ✔
[f4/828fde] process > NFCORE_SAREK:SAREK:CRAM_QC_RECAL:SAMTOOLS_STATS (test)                                                        [100%] 1 of 1 ✔
[fb/a9d66f] process > NFCORE_SAREK:SAREK:CRAM_QC_RECAL:MOSDEPTH (test)                                                              [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:CRAM_TO_BAM_RECAL                                                                          -
[ef/026185] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:STRELKA_SINGLE (test)  [100%] 1 of 1 ✔
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA          -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_GERMLINE_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA_GENOME   -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_TUMOR_ONLY_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:STRELKA_SINGLE       -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_TUMOR_ONLY_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA        -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_TUMOR_ONLY_ALL:BAM_VARIANT_CALLING_SINGLE_STRELKA:MERGE_STRELKA_GENOME -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_SOMATIC_ALL:BAM_VARIANT_CALLING_SOMATIC_STRELKA:STRELKA_SOMATIC        -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_SOMATIC_ALL:BAM_VARIANT_CALLING_SOMATIC_STRELKA:MERGE_STRELKA_INDELS   -
[-        ] process > NFCORE_SAREK:SAREK:BAM_VARIANT_CALLING_SOMATIC_ALL:BAM_VARIANT_CALLING_SOMATIC_STRELKA:MERGE_STRELKA_SNVS     -
[bc/f3f5cf] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:BCFTOOLS_STATS (test)                                             [100%] 1 of 1 ✔
[21/8d4f02] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_COUNT (test)                                        [100%] 1 of 1 ✔
[36/957fba] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_QUAL (test)                                         [100%] 1 of 1 ✔
[70/a8e064] process > NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_SUMMARY (test)                                           [100%] 1 of 1 ✔
[36/e35b1b] process > NFCORE_SAREK:SAREK:CUSTOM_DUMPSOFTWAREVERSIONS (1)                                                            [100%] 1 of 1 ✔
[3f/3c3356] process > NFCORE_SAREK:SAREK:MULTIQC                                                                                    [100%] 1 of 1 ✔
-[nf-core/sarek] Pipeline completed successfully-
Completed at: 09-Jun-2023 13:46:31
Duration    : 1m 50s
CPU hours   : (a few seconds)
Succeeded   : 27
```

The pipeline comes with a number of possible paths and tools that can be used.

Due to the small test data size, unfortunately not everything can be tested from top-to-bottom, but often is done by utilizing the pipeline's `--step` parameter.

For more extensive testing purpose, we have the `test_cache` profile that contain the same data, but on which the path to the reference and input files can be changed using the `--test_data_base` params.

Annotation is generally tested separately from the remaining workflow, since we use references for `C.elegans`, while the remaining tests are run on downsampled human data.

```bash
nextflow run nf-core/sarek -profile test_cache,<container/institute> --outdir results --tools snpeff --step annotation
```

If you are interested in any of the other tests that are run on every code change or would like to run them yourself, you can take a look at `tests/<filename>.yml`.
For each entry the respective nextflow command run and the expected output is specified.

Some of the currently, available test profiles:

| Test profile    | Run command                                                                           |
| :-------------- | :------------------------------------------------------------------------------------ |
| annotation      | `nextflow run main.nf -profile test_cache,annotation,docker --tools snpeff,vep,merge` |
| no_intervals    | `nextflow run main.nf -profile test_cache,no_intervals,docker`                        |
| targeted        | `nextflow run main.nf -profile test_cache,targeted,docker`                            |
| tools_germline  | `nextflow run main.nf -profile test_cache,tools_germline,docker --tools strelka`      |
| tools_tumoronly | `nextflow run main.nf -profile test_cache,tools_tumoronly,docker --tools mutect2`     |
| tools_somatic   | `nextflow run main.nf -profile test_cache,tools_somatic,docker --tools strelka`       |
| trimming        | `nextflow run main.nf -profile test_cache,trim_fastq,docker`                          |
| umi             | `nextflow run main.nf -profile test_cache,umi,docker`                                 |
| use_gatk_spark  | `nextflow run main.nf -profile test_cache,use_gatk_spark,docker`                      |

If you are interested in any of the other profiles that are used, you can take a look at `conf/test/<filename>.config`.

## How can the different steps be used

Sarek can be started at different points in the analysis by setting the parameter `--step`. Once started at a certain point, the pipeline runs through all the following steps without additional intervention. For example when starting from `--step mapping` (set by default) and `--tools strelka,vep`, the input reads will be aligned, duplicate marked, recalibrated, variant called with Strelka, and finally VEP will annotate the called variants.

## Which variant calling tool is implemented for which data type?

This list is by no means exhaustive and it will depend on the specific analysis you would like to run. This is a suggestion based on the individual docs of the tools specifically for human genomes and a garden-variety sequencing run as well as what has been added to the pipeline.

| Tool                                                                                                    | WGS | WES |  Panel |  Germline | Tumor-Only | Somatic (Tumor-Normal) |
| :------------------------------------------------------------------------------------------------------ | :-: | :-: | :----: | :-------: | :--------: | :--------------------: |
| [DeepVariant](https://github.com/google/deepvariant)                                                    |  x  |  x  |   x    |     x     |     -      |           -            |
| [FreeBayes](https://github.com/ekg/freebayes)                                                           |  x  |  x  |   x    |     x     |     x      |           x            |
| [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/5358864757787-HaplotypeCaller) |  x  |  x  |   x    |     x     |     -      |           -            |
| [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2)                 |  x  |  x  |   x    |     -     |     x      |           x            |
| [lofreq](https://github.com/CSB5/lofreq)                                                                |  x  |  x  |   x    |     -     |     x      |           -            |
| [mpileup](https://www.htslib.org/doc/samtools-mpileup.html)                                             |  x  |  x  |   x    |     x     |     x      |           -            |
| [Strelka](https://github.com/Illumina/strelka)                                                          |  x  |  x  |   x    |     x     |     -      |           x            |
| [Manta](https://github.com/Illumina/manta)                                                              |  x  |  x  |   x    |     x     |     x      |           x            |
| [indexcov](https://github.com/brentp/goleft/tree/master/indexcov)                                       |  x  |  -  |   -    |     x     |     -      |           x            |
| [TIDDIT](https://github.com/SciLifeLab/TIDDIT)                                                          |  x  |  x  |   x    |     x     |     x      |           x            |
| [ASCAT](https://github.com/VanLoo-lab/ascat)                                                            |  x  |  x  |   -    |     -     |     -      |           x            |
| [CNVKit](https://cnvkit.readthedocs.io/en/stable/)                                                      |  x  |  x  |   -    |     x     |     x      |           x            |
| [Control-FREEC](https://github.com/BoevaLab/FREEC)                                                      |  x  |  x  |   x    |     -     |     x      |           x            |
| [MSIsensorPro](https://github.com/xjtu-omics/msisensor-pro)                                             |  x  |  x  |   x    |     -     |     -      |           x            |

## How to run ASCAT with whole-exome sequencing data?

ASCAT runs out of the box on whole genome sequencing data using iGenomes resources. While the ASCAT implementation in sarek is capable of running with whole-exome sequencing data, the needed references are currently not provided with the igenomes.config. According to the [developers](https://github.com/VanLoo-lab/ascat/issues/97) of ASCAT, loci and allele files (one file per chromosome) can be downloaded directly from the [Battenberg repository](https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52).

Please note that:

- Row names (for GC and RT correction files) should be `${chr}_${position}` (there is no SNP/probe ID for HTS data).
- All row names in GC and RT correction files should also appear in the loci files
- Loci and allele files must contain the same set of SNPs
- ASCAT developers strongly recommend using a BED file for WES/TS data. This prevents considering SNPs covered by off-target reads that would add noise to log/BAF tracks.
- The total number of GC correction loci in a sample must be at least 10% of the number of loci with logR values. If the number of GC correction loci is too small compared to the total number of loci, ASCAT will throw an error.

From 'Reference files' https://github.com/VanLoo-lab/ascat:

> For WES and targeted sequencing, we recommend using the reference files (loci, allele and logR correction files) as part of the Battenberg package. Because they require a high-resolution input, our reference files for WGS are not suitable for WES and targeted sequencing. For WES, loci and allele files from the Battenberg package can be fed into ascat.prepareHTS. For targeted sequencing, allele files from the Battenberg package can be fed into ascat.prepareTargetedSeq, which will generate cleaned loci and allele files that can be fed into ascat.prepareHTS.

### How to generate ASCAT resources for exome or targeted sequencing

1. Fetch the GC content correction and replication timing (RT) correction files from the [Dropbox links provided by the ASCAT developers](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS) and intersect the SNP coordinates with the exome target coordinates. If the target file has 'chr' prefixes, make a copy with these removed first. Extract the GC and RT information for only the on target SNPs and zip the results.

```bash
sed -e 's/chr//' targets_with_chr.bed > targets.bed

for t in GC RT
do
  unzip ${t}_G1000_hg38.zip

  cut -f 1-3 ${t}_G1000_hg38.txt > ascat_${t}_snps_hg38.txt
  tail -n +2 ascat_${t}_snps_hg38.txt | awk '{ print $2 "\t" $3-1 "\t" $3 "\t" $1 }' > ascat_${t}_snps_hg38.bed
  bedtools intersect -a ascat_${t}_snps_hg38.bed -b targets.bed | awk '{ print $1 "_" $3 }' > ascat_${t}_snps_on_target_hg38.txt

  head -n 1 ${t}_G1000_hg38.txt > ${t}_G1000_on_target_hg38.txt
  grep -f ascat_${t}_snps_on_target_hg38.txt ${t}_G1000_hg38.txt >> ${t}_G1000_on_target_hg38.txt
  zip ${t}_G1000_on_target_hg38.zip ${t}_G1000_on_target_hg38.txt

  rm ${t}_G1000_hg38.zip
done
```

2. Download the Battenberg 1000G loci and alleles files. The steps below follow downloading from the [Battenberg repository at the Oxford University Research Archive](https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52). The files are also available via Dropbox links from the same page as the GC and RT correction files above.

```bash
wget https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/files/rt435gd52w
mv rt345gd52w battenberg.zip
tar xf battenberg.zip

unzip 1000G_loci_hg38_chr.zip
cd 1000G_loci_hg38
mkdir battenberg_alleles_on_target_hg38
mv *allele* battenberg_alleles_on_target_hg38/
mkdir battenberg_loci_on_target_hg38
mv *loci* battenberg_loci_on_target_hg38/
```

3. Copy the `targets_with_chr.bed` and `GC_G1000_on_target_hg38.txt` files into the newly created `battenberg_loci_on_target_hg38` folder before running the next set of steps. ASCAT generates a list of GC correction loci with sufficient coverage in a sample, then intersects that with the list of all loci with tumour logR values in that sample. If the intersection is <10% the size of the latter, it will fail with an error. Because the Battenberg loci/allele sets are very dense, subsetting to on-target regions is still too many loci. This script ensures that all SNPs with GC correction information are included in the loci list, plus a random sample of another 30% of all on target loci. You may need to vary this proportion depending on your set of targets. A good rule of thumb is that the size of your GC correction loci list should be about 15% the size of your total loci list. This allows for a margin of error.

### 'chr'-based versus non 'chr'-based reference

Please note that loci files provided from ASCAT developers (https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WES) are not 'chr'-based (chromosome names are '1', '2', '3', etc. and not 'chr1', 'chr2', 'chr3', etc.). If your BAMs are 'chr'-based, you will need to add 'chr'

```bash
for i in {1..22} X;
  do sed -i 's/^/chr/' G1000_loci_hg19_chr${i}.txt;
done).
```

ASCAT will internally remove 'chr' so the other files (allele, GC correction and RT correction) should not be modified and chrom_names (ascat.prepareHTS) should be c(1:22,'X').

If using ASCAT provided references:

```bash

cd .../G1000_lociAll_hg38_unzipped/G1000_lociAll_hg38

# Function to check and correct 'chr' prefix
check_and_correct_chr_prefix() {
    local file=$1
    local chr_number=$2

    # Check if file exists
    if [ ! -f "$file" ]; then
        echo "Error: File $file not found."
        exit 1
    fi

    # Check first line of the file
    first_line=$(head -n 1 "$file")

    if [[ $first_line == chr${chr_number}* ]]; then
        echo "File $file already has correct 'chr' prefix. No changes needed."
    elif [[ $first_line == chrchr${chr_number}* ]]; then
        echo "File $file has duplicate 'chr' prefix. Correcting..."
        sed -i 's/^chrchr/chr/' "$file"
    elif [[ $first_line == ${chr_number}* ]]; then
        echo "File $file is missing 'chr' prefix. Adding..."
        sed -i 's/^/chr/' "$file"
    else
        echo "Error: Unexpected format in $file. Please check manually."
        exit 1
    fi
}

# Check and correct 'chr' prefix for each loci file
for i in {1..22} X; do
    check_and_correct_chr_prefix "G1000_loci_hg38_chr${i}.txt" "${i}"
done

for i in {1..22} X
do
  # Generate BED file from the tailored loci set
  awk '{ print $1 "\t" $2-1 "\t" $2 }' G1000_loci_hg38_chr${i}.txt > chr${i}.bed

  # Extract relevant GC content data for this chromosome
  grep "^chr${i}_" GC_G1000_on_target_hg38.txt > chr${i}.txt

  # Intersect BED file with target regions to find loci on target
  bedtools intersect -a chr${i}.bed -b targets_with_chr.bed | awk '{ print $1 "_" $3 }' > chr${i}_on_target.txt

  # Calculate the number of lines needed for random sampling (30% of total)
  n=$(wc -l < chr${i}_on_target.txt)
  count=$((n * 3 / 10))

  # Get loci that are both on target and match the GC content data
  grep -xf chr${i}.txt chr${i}_on_target.txt > chr${i}.temp

  # Add random subset of on-target loci to the list
  shuf -n $count chr${i}_on_target.txt >> chr${i}.temp

  # Sort, remove duplicates, and format output
  sort -n -k2 -t '_' chr${i}.temp | uniq | awk 'BEGIN { FS="_" } ; { print $1 "\t" $2 }' > battenberg_loci_on_target_hg38_chr${i}.txt
done

# Compress the resulting loci files into a zip archive
zip battenberg_loci_on_target_hg38.zip battenberg_loci_on_target_hg38_chr*.txt

```

If using Battenberg provided references:

```bash
cd battenberg_loci_on_target_hg38/
rm *chrstring*
rm 1kg.phase3.v5a_GRCh38nounref_loci_chr23.txt
for i in {1..22} X
do
  awk '{ print $1 "\t" $2-1 "\t" $2 }' 1kg.phase3.v5a_GRCh38nounref_loci_chr${i}.txt > chr${i}.bed
  grep "^${i}_" GC_G1000_on_target_hg38.txt | awk '{ print "chr" $1 }' > chr${i}.txt
  bedtools intersect -a chr${i}.bed -b targets_with_chr.bed | awk '{ print $1 "_" $3 }' > chr${i}_on_target.txt
  n=`wc -l chr${i}_on_target.txt | awk '{ print $1 }'`
  count=$((n * 3 / 10))
  grep -xf chr${i}.txt chr${i}_on_target.txt > chr${i}.temp
  shuf -n $count chr${i}_on_target.txt >> chr${i}.temp
  sort -n -k2 -t '_' chr${i}.temp | uniq | awk 'BEGIN { FS="_" } ; { print $1 "\t" $2 }' > battenberg_loci_on_target_hg38_chr${i}.txt
done
zip battenberg_loci_on_target_hg38.zip battenberg_loci_on_target_hg38_chr*.txt
```

4. Extract the alleles for the same set of SNPs. Uses a short R script defined below.

```bash
cd ../battenberg_alleles_on_target_hg38/
rm 1kg.phase3.v5a_GRCh38nounref_allele_index_chr23.txt
for i in {1..22} X
do
  Rscript intersect_ascat_alleles.R ../battenberg_loci_on_target_hg38/battenberg_loci_on_target_hg38_chr${i}.txt \
    1kg.phase3.v5a_GRCh38nounref_allele_index_chr${i}.txt battenberg_alleles_on_target_hg38_chr${i}.txt
done
zip battenberg_alleles_on_target_hg38.zip battenberg_alleles_on_target_hg38_chr*.txt
```

Rscript `intersect_ascat_alleles.R`

```bash
#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

loci = read.table(args[1], header=F, sep="\t", stringsAsFactors=F)
alleles = read.table(args[2], header=T, sep="\t", stringsAsFactors=F)

i = intersect(loci$V2, alleles$position)

out = subset(alleles, alleles$position %in% i)
write.table(out, args[3], col.names=T, row.names=F, quote=F, sep="\t")
```

5. Move or copy all of the zip files you've created to a suitable location. Specify these in your parameters, e.g.

```json
{
  "ascat_alleles": "/path/to/battenberg_alleles_on_target_hg38.zip",
  "ascat_loci": "/path/to/battenberg_loci_on_target_hg38.zip",
  "ascat_loci_gc": "/path/to/GC_G1000_on_target_hg38.zip",
  "ascat_loci_rt": "/path/to/RT_G1000_on_target_hg38.zip"
}
```

## What are the bwa, bwa-mem2 and sentieon bwa mem parameters?

For mapping, sarek follows the parameter suggestions provided in this [paper](https://www.nature.com/articles/s41467-018-06159-4):

`-K 100000000` : for deterministic pipeline results, for more info see [here](https://github.com/CCDG/Pipeline-Standardization/issues/2)

`-Y`: force soft-clipping rather than default hard-clipping of supplementary alignments

In addition, currently the mismatch penalty for reads with tumor status in the sample sheet are mapped with a mismatch penalty of `-B 3`.

## How to manage scatter/gathering (parallelization with-in each sample)

While Nextflow ensures all samples are run in parallel, the pipeline can split input files for each sample into smaller chunks which are processes in parallel.
This speeds up analysis for individual chunks, but might occupy more storage space.

Therefore, the different scatter/gather options can be set by the user:

### Split Fastq files

By default, the input fastq files are split into smaller chunks with FASTP, mapped in parallel, and then merged and duplicate marked. This can be customized by setting the parameter `--split_fastq`.
This parameter determines how many reads are within each split. Setting it to `0` will turn of any splitting and only one mapping process is run per input fastq file.

> FastP creates as many chunks as CPUs are specified (by default 12) and subdivides them further, if the number of reads in a chunk is larger then the value specified in `--split_fastq`. Thus, the parameter `--split_fastq` is an upper bound, e.g. if 1/12th of the Fastq file exceeds the provided value another fastq file will be generated.

### Intervals for Base Quality Score Recalibration and Variantcalling

The pipeline can parallelize base quality score recalibration and variant calling across genomic chunks of roughly similar sizes.
For this, a bed file containing genomic regions of interest is used, it's the intervals file.
By default, the intervals file for WGS used is the one provided by GATK (details [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035889551-When-should-I-restrict-my-analysis-to-specific-intervals-)).
When running targeted analysis, it is recommended to use the bed file containing the targeted regions.

The amount of scatter/gathering can be customized by adjusting the parameter `--nucleotides_per_second`.

> **NB:** The _same_ intervals are processed regardless of the number of groups. The number of groups however determines over how many compute nodes the analysis is scattered on.

The default value is `200000`, increasing this value will _reduce_ the number of groups that are processed in parallel.
Generally, smaller numbers of groups (each group has more regions), the slower the processing, and less storage space is consumed.
In particular, in cloud computing setting it is often advisable to reduce the number of groups to be run in parallel to reduce data staging steps.

## How to create a panel-of-normals for Mutect2

For a detailed tutorial on how to create a panel-of-normals, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132).

## Spark related issues

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

```bash
nextflow run nf-core/sarek --genome null --igenomes_ignore --fasta <custom.fasta>
```

### Overwrite specific reference files

If you don't want to use some of the provided reference genomes, they can be overwritten by either providing a new file or setting the respective file parameter to `false`, if it should be ignored:

Example for using a custom known indels file:

```bash
nextflow run nf-core/sarek --known_indels <my_known_indels.vcf.gz> --genome GRCh38.GATK
```

Example for not using known indels, but all other provided reference file:

```bash
nextflow run nf-core/sarek --known_indels false --genome GRCh38.GATK
```

### Where do the used reference genomes originate from

For GATK.GRCh38 the links for each reference file and the corresponding processes that use them is listed below. For GATK.GRCh37 the files originate from the same sources:

| File                  | Tools                                                                                                                                                                                                                                                                                                                                                                                                                                                          | Origin                                                                                                                | Docs                                                                                 |
| :-------------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------- |
| ascat_alleles         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                          | https://www.dropbox.com/s/uouszfktzgoqfy7/G1000_alleles_hg38.zip                                                      | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci            | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                          | https://www.dropbox.com/s/80cq0qgao8l1inj/G1000_loci_hg38.zip                                                         | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci_gc         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                          | https://www.dropbox.com/s/80cq0qgao8l1inj/G1000_loci_hg38.zip                                                         | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| ascat_loci_rt         | ASCAT                                                                                                                                                                                                                                                                                                                                                                                                                                                          | https://www.dropbox.com/s/xlp99uneqh6nh6p/RT_G1000_hg38.zip                                                           | https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS                   |
| bwa                   | bwa-mem                                                                                                                                                                                                                                                                                                                                                                                                                                                        | `bwa index -p bwa/${fasta.baseName} $fasta`                                                                           |                                                                                      |
| bwamem2               | bwa-mem2                                                                                                                                                                                                                                                                                                                                                                                                                                                       | `bwa-mem2 index -p bwamem2/${fasta} $fasta`                                                                           |                                                                                      |
| dragmap               | DragMap                                                                                                                                                                                                                                                                                                                                                                                                                                                        | `dragen-os --build-hash-table true --ht-reference $fasta --output-directory dragmap`                                  |                                                                                      |
| dbsnp                 | Baserecalibrator, ControlFREEC, GenotypeGVCF, HaplotypeCaller                                                                                                                                                                                                                                                                                                                                                                                                  | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| dbsnp_tbi             | Baserecalibrator, ControlFREEC, GenotypeGVCF, HaplotypeCaller                                                                                                                                                                                                                                                                                                                                                                                                  | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| dict                  | Baserecalibrator(Spark), CNNScoreVariant, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, MarkDulpicates(Spark), MergeVCFs, Mutect2, Variantrecalibrator                                                                                                                                                                                                         | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| fasta                 | ApplyBQSR(Spark), ApplyVQSR, ASCAT, Baserecalibrator(Spark), BWA, BWAMem2, CNNScoreVariant, CNVKit, ControlFREEC, DragMap, DEEPVariant, EnsemblVEP, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, FreeBayes, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, indexcov, interval building, Manta, MarkDuplicates(Spark),MergeVCFs,MSISensorPro, Mutect2, Samtools, SnpEff, Strelka, Tiddit, Variantrecalibrator | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| fasta_fai             | ApplyBQSR(Spark), ApplyVQSR, ASCAT, Baserecalibrator(Spark), BWA, BWAMem2, CNNScoreVariant, CNVKit, ControlFREEC, DragMap, DEEPVariant, EnsemblVEP, EstimateLibraryComplexity, FilterMutectCalls, FilterVariantTranches, FreeBayes, GatherPileupSummaries,GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, indexcov, interval building, Manta, MarkDuplicates(Spark),MergeVCFs,MSISensorPro, Mutect2, Samtools, SnpEff, Strelka, Tiddit, Variantrecalibrator | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle       |
| germline_resource     | GetPileupsummaries,Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                     | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| germline_resource_tbi | GetPileupsummaries,Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                     | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| intervals             | ApplyBQSR(Spark), ASCAT, Baserecalibrator(Spark), BCFTools, CNNScoreVariants, ControlFREEC, Deepvariant, FilterVariantTranches, FreeBayes, GenotypeGVCF, GetPileupSummaries, HaplotypeCaller, Strelka, mpileup, MSISensorPro, Mutect2, VCFTools                                                                                                                                                                                                                | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_indels          | BaseRecalibrator(Spark), FilterVariantTranches                                                                                                                                                                                                                                                                                                                                                                                                                 | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_indels_tbi      | BaseRecalibrator(Spark), FilterVariantTranches                                                                                                                                                                                                                                                                                                                                                                                                                 | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_snps            | BaseRecalibrator(Spark), FilterVariantTranches, VariantRecalibrator                                                                                                                                                                                                                                                                                                                                                                                            | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |                                                                                      |
| known_snps_tbi        | BaseRecalibrator(Spark), FilterVariantTranches, VariantRecalibrator                                                                                                                                                                                                                                                                                                                                                                                            | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) |
| mappability           | ControlFREEC                                                                                                                                                                                                                                                                                                                                                                                                                                                   | http://xfer.curie.fr/get/vyIi4w8EONl/out100m2_hg38.zip                                                                | http://boevalab.inf.ethz.ch/FREEC/tutorial.html                                      |
| pon                   | Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                                        | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- |
| pon_tbi               | Mutect2                                                                                                                                                                                                                                                                                                                                                                                                                                                        | [GATKBundle](https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/) | https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON- |

## How to customise SnpEff and VEP annotation

SNPeff and VEP both require a large resource of files known as a cache.
These are folders composed of multiple gigabytes of files which need to be available for the software to properly function.
To use these, supply the parameters `--vep_cache` and/or `--snpeff_cache` with the locations to the root of the annotation cache folder for each tool.

### Specify the cache location

Params `--snpeff_cache` and `--vep_cache` are used to specify the locations to the root of the annotation cache folder.
The cache will be located within a subfolder with the path `${snpeff_species}.${snpeff_version}` for SnpEff and `${vep_species}/${vep_cache_version}_${vep_genome}` for VEP.
If this directory is missing, Sarek will raise an error.

For example this is a typical folder structure for `GRCh38` and `WBCel235`, with SNPeff cache version 105 and VEP cache version 110:

```text
/data/
├─ snpeff_cache/
│  ├─ GRCh38.105/
│  ├─ WBcel235.105/
├─ vep_cache/
│  ├─ caenorhabditis_elegans/
│  │  ├─ 110_WBCel235/
│  ├─ homo_sapiens/
│  │  ├─ 110_GRCh38/
```

For this example, the parameters `--snpeff_cache /data/snpeff_cache` and `--vep_cache /data/vep_cache` would be used.
Both SnpEff and VEP will figure out internally the path towards the specific cache version / species the annotation should be performed given the parameters specified to Sarek.

### Change cache version and species

By default all is specified in the [igenomes.config](https://github.com/nf-core/sarek/blob/master/conf/igenomes.config) file.
Explanation can be found for all params in the documentation:

- [snpeff_db](https://nf-co.re/sarek/parameters#snpeff_db)
- [vep_genome](https://nf-co.re/sarek/parameters#vep_genome)
- [vep_species](https://nf-co.re/sarek/parameters#vep_species)
- [vep_cache_version](https://nf-co.re/sarek/parameters#vep_cache_version)

With the previous example of `GRCh38`, these are the values that were used for these params:

```bash
snpeff_db         = 'GRCh38.105'
vep_cache_version = '110'
vep_genome        = 'GRCh38'
vep_species       = 'homo_sapiens'
```

### Usage recommendation with AWS iGenomes

The cache for each of these annotation tools has its own structure and is frequently updated, therefore it is kept separate from AWS iGenomes. It is not recommended to put any cache for each of this annotation tools in your local AWS iGenomes folder.

A classical organisation on a shared storage area might be:

```bash
/data/igenomes/
/data/cache/snpeff_cache
/data/cache/vep_cache
```

Which can then be used this way in Sarek:

```bash
nextflow run nf-core/sarek \
    --igenomes_base /data/igenomes/ \
    --snpeff_cache /data/cache/snpeff_cache/ \
    --vep_cache /data/cache/vep_cache/ \
    ...
```

Alternatively the data may be stored on AWS S3 storage, therefore the parameters might be:

```bash
s3://my-reference-data/igenomes/
s3://my-reference-data/cache/snpeff_cache/
s3://my-reference-data/cache/vep_cache/
```

Which can then be used this way in Sarek:

```bash
nextflow run nf-core/sarek \
    --igenomes_base s3://my-reference-data/igenomes/ \
    --snpeff_cache s3://my-reference-data/cache/ensemblvep/ \
    --vep_cache s3://my-reference-data/cache/snpeff/ \
    ...
```

These params can be specified in a config file or in a profile using the params scope, or even in a json or a yaml file using the `-params-file` nextflow option.

Note: we recommend storing each annotation cache in a separate directory so each cache version is handled differently.
This may mean you have many similar directories but will dramatically reduce the storage burden on machines running the SnpEff or VEP process.

### Use annotation-cache for SnpEff and VEP

[Annotation-cache](https://annotation-cache.github.io) is an open AWS registry resource that stores a mirror of some cache files on AWS S3 which can be used with Sarek.
It contains some genome builds which can be found by checking the contents of the S3 bucket.

SNPeff and VEP cache are stored at the following location on S3:

```bash
snpeff_cache = s3://annotation-cache/snpeff_cache/
vep_cache = s3://annotation-cache/vep_cache/
```

The contents of said cache can be listed with the following command using the S3 CLI:

```bash
aws s3 --no-sign-request ls s3://annotation-cache/snpeff_cache
aws s3 --no-sign-request ls s3://annotation-cache/vep_cache/
```

Since both Snpeff and VEP are internally figuring the path towards the specific cache version / species, `annotation-cache` is using an extra set of keys to specify the species and genome build.

Which is handled internally by Sarek.

Please refer to the [annotation-cache documentation](https://annotation-cache.github.io) for more details.

### Use Sarek to download cache and annotate in one go

Both VEP and snpEff come with built-in download functionality to download the cache prior to use.
Sarek includes these as optional processes.
Use the params `--download_cache`, and specify the tool with `--tools` and Sarek will download the relevant cache (`snpeff` and/or `vep`) using their respective download functions.
It is recommended to save the cache somewhere highly accessible for subsequent runs of Sarek, so the cache does not have to be re-downloaded.

Sarek will automatically download the cache using each tools (SnpEff and/or VEP) to your work directory.
And subsequently perform the annotation of VCF files specified as an input in a samplesheet or produced by Sarek.

### Only download cache

Using the params `--build_only_index` allow for only downloading the cache for the specified tools.

### Location for the cache

Cache can be downloaded in the specified `--outdir_cache` location.
Else, it will be downloaded in `cache/` in the specified `--outdir` location.

This command could be used to download the cache for both tools in the specified `--outdir_cache` location:

```bash
nextflow run nf-core/sarek --outdir results --outdir_cache /path_to/my-own-cache --tools vep,snpeff --download_cache --build_only_index --input false
```

This command could be used to point to the recently downloaded cache and run SnpEff and VEP:

```bash
nextflow run nf-core/sarek --outdir results --vep_cache /path_to/my-own-cache/vep_cache --snpeff_cache /path_to/my-own-cache/snpeff_cache --tools vep,snpeff --input samplesheet_vcf.csv
```

Here is an example on how sarek may be used to download the SnpEff cache for Candida auris:

```bash
nextflow run nf-core/sarek --outdir results --outdir_cache /path_to/my-own-cache --tools snpeff --download_cache --build_only_index --input false --snpeff_db _candida_auris_gca_001189475 --step annotate --genome null --igenomes_ignore
```

### Create containers with pre-downloaded cache

nf-core is no longer maintaining containers with pre-downloaded cache. Hosting the cache within the container is not recommended as it can cause a number of problems. Instead we recommned using an external cache. The following is left for legacy reasons.

But for each of these tools, an helper script `build.sh` can be found at the root of the tool folder in the nf-core module repo ([snpeff](https://github.com/nf-core/modules/tree/master/modules/nf-core/snpeff) and [ensemblvep](https://github.com/nf-core/modules/tree/master/modules/nf-core/ensemblvep)), and can be adapted for your usage.

Overwritting the container declaration is then possible to accomodate for the new container.

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

### BCFTOOLS Annotate

It is possible to annotate a VCF file with a custom annotation file using [BCFTOOLS Annotate](https://samtools.github.io/bcftools/bcftools.html#annotate). This can be done by setting adding bcfann to the tools list and setting the following parameters:

- annotations: path to vcf annotation file
- annotations_index: path to vcf annotation index file
- header_lines: path to header lines file

## MultiQC related issues

### Plots for SnpEff are missing

When plots are missing, it is possible that the fasta and the custom SnpEff database are not matching https://pcingola.github.io/SnpEff/se_faq/#error_chromosome_not_found-details.
The SnpEff completes without throwing an error causing nextflow to complete successfully. An indication for the error are these lines in the `.command` files:

```text
ERRORS: Some errors were detected
Error type      Number of errors
ERROR_CHROMOSOME_NOT_FOUND      17522411
```

## Sentieon

[Sentieon](https://www.sentieon.com/) is a commercial solution to process genomics data with high computing efficiency, fast turnaround time, exceptional high accuracy, and 100% consistency.

In particular, Sentieon contains what may be view as speedup version of some standard GATK tools, like bwamem and haplotyper. Sarek contains support for some of the functions in Sentieon. In order to use those functions, the user will need to supply Sarek with a license for Sentieon.

### Setup of Sentieon license

Sentieon supply license in the form of a string-value (a url) or a file. It should be base64-encoded and stored in a nextflow secret named `SENTIEON_LICENSE_BASE64`. If a license string (url) is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 $(echo -n <sentieon_license_string> | base64 -w 0)
```

:::note
<sentieon_license_string> is formatted as `IP:Port` for example: `12.12.12.12:8990`
:::

If a license file is supplied, then the nextflow secret should be set like this:

```bash
nextflow secrets set SENTIEON_LICENSE_BASE64 \$(cat <sentieon_license_file.lic> | base64 -w 0)
```

:::note
If you're looking for documentation on how the nf-core Sentieon GitHub Actions and Sentieon License Server are set up: [Here be dragons.](https://github.com/nf-core/ops/blob/main/pulumi/sentieon_license_server/README.md)
:::

### Available Sentieon functions

Sarek contains the following Sentieon functions from [DnaSeq](https://support.sentieon.com/manual/DNAseq_usage/dnaseq/) : [bwa mem](https://support.sentieon.com/manual/usages/general/#bwa-mem-syntax), [LocusCollector](https://support.sentieon.com/manual/usages/general/#locuscollector-algorithm) + [Dedup](https://support.sentieon.com/manual/usages/general/#dedup-algorithm), [Haplotyper](https://support.sentieon.com/manual/usages/general/#haplotyper-algorithm), [GVCFtyper](https://support.sentieon.com/manual/usages/general/#gvcftyper-algorithm) and [VarCal](https://support.sentieon.com/manual/usages/general/#varcal-algorithm) + [ApplyVarCal](https://support.sentieon.com/manual/usages/general/#applyvarcal-algorithm), so the basic processing of alignment of fastq-files to VCF-files can be done using speedup Sentieon functions.

Sarek also contains the Sentieon functions [DnaScope](https://support.sentieon.com/manual/usages/general/?highlight=dnamodelapply#dnascope-algorithm) and [DNAModelApply](https://support.sentieon.com/manual/usages/general/?highlight=dnamodelapply#dnamodelapply-algorithm).

### Basic usage of Sentieon functions

To use Sentieon's aligner `bwa mem`, set the aligner option `sentieon-bwamem`.
(This can, for example, be done by adding `--aligner sentieon-bwamem` to the `nextflow run` command.)

To use Sentieon's function `Dedup`, specify `sentieon_dedup` as one of the tools.
(This can, for example, be done by adding `--tools sentieon_dedup` to the `nextflow run` command.)

To use Sentieon's function `DNAscope`, specify `sentieon_dnascope` as one of the tools.
This can, for example, be done by adding `--tools sentieon_dnascope` to the `nextflow run` command.
In order to skip Sentieon's variant-filter `DNAModelApply`, one may add `--skip_tools dnascope_filter` to the `nextflow run` command.
Sarek also provides the option `sentieon_dnascope_emit_mode` which can be used to set the [emit-mode](https://support.sentieon.com/manual/usages/general/#dnascope-algorithm) of Sentieon's dnascope.
Sentieon's dnascope can output both a vcf-file and a gvcf-file in the same run; this is achieved by setting `sentieon_dnascope_emit_mode` to `<vcf_emit_mode>,gvcf`, where `<vcf_emit_mode>` is `variant`, `confident` or `all`.

Sentieon's function `Haplotyper` is used in much the same way as `DNAscope`.
To use Sentieon's function `Haplotyper`, specify `sentieon_haplotyper` as one of the tools.
This can, for example, be done by adding `--tools sentieon_haplotyper` to the `nextflow run` command.
In order to skip the GATK-based variant-filter, one may add `--skip_tools haplotyper_filter` to the `nextflow run` command.
Sarek also provides the option `sentieon_haplotyper_emit_mode` which can be used to set the [emit-mode](https://support.sentieon.com/manual/usages/general/#haplotyper-algorithm) of Sentieon's haplotyper.
Sentieon's haplotyper can output both a vcf-file and a gvcf-file in the same run; this is achieved by setting `sentieon_haplotyper_emit_mode` to `<vcf_emit_mode>,gvcf`, where `<vcf_emit_mode>` is `variant`, `confident` or `all`.

To use Sentieon's function `GVCFtyper` along with Sention's version of VQSR (`VarCal` and `ApplyVarCal`) for joint-germline genotyping, specify `sentieon_haplotyper` as one of the tools, set the option `sentieon_haplotyper_emit_mode` to `gvcf`, and add the option `joint_germline`.
This can, for example, be done by adding `--tools sentieon_haplotyper --joint_germline --sentieon_haplotyper_emit_mode gvcf` to the `nextflow run` command.
If `sentieon_dnascope` is chosen instead of `sentieon_haplotyper`, then Sention's version of VQSR is skipped, as recommended by Sentieon.

### Joint germline variant calling

Sentieon's [GVCFtyper](https://support.sentieon.com/manual/usages/general/#gvcftyper-algorithm) does not support the [GenomicsDB](https://gatk.broadinstitute.org/hc/en-us/articles/5358869876891-GenomicsDBImport) datastore format. This means that, in contrast to the GATK based joint germline variant calling subworkflow in Sarek, the Sentieon/DNAseq based joint germline variant calling subworkflow does not use the GenomicsDB datastore format.

### QualCal (BQSR)

Currently, Sentieon's version of BQSR, QualCal, is not available in Sarek. Recent Illumina sequencers tend to provide well-calibrated BQs, so BQSR may not provide much benefit. By default Sarek runs GATK's BQSR; that can be skipped by adding the option `--skip_tools baserecalibrator`.

## Requested resources for the tools

Resource requests are difficult to generalize and are often dependent on input data size. Currently, the number of cpus and memory requested by default were adapted from tests on 5 ICGC paired whole-genome sequencing samples with approximately 40X and 80X depth.
For targeted data analysis, this is overshooting by a lot. In this case resources for each process can be limited by tailoring the request by process name as described [here](#resource-requests). If you are using sarek for a certain data type regulary, and would like to make these requests available to others on your system, an institution-specific, pipeline-specific config file can be added [here](https://github.com/nf-core/configs/tree/master/conf/pipeline/sarek).

## CNV calling with CNVkit

The CNV calling in Sarek implements the approach proposed by [CNVkit](https://cnvkit.readthedocs.io/en/stable/).
It is possible to call CNVs with whole-genome or targeted capture data (exome and amplicons): depending on the sequencing approach, Sarek applies different [settings](https://cnvkit.readthedocs.io/en/stable/nonhybrid.html) as recommended by CNVkit.

### Reference background

Given the nature of this type of CNV calling algorithms, which rely on the detection of variations in the coverage profile, the definition of a background reference in control data is known to improve the calling in targeted and hybrid capture applications. This is to ensure an accurate profiling, especially in the off-target regions.
We recommend creating a background reference with the nf-core pipeline [createpanelrefs](https://nf-co.re/createpanelrefs).

:warning: In creating a coverage reference, one should pay particular attention to:

- the control samples should be processed with the same targeted capture and sequencing technology
- if BAM files are used to compute the background, they should have been processed with the same pipeline used to call the CNVs

### Germline calling

Sarek implements the [recommended germline settings](https://cnvkit.readthedocs.io/en/stable/germline.html), i.e. applying the `--filter ci` option in the CVNkit call step.
However, this is defined at a config level by adding this option to the `ext.args`: the user can therefore choose any desired different approach by changing the arguments in a custom config.

### Somatic calling

The [available options](https://cnvkit.readthedocs.io/en/stable/tumor.html) a user can choose from for tumour analysis depend very much on the specific design being analysed. Sarek therefore doesn't implement any of these choices, i.e. it runs the CNVkit call step with default settings.
We encourage the user to verify whether particular settings might be more appropriate for their data.
