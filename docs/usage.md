# nf-core/sarek: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/sarek/usage](https://nf-co.re/sarek/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

Sarek is a workflow designed to detect variants on whole genome or targeted sequencing data.
Initially designed for Human, and Mouse, it can work on any species with a reference genome.
Sarek can also handle tumour / normal pairs and could include additional relapses.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/sarek --input <sample.tsv> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/sarek
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `2.6.1`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.6.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/sarek`](http://hub.docker.com/r/nfcore/sarek/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/sarek`](http://hub.docker.com/r/nfcore/sarek/)
* `podman`
  * A generic configuration profile to be used with [Podman](https://podman.io/)
  * Pulls software from Docker Hub: [`nfcore/sarek`](https://hub.docker.com/r/nfcore/sarek/)
* `shifter`
  * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
  * Pulls software from Docker Hub: [`nfcore/sarek`](https://hub.docker.com/r/nfcore/sarek/)
* `charliecloud`
  * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
  * Pulls software from Docker Hub: [`nfcore/sarek`](https://hub.docker.com/r/nfcore/sarek/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters
* `test_annotation`
  * A profile with a complete configuration for automated testing
  * Input data is a `VCF` for testing annotation
* `test_use_gatk_spark`
  * A profile with a complete configuration for automated testing
  * Specify `--use_gatk_spark`
* `test_split_fastq`
  * A profile with a complete configuration for automated testing
  * Specify `--split_fastq 500`
* `test_targeted`
  * A profile with a complete configuration for automated testing
  * Include link to a target `BED` file and use `Manta` and `Strelka` for Variant Calling
* `test_tool`
  * A profile with a complete configuration for automated testing
  * Test directly Variant Calling with a specific TSV file and `--step variantcalling`
* `test_trimming`
  * A profile with a complete configuration for automated testing
  * Test trimming options
* `test_umi_qiaseq`
  * A profile with a complete configuration for automated testing
  * Test a specific `UMI` structure
* `test_umi_tso`
  * A profile with a complete configuration for automated testing
  * Test a specific `UMI` structure

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `VEP` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: VEP {
    memory = 32.GB
  }
}
```

Alternatively, to give the workflow both processes `VEP` and `VEPmerge` 32GB of memory, you could use the following config:

```nextflow
process {
  withLabel: VEP {
    memory = 32.GB
  }
}
```

To find the exact name of a process you wish to modify the compute resources, check the live-status of a nextflow run displayed on your terminal or check the nextflow error for a line like so: `Error executing process > 'bwa'`. In this case the name to specify in the custom config file is `bwa`.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
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

`sarek`, our main container is designed using [Conda](https://conda.io/).

[![sarek-docker status](https://img.shields.io/docker/automated/nfcore/sarek.svg)](https://hub.docker.com/r/nfcore/sarek)

Based on [nfcore/base:1.12.1](https://hub.docker.com/r/nfcore/base/tags), it contains:

* **[ASCAT](https://github.com/Crick-CancerGenomics/ascat)** 2.5.2
* **[AlleleCount](https://github.com/cancerit/alleleCount)** 4.0.2
* **[BCFTools](https://github.com/samtools/bcftools)** 1.9
* **[bwa](https://github.com/lh3/bwa)** 0.7.17
* **[bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)** 2.0
* **[CNVkit](https://github.com/etal/cnvkit)** 0.9.6
* **[Control-FREEC](https://github.com/BoevaLab/FREEC)** 11.6
* **[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** 0.11.9
* **[fgbio](https://github.com/fulcrumgenomics/fgbio)** 1.1.0
* **[FreeBayes](https://github.com/ekg/freebayes)** 1.3.2
* **[GATK4-spark](https://github.com/broadinstitute/gatk)** 4.1.7.0
* **[GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)** 1.0
* **[ggplot2](https://github.com/tidyverse/ggplot2)** 3.3.0
* **[HTSlib](https://github.com/samtools/htslib)** 1.9
* **[Manta](https://github.com/Illumina/manta)** 1.6.0
* **[msisensor](https://github.com/ding-lab/msisensor)** 0.5
* **[MultiQC](https://github.com/ewels/MultiQC/)** 1.8
* **[Qualimap](http://qualimap.bioinfo.cipf.es)** 2.2.2d
* **[SAMBLASTER](https://github.com/GregoryFaust/samblaster)** 0.1.24
* **[samtools](https://github.com/samtools/samtools)** 1.9
* **[snpEff](http://snpeff.sourceforge.net/)** 4.3.1t
* **[Strelka2](https://github.com/Illumina/strelka)** 2.9.10
* **[TIDDIT](https://github.com/SciLifeLab/TIDDIT)** 2.7.1
* **[pigz](https://zlib.net/pigz/)** 2.3.4
* **[Trim Galore](https://github.com/FelixKrueger/TrimGalore)** 0.6.5
* **[VCFanno](https://github.com/brentp/vcfanno)** 0.3.2
* **[VCFtools](https://vcftools.github.io/index.html)** 0.1.16
* **[VEP](https://github.com/Ensembl/ensembl-vep)** 99.2

For annotation, the main container can be used, but then cache has to be downloaded, or additional containers are available with cache.

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
