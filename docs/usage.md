# nf-core/sarek: Usage <!-- omit in toc -->

- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Main arguments](#main-arguments)
  - [-profile](#-profile)
  - [--input](#--input)
  - [--split_fastq](#--split_fastq)
  - [--trim_fastq](#--trim_fastq)
  - [--clip_r1](#--clip_r1)
  - [--clip_r2](#--clip_r2)
  - [--three_prime_clip_r1](--three_prime_clip_r1)
  - [--three_prime_clip_r2](--three_prime_clip_r2)
  - [--trim_nextseq](--trim_nextseq)
  - [--save_trimmed](--save_trimmed)
  - [--sample](#--sample)
  - [--sampleDir](#--sampledir)
  - [--annotateVCF](#--annotatevcf)
  - [--no_gvcf](#--no_gvcf)
  - [--noGVCF](#--nogvcf)
  - [--skip_qc](#--skip_qc)
  - [--skipQC](#--skipqc)
  - [--noReports](#--noreports)
  - [--nucleotides_per_second](#--nucleotides_per_second)
  - [--nucleotidesPerSecond](#--nucleotidespersecond)
  - [--step](#--step)
  - [--tools](#--tools)
  - [--sentieon](#--sentieon)
  - [--no_strelka_bp](#--no_strelka_bp)
  - [--noStrelkaBP](#--nostrelkabp)
  - [--no_intervals](#--no_intervals)
  - [--target_bed](#--target_bed)
  - [--targetBED](#--targetbed)
- [Reference genomes](#reference-genomes)
  - [--genome (using iGenomes)](#--genome-using-igenomes)
  - [--ac_loci](#--ac_loci)
  - [--acLoci](#--acloci)
  - [--ac_loci_gc](#--ac_loci_gc)
  - [--acLociGC](#--aclocigc)
  - [--bwa](#--bwa)
  - [--bwaIndex](#--bwaindex)
  - [--chr_dir](#--chr_dir)
  - [--chrDir](#--chrdir)
  - [--chr_length](#--chr_length)
  - [--chrLength](#--chrlength)
  - [--dbsnp](#--dbsnp)
  - [--dbsnp_index](#--dbsnp_index)
  - [--dbsnpIndex](#--dbsnpindex)
  - [--dict](#--dict)
  - [--fasta](#--fasta)
  - [--fasta_fai](#--fasta_fai)
  - [--fastaFai](#--fastafai)
  - [--genomeDict](#--genomedict)
  - [--genomeFile](#--genomefile)
  - [--genomeIndex](#--genomeindex)
  - [--germline_resource](#--germline_resource)
  - [--germlineResource](#--germlineresource)
  - [--germline_resource_index](#--germline_resource_index)
  - [--germlineResourceIndex](#--germlineresourceindex)
  - [--intervals](#--intervals)
  - [--known_indels](#--known_indels)
  - [--knownIndels](#--knownindels)
  - [--known_indels_index](#--known_indels_index)
  - [--knownIndelsIndex](#--knownindelsindex)
  - [--pon](#--pon)
  - [--pon_index](#--pon_index)
  - [--snpeff_db](#--snpeff_db)
  - [--snpeffDb](#--snpeffdb)
  - [--species](#--species)
  - [--vep_cache_version](#--vep_cache_version)
  - [--vepCacheVersion](#--vepcacheversion)
  - [--igenomes_ignore](#--igenomes_ignore)
  - [--igenomesIgnore](#--igenomesignore)
- [Job resources](#job-resources)
  - [Automatic resubmission](#automatic-resubmission)
  - [Custom resource requests](#custom-resource-requests)
- [AWS Batch specific parameters](#aws-batch-specific-parameters)
  - [--awsqueue](#--awsqueue)
  - [--awsregion](#--awsregion)
- [Other command line parameters](#other-command-line-parameters)
  - [--outdir](#--outdir)
- [--publish_dir_mode](#--publish_dir_mode)
- [--publishDirMode](#--publishdirmode)
  - [--sequencing_center](#--sequencing_center)
  - [--email](#--email)
  - [-name](#-name)
  - [-resume](#-resume)
  - [-c](#-c)
  - [--custom_config_version](#--custom_config_version)
  - [--custom_config_base](#--custom_config_base)
  - [--max_memory](#--max_memory)
  - [--single_cpu_mem](#--single_cpu_mem)
  - [--singleCPUMem](#--singlecpumem)
  - [--max_time](#--max_time)
  - [--max_cpus](#--max_cpus)
  - [--plaintext_email](#--plaintext_email)
  - [--monochrome_logs](#--monochrome_logs)
  - [--multiqc_config](#--multiqc_config)
  - [--ascat_ploidy](#--ascat_ploidy)
  - [--ascat_purity](#--ascat_purity)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs.
Thus the Nextflow process must run until the pipeline is finished.
We recommend that you put the process running in the background through `screen` / `tmux` or similar tool.
Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory.
We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/sarek --input sample.tsv -profile docker
```

This will launch the pipeline with the `docker` configuration profile.
See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

The nf-core/sarek pipeline comes with more documentation about running the pipeline, found in the `docs/` directory:

- [Output and how to interpret the results](output.md)
- [Extra Documentation on annotation](annotation.md)

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version.
When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since.
To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/sarek
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data.
This ensures that a specific version of the pipeline code and software are used when you run your pipeline.
If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `2.5.0`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.5.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### -profile

Use this parameter to choose a configuration profile.
Profiles can give configuration presets for different compute environments.
Note that multiple profiles can be loaded, for example: `-profile docker` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

- `awsbatch`
  - A generic configuration profile to be used with AWS Batch.
- `conda`
  - A generic configuration profile to be used with [conda](https://conda.io/docs/)
  - Pulls most software from [Bioconda](https://bioconda.github.io/)
- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from dockerhub: [`nfcore/sarek`](http://hub.docker.com/r/nfcore/sarek/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  - Pulls software from DockerHub: [`nfcore/sarek`](http://hub.docker.com/r/nfcore/sarek/)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### --input

Use this to specify the location of your input TSV file, on `mapping`, `recalibrate` and `variantcalling` steps.
For example:

```bash
--input sample.tsv
```

Multiple TSV files can be specified if the path must be enclosed in quotes

Use this to specify the location to a directory on `mapping` step with a single germline sample only.
For example:

```bash
--input PathToDirectory
```

Use this to specify the location of your VCF input file on `annotate` step.
For example:

```bash
--input sample.vcf
```

Multiple VCF files can be specified if the path must be enclosed in quotes

### --split_fastq

Use the Nextflow [`splitFastq`](https://www.nextflow.io/docs/latest/operator.html#splitfastq) operator to specify how many reads should be contained in the split fastq file.
For example:

```bash
--split_fastq 10000
```

### --sample

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--input`](#--input)

### --sampleDir

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--input`](#--input)

### --annotateVCF

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--input`](#--input)

Multiple VCF files can be specified if the path must be enclosed in quotes

### --no_gvcf

Use this to disable g.vcf from `HaplotypeCaller`.

### --noGVCF

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--no_gvcf`](#--no_gvcf)

### --skip_qc

Use this to disable specific QC and Reporting tools.
Available: `all`, `bamQC`, `BCFtools`, `FastQC`, `MultiQC`, `samtools`, `vcftools`, `versions`
Default: `None`

### --skipQC

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--skip_qc`](#--skip_qc)

### --noReports

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--skipQC`](#--skipQC)

### --nucleotides_per_second

Use this to estimate of how many seconds it will take to call variants on any interval, the default value is `1000` is it's not specified in the `<intervals>.bed` file.

### --nucleotidesPerSecond

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--nucleotides_per_second`](#--nucleotides_per_second)

### --step

Use this to specify the starting step:
Default `mapping`
Available: `mapping`, `recalibrate`, `variantcalling` and `annotate`

### --tools

Use this to specify the tools to run:
Available: `ASCAT`, `ControlFREEC`, `FreeBayes`, `HaplotypeCaller`, `Manta`, `mpileup`, `Mutect2`, `Strelka`, `TIDDIT`

### --sentieon

If [Sentieon](https://www.sentieon.com/) is available, use this to enable it for preprocessing, and variant calling.
Adds the following tools for the [`--tools`](#--tools) options: `DNAseq`, `DNAscope` and `TNscope`.

Please refer to the [nf-core/configs](https://github.com/nf-core/configs#adding-a-new-pipeline-specific-config) repository on how to make a pipeline-specific configuration file based on the [munin-sarek specific configuration file](https://github.com/nf-core/configs/blob/master/conf/pipeline/sarek/munin.config).

Or ask us on the [nf-core Slack](http://nf-co.re/join/slack) on the following channels: [#sarek](https://nfcore.slack.com/channels/sarek) [#configs](https://nfcore.slack.com/channels/configs).

### --no_strelka_bp

Use this not to use `Manta` `candidateSmallIndels` for `Strelka` as Best Practice.

### --noStrelkaBP

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--no_strelka_bp`](#--no_strelka_bp)

### --no_intervals

Disable usage of intervals file, and disable automatic generation of intervals file when none are provided.

### --target_bed

Use this to specify the target BED file for targeted or whole exome sequencing.

### --targetBED

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--target_bed`](#--target_bed)

## Reference genomes

The pipeline config files come bundled with paths to the Illumina iGenomes reference index files.
If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

### --genome (using iGenomes)

There are 2 different species supported by Sarek in the iGenomes references.
To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config).
Genomes that are supported are:

- Homo sapiens
  - `--genome GRCh37` (GATK Bundle)
  - `--genome GRCh38` (GATK Bundle)

- Mus musculus
  - `--genome GRCm38` (Ensembl)

Limited support for:

- Arabidopsis thaliana
  - `--genome TAIR10` (Ensembl)

- Bacillus subtilis 168
  - `--genome EB2` (Ensembl)

- Bos taurus
  - `--genome UMD3.1` (Ensembl)
  - `--genome bosTau8` (UCSC)

- Caenorhabditis elegans
  - `--genome WBcel235` (Ensembl)
  - `--genome ce10` (UCSC)

- Canis familiaris
  - `--genome CanFam3.1`  (Ensembl)
  - `--genome canFam3`  (UCSC)

- Danio rerio
  - `--genome GRCz10`  (Ensembl)
  - `--genome danRer10`  (UCSC)

- Drosophila melanogaster
  - `--genome BDGP6`  (Ensembl)
  - `--genome dm6`  (UCSC)

- Equus caballus
  - `--genome EquCab2`  (Ensembl)
  - `--genome equCab2`  (UCSC)

- Escherichia coli K 12 DH10B
  - `--genome EB1`  (Ensembl)

- Gallus gallus
  - `--genome Galgal4`  (Ensembl)
  - `--genome galgal4`  (UCSC)

- Glycine max
  - `--genome Gm01`  (Ensembl)

- Homo sapiens
  - `--genome hg19`  (UCSC)
  - `--genome hg38`  (UCSC)

- Macaca mulatta
  - `--genome Mmul_1`  (Ensembl)

- Mus musculus
  - `--genome mm10`  (Ensembl)

- Oryza sativa japonica
  - `--genome IRGSP-1.0`  (Ensembl)

- Pan troglodytes
  - `--genome CHIMP2.1.4`  (Ensembl)
  - `--genome panTro4`  (UCSC)

- Rattus norvegicus
  - `--genome Rnor_6.0`  (Ensembl)
  - `--genome rn6`  (UCSC)

- Saccharomyces cerevisiae
  - `--genome R64-1-1`  (Ensembl)
  - `--genome sacCer3`  (UCSC)

- Schizosaccharomyces pombe
  - `--genome EF2`  (Ensembl)

- Sorghum bicolor
  - `--genome Sbi1`  (Ensembl)

- Sus scrofa
  - `--genome Sscrofa10.2`  (Ensembl)
  - `--genome susScr3`  (UCSC)

- Zea mays
  - `--genome AGPv3`  (Ensembl)

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource.
See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh38' {
      ac_loci                  = '<path to the acLoci file>'
      ac_lociGC                = '<path to the acLociGC file>'
      bwa                      = '<path to the bwa indexes>'
      chr_dir                  = '<path to the chromosomes folder>'
      chr_length               = '<path to the chromosomes lenght file>'
      dbsnp                    = '<path to the dbsnp file>'
      dbsnp_index              = '<path to the dbsnp index>'
      dict                     = '<path to the dict file>'
      fasta                    = '<path to the fasta file>'
      fasta_fai                = '<path to the fasta index>'
      germline_resource        = '<path to the germlineResource file>'
      germline_resource_index  = '<path to the germlineResource index>'
      intervals                = '<path to the intervals file>'
      known_indels             = '<path to the knownIndels file>'
      known_indels_index       = '<path to the knownIndels index>'
      snpeff_db                = '<version of the snpEff DB>'
      species                  = '<species>'
      vep_cache_version        = '<version of the VEP cache>'
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### --ac_loci

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--ac_loci '[path to the acLoci file]'
```

### --acLoci

> :warning: This params has been removed.
> Please check: [`--ac_loci`](#--ac_loci)

### --ac_loci_gc

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--ac_loci_gc '[path to the acLociGC file]'
```

### --acLociGC

> :warning: This params has been removed.
> Please check: [`--ac_loci_gc`](#--ac_loci_gc)

### --bwa

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--bwa '[path to the bwa indexes]'
```

### --ascat_ploidy

Use this parameter to overwrite default behavior from ASCAT regarding ploidy.

### --ascat_purity

Use this parameter to overwrite default behavior from ASCAT regarding purity.

### --bwaIndex

> :warning: This params has been removed.
> Please check: [`--bwa`](#--bwa)

### --chr_dir

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--chr_dir '[path to the Chromosomes folder]'
```

### --chrDir

> :warning: This params has been removed.
> Please check: [`--chr_dir`](#--chr_dir)

### --chr_length

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--chr_length '[path to the Chromosomes length file]'
```

### --chrLength

> :warning: This params has been removed.
> Please check: [`--chr_length`](#--chr_length)

### --dbsnp

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--dbsnp '[path to the dbsnp file]'
```

### --dbsnp_index

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--dbsnp_index '[path to the dbsnp index]'
```

### --dbsnpIndex

> :warning: This params has been removed.
> Please check: [`--dbsnp_index`](#--dbsnp_index)

### --dict

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--dict '[path to the dict file]'
```

### --fasta

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta '[path to the reference fasta file]'
```

### --fasta_fai

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta_fai '[path to the reference index]'
```

### --fastaFai

> :warning: This params has been removed.
> Please check: [`--dbsnp_index`](#--dbsnp_index)

### --genomeDict

> :warning: This params has been removed.
> Please check: [`--dict`](#--dict)

### --genomeFile

> :warning: This params has been removed.
> Please check: [`--fasta`](#--fasta)

### --genomeIndex

> :warning: This params has been removed.
> Please check: [`--fastaFai`](#--fastaFai)

### --germline_resource

The [germline resource VCF file](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php#--germline-resource) (bgzipped and tabixed) needed by GATK4 Mutect2 is a collection of calls that are likely present in the sample, with allele frequencies.
The AF info field must be present.
You can find a smaller, stripped gnomAD VCF file (most of the annotation is removed and only calls signed by PASS are stored) in the iGenomes Annotation/GermlineResource folder.
To add your own germline resource supply

```bash
--germline_resource '[path to my resource.vcf.gz]'
```

### --germlineResource

> :warning: This params has been removed.
> Please check: [`--germline_resource`](#--germline_resource)

### --germline_resource_index

Tabix index of the germline resource specified at [`--germline_resource`](#--germline_resource).
To add your own germline resource supply

```bash
--germline_resource_index '[path to my resource.vcf.gz.idx]'
```

### --germlineResourceIndex

> :warning: This params has been removed.
> Please check: [`--germline_resource_index`](#--germline_resource_index)

### --intervals

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--intervals '[path to the intervals file]'
```

### --known_indels

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--known_indels '[path to the knownIndels file]'
```

### --knownIndels

> :warning: This params has been removed.
> Please check: [`--known_indels`](#--known_indels)

### --known_indels_index

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--known_indels_index '[path to the knownIndels index]'
```

### --knownIndelsIndex

> :warning: This params has been removed.
> Please check: [`--known_indels_index`](#--known_indels_index)

### --pon

When a panel of normals [PON](https://gatkforums.broadinstitute.org/gatk/discussion/24057/how-to-call-somatic-mutations-using-gatk4-mutect2#latest) is defined, you will get filtered somatic calls as a result.
Without PON, there will be no calls with PASS in the INFO field, only an _unfiltered_ VCF is written.
It is recommended to make your own panel-of-normals, as it depends on sequencer and library preparation.
For tests in iGenomes there is a dummy PON file in the Annotation/GermlineResource directory, but it _should not be used_ as a real panel-of-normals file.
Provide your PON by:

```bash
--pon '[path to the PON VCF]'
```

If the PON file is bgzipped, there has to be a tabixed index file at the same directory.

### --pon_index

Tabix index of the panel-of-normals bgzipped VCF file.

### --snpeff_db

If you prefer, you can specify the DB version when you run the pipeline:

```bash
--snpeff_db '[version of the snpEff DB]'
```

### --snpeffDb

> :warning: This params has been removed.
> Please check: [`--snpeff_db`](#--snpeff_db)

### --species

This specifies the species used for running VEP annotation. For human data, this needs to be set to `homo_sapiens`, for mouse data `mus_musculus` as the annotation needs to know where to look for appropriate annotation references. If you use iGenomes or a local resource with `genomes.conf`, this has already been set for you appropriately.

### --vep_cache_version

If you prefer, you can specify the cache version when you run the pipeline:

```bash
--vep_cache_version '[version of the VEP cache]'
```

### --vepCacheVersion

> :warning: This params has been removed.
> Please check: [`--vep_cache_version`](#--vep_cache_version)

### --igenomes_ignore

Do not load `igenomes.config` when running the pipeline.
You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.

```bash
--igenomes_ignore
```

### --igenomesIgnore

> :warning: This params has been removed.
> Please check: [`--igenomes_ignore`](#--igenomes_ignore)

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time.
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original).
If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file.
See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository.
Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below).
You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-core-invite.herokuapp.com/).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration.
Please use the `-awsbatch` profile and then specify all of the following parameters.

### --awsqueue

The JobQueue that you intend to use on AWS Batch.

### --awsregion

The AWS region to run your job in.
Default is set to `eu-west-1` but can be adjusted to your needs.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### --outdir

The output directory where the results will be saved.
Default: `results/

## --publish_dir_mode

The file publishing method.
Available: `symlink`, `rellink`, `link`, `copy`, `copyNoFollow`, `move`
Default: `copy`

## --publishDirMode

> :warning: This params is deprecated -- it will be removed in a future release.
> Please check: [`--publish_dir_mode`](#--publish_dir_mode)

### --sequencing_center

The sequencing center that will be used in the BAM CN field

### --email

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits.
If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### -name

Name for the pipeline run.
If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### -resume

Specify this when restarting a pipeline.
Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`.
Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### -c

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### --custom_config_version

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`.
This was implemented for reproducibility purposes.
Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### --custom_config_base

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet.
If you don't need them, then this is not a problem.
If you do need them, you should download the files from the repo and tell nextflow where to find them with the `custom_config_base` option.
For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### --max_memory

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit eg. `--max_memory '8.GB'`

### --single_cpu_mem

Use to set memory for a single CPU.
Should be a string in the format integer-unit eg. `--single_cpu_mem '8.GB'`

### --singleCPUMem

> :warning: This params has been removed.
> Please check: [`--single_cpu_mem`](#--single_cpu_mem)

### --max_time

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit eg. `--max_time '2.h'`

### --max_cpus

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit eg. `--max_cpus 1`

### --plaintext_email

Set to receive plain-text e-mails instead of HTML formatted.

### --monochrome_logs

Set to disable colourful command line output and live life in monochrome.

### --multiqc_config

Specify a path to a custom MultiQC configuration file.
