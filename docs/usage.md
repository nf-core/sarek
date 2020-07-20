# nf-core/sarek: Usage <!-- omit in toc -->

- [Introduction](#introduction)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Main arguments](#main-arguments)
  - [-profile](#-profile)
  - [--input](#--input)
  - [--step](#--step)
  - [--help](#--help)
  - [--no_intervals](#--no_intervals)
  - [--nucleotides_per_second](#--nucleotides_per_second)
  - [--sentieon](#--sentieon)
  - [--skip_qc](#--skip_qc)
  - [--target_bed](#--target_bed)
  - [--tools](#--tools)
- [Modify fastqs (trim/split)](#modify-fastqs-trimsplit)
  - [--trim_fastq](#--trim_fastq)
  - [--clip_r1](#--clip_r1)
  - [--clip_r2](#--clip_r2)
  - [--three_prime_clip_r1](#--three_prime_clip_r1)
  - [--three_prime_clip_r2](#--three_prime_clip_r2)
  - [--trim_nextseq](#--trim_nextseq)
  - [--save_trimmed](#--save_trimmed)
  - [--split_fastq](#--split_fastq)
- [Preprocessing](#preprocessing)
  - [--markdup_java_options](#--markdup_java_options)
  - [--no_gatk_spark](#--no_gatk_spark)
  - [--save_bam_mapped](#--save_bam_mapped)
  - [--skip_markduplicates](#--skip_markduplicates)
- [Variant Calling](#variant-calling)
  - [--ascat_ploidy](#--ascat_ploidy)
  - [--ascat_purity](#--ascat_purity)
  - [--cf_coeff](#--cf_coeff)
  - [--cf_ploidy](#--cf_ploidy)
  - [--cf_window](#--cf_window)
  - [--no_gvcf](#--no_gvcf)
  - [--no_strelka_bp](#--no_strelka_bp)
  - [--umi](#--umi)
  - [--read_structure1](#--read_structure1)
  - [--read_structure2](#--read_structure2)
  - [--pon](#--pon)
  - [--pon_index](#--pon_index)
- [Annotation](#annotation)
  - [--annotate_tools](#--annotate_tools)
  - [--annotation_cache](#--annotation_cache)
  - [--snpeff_cache](#--snpeff_cache)
  - [--vep_cache](#--vep_cache)
  - [--cadd_cache](#--cadd_cache)
  - [--cadd_indels](#--cadd_indels)
  - [--cadd_indels_tbi](#--cadd_indels_tbi)
  - [--cadd_wg_snvs](#--cadd_wg_snvs)
  - [--cadd_wg_snvs_tbi](#--cadd_wg_snvs_tbi)
  - [--genesplicer](#--genesplicer)
- [Reference genomes](#reference-genomes)
  - [--genome (using iGenomes)](#--genome-using-igenomes)
  - [--igenomes_base](#--igenomes_base)
  - [--igenomes_ignore](#--igenomes_ignore)
  - [--genomes_base](#--genomes_base)
  - [--save_reference](#--save_reference)
  - [--ac_loci](#--ac_loci)
  - [--ac_loci_gc](#--ac_loci_gc)
  - [--bwa](#--bwa)
  - [--chr_dir](#--chr_dir)
  - [--chr_length](#--chr_length)
  - [--dbsnp](#--dbsnp)
  - [--dbsnp_index](#--dbsnp_index)
  - [--dict](#--dict)
  - [--fasta](#--fasta)
  - [--fasta_fai](#--fasta_fai)
  - [--germline_resource](#--germline_resource)
  - [--germline_resource_index](#--germline_resource_index)
  - [--intervals](#--intervals)
  - [--known_indels](#--known_indels)
  - [--known_indels_index](#--known_indels_index)
  - [--mappability](#--mappability)
  - [--snpeff_db](#--snpeff_db)
  - [--species](#--species)
  - [--vep_cache_version](#--vep_cache_version)
- [Other command line parameters](#other-command-line-parameters)
  - [--outdir](#--outdir)
  - [--publish_dir_mode](#--publish_dir_mode)
  - [--sequencing_center](#--sequencing_center)
  - [--multiqc_config](#--multiqc_config)
  - [--monochrome_logs](#--monochrome_logs)
  - [--email](#--email)
  - [--email_on_fail](#--email_on_fail)
  - [--plaintext_email](#--plaintext_email)
  - [--max_multiqc_email_size](#--max_multiqc_email_size)
  - [-name](#-name)
  - [-resume](#-resume)
  - [-c](#-c)
  - [--custom_config_version](#--custom_config_version)
  - [--custom_config_base](#--custom_config_base)
- [Job resources](#job-resources)
  - [Automatic resubmission](#automatic-resubmission)
  - [Custom resource requests](#custom-resource-requests)
  - [--max_memory](#--max_memory)
  - [--max_time](#--max_time)
  - [--max_cpus](#--max_cpus)
  - [--single_cpu_mem](#--single_cpu_mem)
- [AWSBatch specific parameters](#awsbatch-specific-parameters)
  - [--awsqueue](#--awsqueue)
  - [--awsregion](#--awsregion)
  - [--awscli](#--awscli)
- [Deprecated params](#deprecated-params)
  - [--annotateVCF](#--annotatevcf)
  - [--noGVCF](#--nogvcf)
  - [--noReports](#--noreports)
  - [--noStrelkaBP](#--nostrelkabp)
  - [--nucleotidesPerSecond](#--nucleotidespersecond)
  - [--publishDirMode](#--publishdirmode)
  - [--sample](#--sample)
  - [--sampleDir](#--sampledir)
  - [--skipQC](#--skipqc)
  - [--targetBED](#--targetbed)

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
nextflow run nf-core/sarek --input <sample.tsv> -profile docker
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

First, go to the [nf-core/sarek releases page](https://github.com/nf-core/sarek/releases) and find the latest version number - numeric only (eg. `2.6`).
Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.6`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### -profile

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](http://docker.com/)
  - Pulls software from dockerhub: [`nfcore/sarek`](http://hub.docker.com/r/nfcore/sarek/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  - Pulls software from DockerHub: [`nfcore/sarek`](http://hub.docker.com/r/nfcore/sarek/)
- `conda`
  - Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  - A generic configuration profile to be used with [conda](https://conda.io/docs/)
  - Pulls most software from [Bioconda](https://bioconda.github.io/)
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_annotation`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_no_gatk_spark`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_split_fastq`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_targeted`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_tool`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_trimming`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### --input

Use this to specify the location of your input TSV file.
For example:
TSV file should correspond to the correct step, see [`--step`](#--step) and [input](input.md) documentation for more information

```bash
--input <sample.tsv>
```

Multiple TSV files can be specified, using a [glob path](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob), if enclosed in quotes.

Use this to specify the location to a directory with fastq files for the `mapping` step of single germline samples only.
For example:

```bash
--input </path/to/directory>
```

Use this to specify the location of your VCF input file on `annotate` step.
For example:

```bash
--input <sample.vcf.gz>
```

Multiple VCF files can be specified, using a [glob path](https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob), if enclosed in quotes.

### --step

Use this to specify the starting step:
Default `mapping`
Available: `mapping`, `prepare_recalibration`, `recalibrate`, `variant_calling`, `annotate`, `Control-FREEC`

### --help

Will display the help message

### --no_intervals

Disable usage of intervals file, and disable automatic generation of intervals file when none are provided.

### --nucleotides_per_second

Use this to estimate of how many seconds it will take to call variants on any interval, the default value is `1000` is it's not specified in the [`intervals`](#--intervals) file.

### --sentieon

If [Sentieon](https://www.sentieon.com/) is available, use this to enable it for preprocessing, and variant calling.
Adds the following tools for the [`--tools`](#--tools) options: `DNAseq`, `DNAscope` and `TNscope`.

More information in the [sentieon](sentieon.md) documentation.

### --skip_qc

Use this to disable specific QC and Reporting tools.
Multiple tools can be specified, separated by commas.
Available: `all`, `bamQC`, `BaseRecalibrator`, `BCFtools`, `Documentation`, `FastQC`, `MultiQC`, `samtools`, `vcftools`, `versions`
Default: `None`

### --target_bed

Use this to specify the target BED file for targeted or whole exome sequencing.

### --tools

Use this parameter to specify the variant calling and annotation tools to be used.
Multiple tools can be specified, separated by commas.
For example:

```bash
--tools 'Strelka,mutect2,SnpEff'
```

Available variant callers: `ASCAT`, `ControlFREEC`, `FreeBayes`, `HaplotypeCaller`, `Manta`, `mpileup`, `MSIsensor`, `Mutect2`, `Strelka`, `TIDDIT`.

> **WARNING** Not all variant callers are available for both germline and somatic variant calling.
For more details please check the [variant calling](variant_calling.md) extra documentation.

Available annotation tools: `VEP`, `SnpEff`, `merge`. For more details, please check the [annotation](annotation.md) extra documentation.

## Modify fastqs (trim/split)

### --trim_fastq

Use this to perform adapter trimming with [Trim Galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

### --clip_r1

Instructs Trim Galore to remove a number of bp from the 5' end of read 1 (or single-end reads).
This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at the 5' end.

### --clip_r2

Instructs Trim Galore to remove a number of bp from the 5' end of read 2 (paired-end reads only).
This may be useful if the qualities were very poor, or if there is some sort of unwanted bias at the 5' end.

### --three_prime_clip_r1

Instructs Trim Galore to remove a number of bp from the 3' end of read 1 (or single-end reads) AFTER adapter/quality trimming has been performed.
This may remove some unwanted bias from the 3' end that is not directly related to adapter sequence or basecall quality.

### --three_prime_clip_r2

Instructs Trim Galore to remove a number of bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed.
This may remove some unwanted bias from the 3' end that is not directly related to adapter sequence or basecall quality.

### --trim_nextseq

This enables the option `--nextseq-trim=3'CUTOFF` within `Cutadapt`, which will set a quality cutoff (that is normally given with `-q` instead), but qualities of G bases are ignored.
This trimming is in common for the NextSeq and NovaSeq-platforms, where basecalls without any signal are called as high-quality G bases.

### --save_trimmed

Option to keep trimmed FASTQs

### --split_fastq

Use the Nextflow [`splitFastq`](https://www.nextflow.io/docs/latest/operator.html#splitfastq) operator to specify how many reads should be contained in the split fastq file.
For example:

```bash
--split_fastq 10000
```

## Preprocessing

### --markdup_java_options

To control the java options necessary for the GATK `MarkDuplicates` process, you can set this parameter.
Default: "-Xms4000m -Xmx7g"
For example:

```bash
--markdup_java_options "-Xms4000m -Xmx7g"
```

### --no_gatk_spark

Use this to disable usage of GATK Spark implementation of their tools in local mode.

### --save_bam_mapped

Will save mapped BAMs.

### --skip_markduplicates

Will skip MarkDuplicates. This params will also save the mapped BAMS, to enable restart from step `prepare_recalibration`

## Variant Calling

### --ascat_ploidy

Use this parameter to overwrite default behavior from ASCAT regarding ploidy.
Requires that [`--ascat_purity`](#--ascat_purity) is set

### --ascat_purity

Use this parameter to overwrite default behavior from ASCAT regarding purity.
Requires that [`--ascat_ploidy`](#--ascat_ploidy) is set

### --cf_coeff

Control-FREEC `coefficientOfVariation`
Default: 0.015

### --cf_ploidy

Control-FREEC `ploidy`
Default: 2

### --cf_window

Control-FREEC `window size`
Default: Disabled

### --no_gvcf

Use this to disable g.vcf output from `GATK HaplotypeCaller`.

### --no_strelka_bp

Use this not to use `Manta` `candidateSmallIndels` for `Strelka` (not recommended by Broad Institute's Best Practices).

### --pon

When a panel of normals [PON](https://gatkforums.broadinstitute.org/gatk/discussion/24057/how-to-call-somatic-mutations-using-gatk4-mutect2#latest) is defined, it will be use to filter somatic calls.
Without PON, there will be no calls with PASS in the INFO field, only an _unfiltered_ VCF is written.
It is recommended to make your own panel-of-normals, as it depends on sequencer and library preparation.
For tests in iGenomes there is a dummy PON file in the Annotation/GermlineResource directory, but it _should not be used_ as a real panel-of-normals file.
Provide your PON by:

```bash
--pon </path/to/PON.vcf.gz>
```

PON file should be bgzipped.

### --pon_index

Tabix index of the panel-of-normals bgzipped VCF file.
If none provided, will be generated automatically from the panel-of-normals bgzipped VCF file.

## Annotation

### --annotate_tools

Specify from which tools Sarek should look for VCF files to annotate, only for step `Annotate`.
Available: `HaplotypeCaller`, `Manta`, `Mutect2`, `Strelka`, `TIDDIT`
Default: `None`

### --annotation_cache

Enable usage of annotation cache, and disable usage of already built containers within Sarek.
For more information, follow the [annotation guidelines](annotation.md#using-downloaded-cache).

### --snpeff_cache

To be used conjointly with [`--annotation_cache`](#--annotation_cache), specify the cache snpEff directory:

```bash
--snpeff_cache </path/to/snpeff_cache>
```

### --vep_cache

To be used conjointly with [`--annotation_cache`](#--annotation_cache), specify the cache VEP directory:

```bash
--vep_cache </path/to/vep_cache>
```

### --cadd_cache

Enable CADD cache.

### --cadd_indels

Path to CADD InDels file.

### --cadd_indels_tbi

Path to CADD InDels index.

### --cadd_wg_snvs

Path to CADD SNVs file.

### --cadd_wg_snvs_tbi

Path to CADD SNVs index.

### --genesplicer

Enable genesplicer within VEP.

### --umi

If provided, UMIs steps will be run to extract and annotate the reads with UMIs and create consensus reads: this part of the pipeline uses *FGBIO* to convert the fastq files into a unmapped BAM, where reads are tagged with the UMIs extracted from the fastq sequences. In order to allow the correct tagging, the UMI sequence must be contained in the read sequence itself, and not in the FASTQ name.
Following this step, the uBam is aligned and reads are then grouped based on mapping position and UMI tag.
Finally, reads in the same groups are collapsed to create a consensus read. To create consensus, we have chosen to use the *adjacency method* [ref](https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/).
In order for the correct tagging to be performed, a read structure needs to be  specified as indicated below.

### --read_structure1

When processing UMIs, a read structure should always be provided for each of the fastq files, to allow the correct annotation of the bam file. If the read does not contain any UMI, the structure will be +T (i.e. only template of any length). The read structure follows a format adopted by different tools, and described [here](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)

### --read_structure2

When processing UMIs, a read structure should always be provided for each of the fastq files, to allow the correct annotation of the bam file. If the read does not contain any UMI, the structure will be +T (i.e. only template of any length). The read structure follows a format adopted by different tools, and described [here](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)

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
    '<GENOME>' {
      ac_loci                  = '</path/to/reference.loci>'
      ac_loci_gc               = '</path/to/ac_loci_gc file>'
      bwa                      = '</path/to/bwa indexes>'
      chr_dir                  = '</path/to/chromosomes/>'
      chr_length               = '</path/to/reference.len>'
      dbsnp                    = '</path/to/dbsnp.vcf.gz>'
      dbsnp_index              = '</path/to/dbsnp.vcf.gz.tbi>'
      dict                     = '</path/to/reference.dict>'
      fasta                    = '</path/to/reference.fasta>'
      fasta_fai                = '</path/to/reference.fasta.fai>'
      germline_resource        = '</path/to/germline_resource.vcf.gz>'
      germline_resource_index  = '</path/to/germline_resource.vcf.gz.tvi>'
      intervals                = '</path/to/reference.intervals>'
      known_indels             = '</path/to/known_indels.vcf.gz>'
      known_indels_index       = '</path/to/known_indels.vcf.gz.tbi>'
      mappability              = '</path/to/reference.gem>'
      snpeff_db                = '<snpEff DB>'
      species                  = '<species>'
      vep_cache_version        = '<VEP cache version'
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

### --igenomes_base

Specify base path to AWS iGenomes
Default: `s3://ngi-igenomes/igenomes/`

### --igenomes_ignore

Do not load `igenomes.config` when running the pipeline.
You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.
This option will load the `genomes.config` file instead.
You can then specify the `--genome custom` and specify any reference file on the command line or within a config file.

```bash
--igenomes_ignore
```

### --genomes_base

Specify base path to reference genome

### --save_reference

Enable saving reference indexes and other files built within Sarek.

```bash
--save_reference
```

### --ac_loci

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--ac_loci <path/to/reference.loci>
```

### --ac_loci_gc

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--ac_loci_gc <path/to/reference.loci.gc>
```

### --bwa

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

> If none provided, will be generated automatically from the fasta reference.

```bash
--bwa <path/to/BWA/indexes>
```

### --chr_dir

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--chr_dir <path/to/chromosomes/>
```

### --chr_length

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--chr_length <path/to/chromosomes.len>
```

### --dbsnp

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--dbsnp <path/to/dbsnp.vcf.gz>
```

### --dbsnp_index

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

> If none provided, will be generated automatically from the fasta reference.

```bash
--dbsnp_index <path/to/dbsnp.vcf.gz.tbi>
```

### --dict

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

> If none provided, will be generated automatically from the fasta reference.

```bash
--dict <path/to/reference.dict>
```

### --fasta

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta <path/to/reference.fasta>
```

### --fasta_fai

> If none provided, will be generated automatically from the fasta reference.

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--fasta_fai <path/to/reference.fasta.fai>
```

### --germline_resource

The [germline resource VCF file](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php#--germline-resource) (bgzipped and tabixed) needed by GATK4 Mutect2 is a collection of calls that are likely present in the sample, with allele frequencies.
The AF info field must be present.
You can find a smaller, stripped gnomAD VCF file (most of the annotation is removed and only calls signed by PASS are stored) in the iGenomes Annotation/GermlineResource folder.
If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--germline_resource </path/to/resource.vcf.gz>
```

### --germline_resource_index

Tabix index of the germline resource specified at [`--germline_resource`](#--germline_resource).
If you prefer, you can specify the full path to your reference genome when you run the pipeline:

> If none provided, will be generated automatically from the fasta reference.

```bash
--germline_resource_index </path/to/resource.vcf.gz.tbi>
```

### --intervals

Used to speed up Preprocessing and/or Variant Calling, for more information, read the [intervals section in the extra documentation on reference](reference.md#Intervals).

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

> If none provided, will be generated automatically from the fasta reference.
> Use [--no_intervals](#--no_intervals) to disable automatic generation

```bash
--intervals <path/to/reference.intervals>
```

### --known_indels

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--known_indels <path/to/known_indels.vcf.gz>
```

### --known_indels_index

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

> If none provided, will be generated automatically from the fasta reference.

```bash
--known_indels_index  <path/to/known_indels.vcf.gz.tbi>
```

### --mappability

If you prefer, you can specify the full path to your reference genome when you run the pipeline:

```bash
--mappability <path/to/reference.gem>
```

### --snpeff_db

If you prefer, you can specify the DB version when you run the pipeline:

```bash
--snpeff_db <snpEff DB>
```

### --species

This specifies the species used for running VEP annotation. For human data, this needs to be set to `homo_sapiens`, for mouse data `mus_musculus` as the annotation needs to know where to look for appropriate annotation references. If you use iGenomes or a local resource with `genomes.conf`, this has already been set for you appropriately.

### --vep_cache_version

If you prefer, you can specify the cache version when you run the pipeline:

```bash
--vep_cache_version <VEP cache version>
```

## Other command line parameters

### --outdir

The output directory where the results will be saved.
Default: `results/`

### --publish_dir_mode

The file publishing method.
Available: `symlink`, `rellink`, `link`, `copy`, `copyNoFollow`, `move`
Default: `copy`

### --sequencing_center

The sequencing center that will be used in the BAM CN field

### --multiqc_config

Specify a path to a custom MultiQC configuration file.

### --monochrome_logs

Set to disable colourful command line output and live life in monochrome.

### --email

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits.
If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### --email_on_fail

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### --plaintext_email

Set to receive plain-text e-mails instead of HTML formatted.

### --max_multiqc_email_size

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

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

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

### --max_memory

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit eg. `--max_memory '8.GB'`

### --max_time

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit eg. `--max_time '2.h'`

### --max_cpus

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit eg. `--max_cpus 1`

### --single_cpu_mem

Use to set memory for a single CPU.
Should be a string in the format integer-unit eg. `--single_cpu_mem '8.GB'`

## AWSBatch specific parameters

Running the pipeline on AWSBatch requires a couple of specific parameters to be set according to your AWSBatch configuration.
Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### --awsqueue

The JobQueue that you intend to use on AWSBatch.

### --awsregion

The AWS region to run your job in.
Default is set to `eu-west-1` but can be adjusted to your needs.

### --awscli

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI.
Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Deprecated params

> **WARNING** These params are deprecated -- They will be removed in a future release.

### --annotateVCF

> Please check: [`--input`](#--input)

### --noGVCF

> Please check: [`--no_gvcf`](#--no_gvcf)

### --noReports

> Please check: [`--skipQC`](#--skipQC)

### --noStrelkaBP

> Please check: [`--no_strelka_bp`](#--no_strelka_bp)

### --nucleotidesPerSecond

> Please check: [`--nucleotides_per_second`](#--nucleotides_per_second)

### --publishDirMode

> Please check: [`--publish_dir_mode`](#--publish_dir_mode)

### --sample

> Please check: [`--input`](#--input)

### --sampleDir

> Please check: [`--input`](#--input)

### --skipQC

> Please check: [`--skip_qc`](#--skip_qc)

### --targetBED

> Please check: [`--target_bed`](#--target_bed)
