# Try the workflow on UPPMAX clusters

This small tutorial will explain to you how to run CAW on a small sample test data on a Swedish UPPMAX cluster. Some variables are specific, but it can be easily modified to suit any clusters.

## Make a test directory

```bash
mkdir test_CAW
cd test_CAW
```

## Test the workflow on a test tiny set

This workflow itself needs no installation. Nextflow will automatically fetch it from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.

```bash
nextflow run SciLifeLab/CAW --test
```

If you're using a Swedish UPPMAX cluster, don't forget to provide your project ID.

```bash
nextflow run SciLifeLab/CAW --test --project <UPPMAX Project ID>
```

# Other possibilities for advanced users

## Clone the repository and test the workflow on a test tiny set

You can download the repository yourself from GitHub and run them directly:

```bash
git clone https://github.com/SciLifeLab/CAW test_CAW
cd test_CAW
nextflow run main.nf --test
```

## Load Nextflow

If you're running on a Swedish UPPMAX cluster you can load Nextflow as an environment module:

```bash
module load Nextflow
```

Environnement variables are set up each time the module is loaded, so you might want to set them up after loading the module.

## Running tests in interactive mode on milou

You can try the test data by changing to the interactive mode on milou and run the test tiny set like:

```
$ interactive -A <UPPMAX Project ID> -p node
[ ... login messages ... ]
$ nextflow run main.nf -profile interactive --test --project <UPPMAX Project ID>
```

For more tests, see [Usage documentation](USAGE.md#test)

# Try the workflow on a desktop/laptop using Docker

You'll need to install [Nextflow](https://www.nextflow.io/) and  [Docker](https://www.docker.com/).

## Download all the Reference files
You will also need reference files.
We are currently working on an easier way to handle reference files, but this is the current process:

We are using the [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) __b37__ for most of the files.
The instructions are pretty clear on how to dowload files.

The following files need to be downloaded (and unzip):
- '1000G_phase1.indels.b37.vcf'
- '1000G_phase1.indels.b37.vcf.idx'
- 'dbsnp_138.b37.vcf'
- 'dbsnp_138.b37.vcf.idx'
- 'human_g1k_v37_decoy.dict.gz'
- 'human_g1k_v37_decoy.fasta.fai.gz'
- 'human_g1k_v37_decoy.fasta.gz'
- 'Mills_and_1000G_gold_standard.indels.b37.vcf'
- 'Mills_and_1000G_gold_standard.indels.b37.vcf.idx'

This file is on our repo so it will be easy ;-)
- '[centromeres.list](https://raw.githubusercontent.com/SciLifeLab/CAW/master/repeats/centromeres.list)'

These files are stored using GIT-LFS, on (CAW-References)[https://github.com/MaxUlysse/CAW-References], instructions on on how to set up git-lfs are on [git-lfs.github.com](https://git-lfs.github.com/).

- '1000G_phase3_20130502_SNP_maf0.3.loci'
- 'b37_cosmic_v74.noCHR.sort.4.1.vcf'
- 'b37_cosmic_v74.noCHR.sort.4.1.vcf.idx'

The last files are made when indexing `human_g1k_v37_decoy.fasta` with `bwa`.
- 'human_g1k_v37_decoy.fasta.pac'
- 'human_g1k_v37_decoy.fasta.amb'
- 'human_g1k_v37_decoy.fasta.ann'
- 'human_g1k_v37_decoy.fasta.bwt'
- 'human_g1k_v37_decoy.fasta.sa'

You can do that by using your own `bwa`:
```bash
bwa index human_g1k_v37_decoy.fasta
```

Or using a Docker image containing it:
```bash
docker run -v `pwd`:/tmp -w /tmp maxulysse/mapreads:1.0 bwa index human_g1k_v37_decoy.fasta
```
The next step would be to clone the repo, and edit `configuration/docker-test.config` to change `PathToReferences` with your own path to the reference files.

And then a:
```bash
nextflow run SciLifeLab/CAW --test -profile docker
```
Should work without a problem.
