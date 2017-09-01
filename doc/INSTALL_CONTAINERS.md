# Installation using containers (Docker or Singularity)

To use this pipeline, you need to have a working version of Nextflow installed, References files and Docker or Singularity to facilitate the use of other tools. You can use a small reference genome as testing.

- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [Reference files documentation](REFERENCES.md)
- See the [Install Docker documentation](https://docs.docker.com/engine/installation/linux/ubuntu/#install-docker)
- See the [Install Singularity documentation](http://singularity.lbl.gov/install-linux)

## Installation

This workflow itself needs no installation. Nextflow will automatically fetch CAW from GitHub when launched if SciLifeLab/CAW is specified as the workflow name. So you can directly use CAW on your own machine.

```bash
# Make a test directory
mkdir test_CAW
cd test_CAW

# Build the smallGRCh37 reference using Docker
nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37 -profile dockerTest

# Test the workflow on a test tiny set using Docker
nextflow run SciLifeLab/CAW --test --genome smallGRCh37 -profile dockerTest

# Build the smallGRCh37 reference using Singularity
nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37 -profile singularityTest

# Test the workflow on a test tiny set using Singularity
nextflow run SciLifeLab/CAW --test --genome smallGRCh37 -profile singularityTest
```

## Update

To update CAW, it's also very simple:

```bash
# Update CAW
nextflow pull SciLifeLab/CAW
```

## Use Singularity

If you plan to use the automatic pull of Singularity images, you can use the [`singularity.config`](../configuration/singularity.config) configuration file. You can also set up the Nextflow environnement variable `NXF_SINGULARITY_CACHEDIR` to choose where to store them.

For example
```bash
export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity
```

If needed (for example, if no direct internet access on a secure machine), you can specify a pathway to Singularity containers that you have downloaded before hand and transfered on the right machine. You can specify the directory to download the containers with `--singularityPublishDir`.

```
# Download the singularity containers
nextflow run SciLifeLab/CAW-containers --singularity --containers bcftools,concatvcf,fastqc,freebayes,gatk,htslib,igvtools,mapreads,multiqc,mutect1,picard,qualimap,runallelecount,runascat,runconvertallelecounts,runmanta,samtools,snpeffgrch37,snpeffgrch38,strelka,vepgrch37,vepgrch38 --singularityPublishDir containers/
nextflow run SciLifeLab/CAW-containers --singularity --containers gatk --tag 1.0 --singularityPublishDir containers/
```

And then you can use the [`singularity-download.config`](../configuration/singularity-download.config) configuration file.

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
