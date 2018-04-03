# Installation

This small tutorial will explain to you how to install and run Sarek on a small sample test data on any POSIX compatible system (Linux, Solaris, OS X, etc).

To use this pipeline, you need to have a working version of Nextflow installed, Reference files and Docker or Singularity to facilitate the use of other tools. You can use a small reference genome as testing

- See the [Install Nextflow documentation](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- See the [Reference files documentation](REFERENCES.md)
- See the [Install Docker documentation](https://docs.docker.com/engine/installation/linux/ubuntu/#install-docker)
- See the [Install Singularity documentation](http://singularity.lbl.gov/install-linux)

## Installation

This workflow itself needs little installation.

Nextflow will automatically fetch Sarek from GitHub when launched if `SciLifeLab/Sarek` is specified as the workflow name.

Sarek use Singularity containers to package all the different tools.

If you plan to use the automatic pull of Singularity images, you can use the [`singularity.config`](../configuration/singularity.config) configuration file. You can also set up the Nextflow environnement variable `NXF_SINGULARITY_CACHEDIR` to choose where to store them.

For example
```bash
export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity
```

Docker can also be used as a container technology.

## Test Sarek with small dataset and small reference

The following tutorial explain how to run Sarek on a small dataset using a small reference.

```bash
# Connect to your system
> ssh -AX [USER]@[system]
# Or just open a terminal

# Install Nextflow
> curl -s https://get.nextflow.io | bash

# Move Nextflow into your $PATH
> mv nextflow $HOME/bin

# Set up the cache directory for Singularity Images if needed
> export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity

# make a test directory
> mkdir test_Sarek
> cd test_Sarek

# Download and build the smallGRCh37 reference using Singularity
> nextflow run SciLifeLab/Sarek/buildReferences.nf --download --genome smallGRCh37 -profile singularity

# Or download and build the smallGRCh37 reference using Docker
> nextflow run SciLifeLab/Sarek/buildReferences.nf --download --genome smallGRCh37 -profile docker

# Test Sarek on a test tiny set using Singularity
> nextflow run SciLifeLab/Sarek/main.nf --test --genome smallGRCh37 --noReports -profile singularity

# Or test Sarek on a test tiny set using Docker
> nextflow run SciLifeLab/Sarek/main.nf --test --genome smallGRCh37 --noReports -profile docker
```

## Update

To update Sarek, it's also very simple:


```bash
# Connect to your system
> ssh -AX [USER]@[system]REFERENCES
# Or just open a terminal

# Update Sarek
> nextflow pull SciLifeLab/Sarek
```

## Run Sarek on real data

Follow the [references documentation](REFERENCES.md) on how to download/build the references files.

Follow the [configuration and profile documentation](CONFIG.md) on how to modify and use the configuration files and profiles.

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
