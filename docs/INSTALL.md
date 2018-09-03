# Installation

This small tutorial will explain to you how to install and run Sarek on a small sample test data on any POSIX compatible system (Linux, Solaris, OS X, etc).

To use this pipeline, you need to have a working version of Nextflow installed, Reference files and Docker or Singularity as a container engine.
You can use a small reference genome as testing.

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

You can [Test Sarek with small dataset and small reference](https://github.com/SciLifeLab/Sarek/blob/master/docs/TESTS.md)

## Update

To update Sarek, it's also very simple:


```bash
# Connect to your system
> ssh -AX [USER]@[system]REFERENCES
# Or just open a terminal

# Update Sarek
> nextflow pull SciLifeLab/Sarek
```

## Running Sarek with real data

Follow the [references documentation](REFERENCES.md) on how to download/build the references files.

Follow the [configuration and profile documentation](CONFIG.md) on how to modify and use the configuration files and profiles.
