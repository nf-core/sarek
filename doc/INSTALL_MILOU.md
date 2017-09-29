# Installation on `milou`

This small tutorial will explain to you how to install and run CAW on a small sample test data on the Swedish UPPMAX cluster `milou`.

For more information about `milou`, follow the [milou user guide](https://www.uppmax.uu.se/support/user-guides/milou-user-guide/).

This workflow itself needs little installation.

You need to install [Nextflow][nextflow-link] and put it somewhere in your `$PATH`

The Reference files are already stored in `milou`.

Nextflow will automatically fetch CAW from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.

CAW use Singularity containers to package all the different tools.

On `milou` you do have the possibility to use the automatic pull of the containers.

You can choose a specific location to store these, otherwise they will be stored in the directory where you're running CAW.

## Test CAW with small dataset and small reference

For more information, follow the [reference files documentation](REFERENCES.md). The following tutorial explain how to run CAW on a small dataset using a small reference.

```bash
# Connect to milou
> ssh -AX [USER]@milou.uppmax.uu.se

# Install Nextflow
> cd $HOME
> curl -s https://get.nextflow.io | bash

# Move Nextflow into your $PATH
> mv nextflow $HOME/bin

# Connect to an interactive session
> interactive -A [PROJECT] -p node

# Set up this two variables as asked by UPPMAX
> export NXF_LAUNCHBASE=$SNIC_TMP
> export NXF_TEMP=$SNIC_TMP

# Set up the cache directory for Singularity Images
> export NXF_SINGULARITY_CACHEDIR=$HOME/.singularity

# make a test directory
> mkdir test_CAW
> cd test_CAW

# Build the smallGRCh37 reference
> nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37 --project [PROJECT] -profile download

# Test the workflow on a test tiny set
> nextflow run SciLifeLab/CAW --test --genome smallGRCh37 --noReports --project [PROJECT] -profile download
```

## Update CAW

```bash
# Connect to Milou
> ssh -AX [USER]@milou.uppmax.uu.se

# Update CAW
> nextflow pull SciLifeLab/CAW
```

## Use CAW with slurm

To use CAW on Milou you will need to use the `slurm` profile.

```bash
# Connect to Milou
> ssh -AX [USER]@milou.uppmax.uu.se

# Run the workflow directly on the login node
> nextflow run SciLifeLab/CAW --sample [FILE.TSV] --genome [GENOME] --project [PROJECT] -profile slurm-download
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
