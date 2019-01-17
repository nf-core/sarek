# Installation on `rackham`

This small tutorial will explain to you how to install and run Sarek on a small sample test data on the Swedish UPPMAX cluster `rackham`.

For more information about `rackham`, follow the [rackham user guide](https://www.uppmax.uu.se/support/user-guides/rackham-user-guide/).

This workflow itself needs little installation.

You need to install [Nextflow][nextflow-link] and put it somewhere in your `$PATH`

Sarek use Singularity containers to package all the different tools.

All containers, and all Reference files are already stored on UPPMAX.

Nextflow will automatically fetch Sarek from GitHub when launched if `SciLifeLab/Sarek` is specified as the workflow name.


## Test Sarek with small dataset and small reference

For more information, follow the [reference files documentation](REFERENCES.md).
The following tutorial explain how to run Sarek on a small dataset using a small reference.

```bash
# Connect to rackham
> ssh -AX [USER]@rackham.uppmax.uu.se

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
> mkdir test_Sarek
> cd test_Sarek

# Build the smallGRCh37 reference
> nextflow run SciLifeLab/Sarek/buildReferences.nf --download --genome smallGRCh37 --project [PROJECT]

# Test the workflow on a test tiny set
> nextflow run SciLifeLab/Sarek --test --genome smallGRCh37 --noReports --project [PROJECT]
```

## Update Sarek

```bash
# Connect to rackham
> ssh -AX [USER]@rackham.uppmax.uu.se

# Update Sarek
> nextflow pull SciLifeLab/Sarek
```

## Use Sarek with slurm

To use Sarek on rackham you will need to use the `slurm` profile.

```bash
# Connect to rackham
> ssh -AX [USER]@rackham.uppmax.uu.se

# Run the workflow directly on the login node
> nextflow run SciLifeLab/Sarek/main.nf --project [PROJECT] -profile slurm
```
