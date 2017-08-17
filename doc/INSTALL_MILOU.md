# Installation on Milou

This small tutorial will explain to you how to run CAW on a small sample test data on the Swedish UPPMAX cluster Milou.

Some variables are specific, but it can be easily modified to suit any clusters.

For more information about Milou, follow the [Milou user guide](https://www.uppmax.uu.se/support/user-guides/milou-user-guide/).

This workflow itself needs no installation.

You just need to load the correct modules: `bioinfo-tools` and `Nextflow`.

The Reference files are already stored in Milou.

Nextflow will automatically fetch CAW from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.

So you can directly use CAW on Milou.

## Test CAW with small dataset and small reference

For more information, follow the [genomes files documentation](GENOMES.md). The following tutorial explain how to run CAW on a small dataset using a small reference.

```bash
# Connect to Milou
ssh -AX [USER]@milou.uppmax.uu.se

# Load modules
module load bioinfo-tools Nextflow

# Set up this two variables as asked by UPPMAX
export NXF_LAUNCHBASE=$SNIC_TMP
export NXF_TEMP=$SNIC_TMP

# make a test directory
mkdir test_CAW
cd test_CAW

# Connect to an interactive session
$ interactive -A [PROJECT] -p node

# Build the smallGRCh37 reference
nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37 --project [PROJECT]

# Test the workflow on a test tiny set
nextflow run SciLifeLab/CAW --test --genome smallGRCh37 --project [PROJECT]
```

## Update CAW

```bash
# Connect to Milou
ssh -AX [USER]@milou.uppmax.uu.se

# Load modules
module load bioinfo-tools Nextflow

# Set up this two variables as asked by UPPMAX
export NXF_LAUNCHBASE=$SNIC_TMP
export NXF_TEMP=$SNIC_TMP

# Update CAW
nextflow pull SciLifeLab/CAW
```

## Use CAW with slurm

To use CAW on Milou you will need to use the `slurm` profile.

```bash
# Connect to Milou
ssh -AX [USER]@milou.uppmax.uu.se

# Load modules
module load bioinfo-tools Nextflow

# Set up this two variables as asked by UPPMAX
export NXF_LAUNCHBASE=$SNIC_TMP
export NXF_TEMP=$SNIC_TMP

# Run the workflow directly on the login node
nextflow run SciLifeLab/CAW --sample [FILE.TSV] --genome [GENOME] --project [PROJECT] -profile slurm
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
