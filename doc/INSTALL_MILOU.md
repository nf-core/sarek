# Installation on Milou

This small tutorial will explain to you how to run CAW on a small sample test data on the Swedish UPPMAX cluster Milou.

Some variables are specific, but it can be easily modified to suit any clusters.

For more information about Milou, follow the [Milou user guide](https://www.uppmax.uu.se/support/user-guides/milou-user-guide/).

This workflow itself needs little installation.

You just need install [Nextflow][nextflow-link] and put it somewhere in your `$PATH`

The Reference files are already stored in Milou.

Nextflow will automatically fetch CAW from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.

So you can directly use CAW on Milou.

## Test CAW with small dataset and small reference

For more information, follow the [genomes files documentation](GENOMES.md). The following tutorial explain how to run CAW on a small dataset using a small reference.

```bash
# Connect to Milou
ssh -AX [USER]@milou.uppmax.uu.se

# Install Nextflow
cd $HOME
curl -s https://get.nextflow.io | bash
mv nextflow $HOME/bin

# Connect to an interactive session
$ interactive -A [PROJECT] -p node

# Set up this two variables as asked by UPPMAX
export NXF_LAUNCHBASE=$SNIC_TMP
export NXF_TEMP=$SNIC_TMP

# make a test directory
mkdir test_CAW
cd test_CAW

# Build the smallGRCh37 reference
nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37 --project [PROJECT]

# Test the workflow on a test tiny set
nextflow run SciLifeLab/CAW --test --genome smallGRCh37 --noReports --project [PROJECT]
```

## Update CAW

```bash
# Connect to Milou
ssh -AX [USER]@milou.uppmax.uu.se

# Update CAW
nextflow pull SciLifeLab/CAW
```

## Use CAW with slurm

To use CAW on Milou you will need to use the `slurm` profile.

```bash
# Connect to Milou
ssh -AX [USER]@milou.uppmax.uu.se

# Run the workflow directly on the login node
nextflow run SciLifeLab/CAW --sample [FILE.TSV] --genome [GENOME] --project [PROJECT] -profile slurm
```

## Use CAW with slurm and Singularity

- See the [Install with containers documentation](INSTALL_CONTAINERS.md)
To use CAW on Milou you will need to use the `singularitySlurm` profile.

```bash
# Connect to Milou
ssh -AX [USER]@milou.uppmax.uu.se

#
# Download the singularity containers
nextflow run SciLifeLab/CAW-containers --singularity --containers bcftools,concatvcf,fastqc,freebayes,gatk,htslib,igvtools,mapreads,multiqc,mutect1,picard,qualimap,runallelecount,runascat,runconvertallelecounts,runmanta,samtools,snpeffgrch37,snpeffgrch38,strelka,vepgrch37,vepgrch38 --singularityPublishDir containers/
nextflow run SciLifeLab/CAW-containers --singularity --containers gatk --tag 1.0 --singularityPublishDir containers/

# Run the workflow directly on the login node
nextflow run SciLifeLab/CAW --sample [FILE.TSV] --genome [GENOME] --project [PROJECT] -profile singularitySlurm
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
