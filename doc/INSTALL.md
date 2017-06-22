# Installation

To use this pipeline, you need to have a working version of Nextflow installed, and References files.
You can use a small reference genome as testing.
Nextflow can also be use conjointly with Docker to facilitate the use of other tools.

- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [Reference files documentation](REFERENCES.md)

This small tutorial will explain to you how to run CAW on a small sample test data on Swedish UPPMAX clusters.
Some variables are specific, but it can be easily modified to suit any clusters.

## How to install guidelines

1. [On Milou](#on-milou)
2. [On Bianca](#on-bianca)
3. [With Docker](#with-docker)

## On Milou

For more information about Milou, follow the [Milou user guide](https://www.uppmax.uu.se/support/user-guides/milou-user-guide/).
This workflow itself needs no installation.
You just need to load the correct modules: `bioinfo-tools` and `Nextflow`.
Nextflow will automatically fetch CAW from GitHub when launched if `SciLifeLab/CAW` is specified as the workflow name.
So you can directly use CAW on Milou.

### Test CAW with small dataset and small reference

For more information, follow the [genomes files documentation](GENOMES.md).
The following tutorial explain how to run CAW on a small dataset using a small reference.

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

### Update CAW
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

### Use CAW with slurm
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

## On Bianca

For more information about Bianca, follow the [Bianca user guide](http://uppmax.uu.se/support/user-guides/bianca-user-guide/).
Bianca is made for sensitive data, so it's quite more complicated to access data from outside.
So CAW will have to be installed and updated manually.
You can either download CAW on your computer or on Milou, make an archive, and send it to Bianca using FileZilla or sftp given your preferences.


```bash
# Clone the repository
> git clone https://github.com/SciLifeLab/CAW.git
> cd CAW

# It is also possible to checkout a specific version using
> git checkout <branch, tag or commit>

# Make an archive to send to Bianca
> ./scripts/makeSnapshot.sh

# You will get this message in your terminal
Wrote CAW-[snapID].tar.gz

# Send the tar to Bianca (here using sftp)
# For FileZilla follow the Bianca user guide
> sftp [USER]-[PROJECT]@bianca-sftp.uppmax.uu.se:[USER]-[PROJECT]
> put CAW-[snapID].tar.gz

# The archive will be in the wharf folder in your user home on your Bianca project
s
# Connect to Bianca (Connect to Milou first if needed)
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your project
> cd /castor/project/proj_nobackup

# Make and go into a CAW directoy (where you will store all CAW versions)
> mkdir CAW
> cd CAW

# Copy the tar from wharf to the project
> cp /castor/project/proj_nobackup/wharf/[USER]/[USER]-[PROJECT]/CAW-[snapID].tgz /castor/project/proj_nobackup/CAW

# extract CAW
> tar -xvzf CAW-[snapID].tgz

# Make a symbolic link to the extracted repository
> ln -s CAW-[snapID] default
```

The idea is to have every member of your project to be able to use the same CAW version at the same time.
So every member of the project who wants to use CAW will need to do

```bash
# Connect to Bianca (Connect to Milou first if needed)
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your user directory
> cd /home/[USER]

# Make a symbolic link to the default CAW
> ln -s /castor/project/proj_nobackup/CAW/default CAW
```

### Update CAW

Repeat the same steps as for installing CAW, and once the tar has been extracted, you can replace the link.

```bash
# Connect to Bianca (Connect to Milou first if needed)
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to the CAW directory in your project
> cd /castor/project/proj_nobackup/CAW

# Remove link
> rm default

# Link to new CAW version
> ln -s CAW-[NEWsnapID] default
```

You can for example keep a `default` version that you are sure is working, an make a link for a `testing` or `development`


## With Docker

You will need to have Docker and Nextflow installed.

- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [Install Docker documentation](https://docs.docker.com/engine/installation/linux/ubuntu/#install-docker)

This workflow itself needs no installation. Nextflow will automatically fetch CAW from GitHub when launched if SciLifeLab/CAW is specified as the workflow name. So you can directly use CAW on your computer.

```bash
# Make a test directory
mkdir test_CAW
cd test_CAW

# Build the smallGRCh37 reference
nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37
  -profile docker

# Test the workflow on a test tiny set
nextflow run SciLifeLab/CAW --test --genome smallGRCh37 -profile docker
```

To update CAW, it's also very simple:
```bash
# Update CAW
nextflow pull SciLifeLab/CAW
```

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link] [![](images/NGI-final-small.png "NGI")][ngi-link]

[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: http://www.scilifelab.se/
