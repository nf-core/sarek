# Installation on a secure cluster

This small tutorial will explain to you how to install and run nf-core/sarek on a small sample test data on the Swedish UPPMAX cluster `bianca` made for sensitive data.
It can be followed to install on any similar secure cluster.

For more information about `bianca`, follow the [`bianca` user guide](http://uppmax.uu.se/support/user-guides/bianca-user-guide/).
For more information about using Singularity with UPPMAX, follow the [Singularity UPPMAX guide](https://www.uppmax.uu.se/support-sv/user-guides/singularity-user-guide/).

## Install Nextflow

```bash
# Connect to rackham
> ssh -AX [USER]@rackham.uppmax.uu.se
# Or just open a terminal

# Download the all Nextflow bundle
> wget https://github.com/nextflow-io/nextflow/releases/download/v[xx.yy.zz]/nextflow-[xx.yy.zz]-all

# Send to bianca (here using sftp)
# For FileZilla follow the bianca user guide
> sftp [USER]-[PROJECT]@bianca-sftp.uppmax.uu.se:[USER]-[PROJECT]
> put nextflow-[xx.yy.zz]-all

# Exit sftp
> exit

# Connect to bianca
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your project
> cd /castor/project/proj_nobackup

# Make directory for Nextflow
> mkdir tools
> mkdir tools/nextflow

# Move Nextflow from wharf to its directory
> mv /castor/project/proj_nobackup/wharf/[USER]/[USER]-[PROJECT]/nextflow-[xx.yy.zz]-all /castor/project/proj_nobackup/tools/nextflow

# Establish permission
> chmod a+x /castor/project/proj_nobackup/tools/nextflow/nextflow-[xx.yy.zz]-all

# If you want other people to use it
# Be sure that your group has rights to the directory as well

> chown -R .[PROJECT] /castor/project/proj_nobackup/tools/nextflow/nextflow-[xx.yy.zz]-all

# Make a link to it
> ln -s /castor/project/proj_nobackup/tools/nextflow/nextflow-[xx.yy.zz]-all /castor/project/proj_nobackup/tools/nextflow/nextflow

# And everytime you're launching Nextflow, don't forget to export the following ENV variables
# Or add them to your .bashrc file
> export NXF_HOME=/castor/project/proj/nobackup/tools/nextflow/
> export PATH=${NXF_HOME}:${PATH}
> export NXF_TEMP=$SNIC_TMP
> export NXF_LAUNCHER=$SNIC_TMP
> export NXF_SINGULARITY_CACHEDIR=/sw/data/uppnex/ToolBox/sarek
```

## Install nf-core/sarek

nf-core/sarek use Singularity containers to package all the different tools.
All containers, and all Reference files are already stored on UPPMAX.

As `bianca` is secure, no direct download is available, so nf-core/sarek will have to be installed and updated manually.

You can either download nf-core/sarek on your computer or on `rackham`, make an archive, and send it to `bianca` using `FileZilla` or `sftp` given your preferences.

```bash
# Connect to rackham
> ssh -AX [USER]@rackham.uppmax.uu.se
# Or just open a terminal

# Clone the repository
> git clone https://github.com/nf-core/sarek.git

# Go to the newly created directory
> cd sarek

# It is also possible to checkout a specific version using
> git checkout <branch, tag or commit>

# You also include the nf-core/test-datasets and nf-core/configs using git-archive-all
# Install pip
> curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
> python get-pip.py

# If it fails due to permission, you could consider using
> python get-pip.py --user

# Install git-archive-all using pip
> pip install git-archive-all
# If you used --user before, you might want to do that here too
> pip install git-archive-all --user
> ./scripts/make_snapshot.sh --include-test-data --include-configs

# Or you can just include nf-core/sarek:
> ./scripts/make_snapshot.sh

# You will get this message in your terminal
Wrote sarek-[snapID].tar.gz

# Send the archive to bianca (here using sftp)
# For FileZilla follow the bianca user guide
> sftp [USER]-[PROJECT]@bianca-sftp.uppmax.uu.se:[USER]-[PROJECT]
> put sarek-[snapID].tar.gz
> exit

# The archive will be in the wharf folder in your user home on your bianca project

# Connect to bianca
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your project
> cd /castor/project/proj_nobackup

# Make and go into a nf-core/sarek directoy (where you will store all nf-core/sarek versions)
> mkdir sarek
> cd sarek

# Copy the tar from wharf to the project
> cp /castor/project/proj_nobackup/wharf/[USER]/[USER]-[PROJECT]/sarek-[snapID].tgz /castor/project/proj_nobackup/sarek

# extract the archive. Also remember to extract the containers you uploaded.
> tar xvzf sarek-[snapID].tgz

# If you want other people to use it,
# Be sure that your group has rights to the directory as well
> chown -R .[PROJECT] sarek-[snapID]

# Make a symbolic link to the extracted repository
> ln -s sarek-[snapID] default
```

The principle is to have every member of your project to be able to use the same nf-core/sarek version at the same time.
So every member of the project who wants to use nf-core/sarek will need to do:

```bash
# Connect to bianca
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your user directory
> cd /home/[USER]

# Make a symbolic link to the default nf-core/sarek
> ln -s /castor/project/proj_nobackup/sarek/default sarek
```

Singularity images for Sarek are available on Uppmax in `/sw/data/uppnex/ToolBox/sarek`.
Sometimes Nextflow needs write access to the image folder, and if so the images needs to be copied to a location with write permission, for example in a subfolder of your project folder.

```bash
#Create a folder for the singularity images somewhere in your project:
mkdir sarek_simg

#Copy the relevant singularity image from the write protected folder on Uppmax to the folder where you have write permission:
cp /sw/data/uppnex/ToolBox/sarek/nfcore-sarek-dev.img /path/to/your/sarek_simg/.

#Update the ENV parameter NXF_SINGULARITY_CACHEDIR
export NXF_SINGULARITY_CACHEDIR=/path/to/your/sarek_simg
```

And then nf-core/sarek can be used with:

```bash
> nextflow run ~/sarek/main.nf -profile uppmax --custom_config_base ~/sarek/configs --project [PROJECT] --genome [GENOME ASSEMBLY] ...
```

This command worked on Bianca 20190906:

```bash
>screen -S SAMPLE /path/to/nextflow run /path/to/sarek/main.nf -profile uppmax --project PROJID --sample SAMPLE.tsv --genome GRCh37 --genomes_base /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37  --step variantcalling --tools ASCAT --igenomesIgnore

#To detach screen:
ctrl-A-D
```

This is an example of how to run sarek with the tool Manta and the genome assembly version GRCh38:

```bash
> nextflow run ~/sarek/main.nf -profile uppmax --custom_config_base ~/sarek/configs --project [PROJECT] --tools Manta --sample [SAMPLE.TSV] --genome GRCh38
```

## Update nf-core/sarek

Repeat the same steps as for installing nf-core/sarek, and once the tar has been extracted, you can replace the link.

```bash
# Connect to bianca (Connect to rackham first if needed)
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to the sarek directory in your project
> cd /castor/project/proj_nobackup/sarek

# Remove link
> rm default

# Link to new nf-core/sarek version
> ln -s sarek-[NEWsnapID] default
```

You can for example keep a `default` version that you are sure is working, an make a link for a `testing` or `development`
