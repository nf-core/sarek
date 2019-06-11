# Installation on `bianca`

This small tutorial will explain to you how to install and run Sarek on a small sample test data on the Swedish UPPMAX cluster `bianca` made for sensitive data.

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
```

## Install Sarek

Sarek use Singularity containers to package all the different tools.
All containers, and all Reference files are already stored on UPPMAX.

As `bianca` is secure, no direct download is available, so Sarek will have to be installed and updated manually.

You can either download Sarek on your computer or on `rackham`, make an archive, and send it to `bianca` using `FileZilla` or `sftp` given your preferences.

```bash
# Connect to rackham
> ssh -AX [USER]@rackham.uppmax.uu.se
# Or just open a terminal

# Clone the repository
> git clone https://github.com/SciLifeLab/Sarek.git

# If you want to include the test data, you should use --recursive
> git clone --recursive https://github.com/SciLifeLab/Sarek.git

# Go to the newly created directory
> cd Sarek

# It is also possible to checkout a specific version using
> git checkout <branch, tag or commit>

# Use our script to make an archive to send to bianca
> ./scripts/makeSnapshot.sh

# Or you can also include the test data in this archive using git-archive-all
# Install pip
> curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
> python get-pip.py

# If it fails due to permission, you could consider using
> python get-pip.py --user

# Install git-archive-all using pip
> pip install git-archive-all
# If you used --user before, you might want to do that here too
> pip install git-archive-all --user
> ./scripts/makeSnapshot.sh --include-test-data

# You will get this message in your terminal
Wrote Sarek-[snapID].tar.gz

# Send the tar to bianca (here using sftp)
# For FileZilla follow the bianca user guide
> sftp [USER]-[PROJECT]@bianca-sftp.uppmax.uu.se:[USER]-[PROJECT]
> put Sarek-[snapID].tar.gz
> exit

# The archive will be in the wharf folder in your user home on your bianca project

# Connect to bianca
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your project
> cd /castor/project/proj_nobackup

# Make and go into a Sarek directoy (where you will store all Sarek versions)
> mkdir Sarek
> cd Sarek

# Copy the tar from wharf to the project
> cp /castor/project/proj_nobackup/wharf/[USER]/[USER]-[PROJECT]/Sarek-[snapID].tgz /castor/project/proj_nobackup/Sarek

# extract Sarek. Also remember to extract the containers you uploaded.
> tar xvzf Sarek-[snapID].tgz

# If you want other people to use it
# Be sure that your group has rights to the directory as well
> chown -R .[PROJECT] Sarek-[snapID]

# Make a symbolic link to the extracted repository
> ln -s Sarek-[snapID] default
```

The principle is to have every member of your project to be able to use the same Sarek version at the same time. So every member of the project who wants to use Sarek will need to do:

```bash
# Connect to bianca
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your user directory
> cd /home/[USER]

# Make a symbolic link to the default Sarek
> ln -s /castor/project/proj_nobackup/Sarek/default Sarek
```

And then Sarek can be used with:

```bash
> nextflow run ~/Sarek/main.nf -profile slurm --project [PROJECT] --genome [GENOME ASSEMBLY] --genome_base [PATH TO REFERENCE FILES] --containerPath [PATH TO CONTAINERS] ...
```

This is an example of how to run Sarek Somaic with the tool Ascat and the genome assembly version GRCh37:

```bash
> nextflow run ~/Sarek/somaticVC.nf -profile slurm --project [PROJECT] --tools ascat --sample [SAMPLE.TSV] --genome GRCh37 --genome_base /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37 --containerPath ~/Sarek/containers
```

## Update Sarek

Repeat the same steps as for installing Sarek, and once the tar has been extracted, you can replace the link.

```bash
# Connect to bianca (Connect to rackham first if needed)
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to the Sarek directory in your project
> cd /castor/project/proj_nobackup/Sarek

# Remove link
> rm default

# Link to new Sarek version
> ln -s Sarek-[NEWsnapID] default
```

You can for example keep a `default` version that you are sure is working, an make a link for a `testing` or `development`
