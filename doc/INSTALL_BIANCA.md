# Installation on `bianca`

This small tutorial will explain to you how to install and run Sarek on a small sample test data on the Swedish UPPMAX cluster `bianca` made for sensitive data.

For more information about `bianca`, follow the [`bianca` user guide](http://uppmax.uu.se/support/user-guides/bianca-user-guide/).
For more information about using Singularity with UPPMAX, follow the [Singularity UPPMAX guide](https://www.uppmax.uu.se/support-sv/user-guides/singularity-user-guide/).

Sarek use Singularity containers to package all the different tools.

As `bianca` is secure, no direct download is available, so Sarek and the Singularity containers will have to be installed and updated manually.

You can either download Sarek and the containers on your computer (you will need Nextflow and Singularity for that) or on `rackham`, make an archive, and send it to `bianca` using `FileZilla` or `sftp` given your preferences.

All Reference files are already stored in `bianca`.

```bash
# Connect to rackham
> ssh -AX [USER]@rackham.uppmax.uu.se
# Or just open a terminal

# Clone the repository
> git clone https://github.com/SciLifeLab/Sarek.git
> cd Sarek

# It is also possible to checkout a specific version using
> git checkout <branch, tag or commit>

# Use our script to make an archive to send to bianca
> ./scripts/makeSnapshot.sh

# You will get this message in your terminal
Wrote Sarek-[snapID].tar.gz

# Send the tar to bianca (here using sftp)
# For FileZilla follow the bianca user guide
> sftp [USER]-[PROJECT]@bianca-sftp.uppmax.uu.se:[USER]-[PROJECT]
> put Sarek-[snapID].tar.gz

# To get the containers
# This script will need Singularity and Nextflow installed
> ./scripts/do_all.sh --pull --tag <VERSION>

# Send the containers to bianca using the same method
# They will be in the containers/ directory as .img files

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

# extract Sarek
> tar -xvzf Sarek-[snapID].tgz

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
> nextflow run ~/Sarek/main.nf -profile slurm --project [PROJECT] ...
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

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
