# Installation on Bianca

This small tutorial will explain to you how to run CAW on a small sample test data on the Swedish UPPMAX cluster Bianca made for sensitive data.

Some variables are specific, but it can be easily modified to suit any clusters.

For more information about Bianca, follow the [Bianca user guide](http://uppmax.uu.se/support/user-guides/bianca-user-guide/).

As Bianca is secure, no direct download is available, so CAW will have to be installed and updated manually.

You can either download CAW on your computer or on Milou, make an archive, and send it to Bianca using FileZilla or sftp given your preferences.

All Reference files are already stored in Bianca.

```bash
# Clone the repository
> git clone https://github.com/SciLifeLab/CAW.git
> cd CAW

# It is also possible to checkout a specific version using
> git checkout <branch, tag or commit>

# Use our script to make an archive to send to Bianca
> ./scripts/makeSnapshot.sh

# You will get this message in your terminal
Wrote CAW-[snapID].tar.gz

# Send the tar to Bianca (here using sftp)
# For FileZilla follow the Bianca user guide
> sftp [USER]-[PROJECT]@bianca-sftp.uppmax.uu.se:[USER]-[PROJECT]
> put CAW-[snapID].tar.gz

# The archive will be in the wharf folder in your user home on your Bianca project

# Connect to Bianca
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

The principle is to have every member of your project to be able to use the same CAW version at the same time. So every member of the project who wants to use CAW will need to do:

```bash
# Connect to Bianca
> ssh -A [USER]-[PROJECT]@bianca.uppmax.uu.se

# Go to your user directory
> cd /home/[USER]

# Make a symbolic link to the default CAW
> ln -s /castor/project/proj_nobackup/CAW/default CAW
```

And then CAW can be used with:

```bash
nextflow run ~/CAW/main.nf ...
```

## Update CAW

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

--------------------------------------------------------------------------------

[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
