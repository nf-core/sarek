# Installation using Docker

To use this pipeline, you need to have a working version of Nextflow installed, References files and Docker to facilitate the use of other tools.
You can use a small reference genome as testing.

- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [Reference files documentation](REFERENCES.md)
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
