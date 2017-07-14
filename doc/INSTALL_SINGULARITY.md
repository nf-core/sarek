# Installation using Singularity

To use this pipeline, you need to have a working version of Nextflow installed, References files and Singularity to facilitate the use of other tools. You can use a small reference genome as testing.

- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [Reference files documentation](REFERENCES.md)
- See the [Install Singularity documentation](http://singularity.lbl.gov/install-linux)

This workflow itself needs little installation. Clone the CAW repository and download/pull the singularity containers. And then you can directly use CAW on your computer.

```bash
# Make a test directory
mkdir test_CAW
cd test_CAW

# Download the singularity containers
nextflow run SciLifeLab/CAW-containers --singularity --containers bcftools,concatvcf,fastqc,freebayes,gatk,htslib,igvtools,mapreads,multiqc,picard,qualimap,runallelecount,runascat,runconvertallelecounts,runmanta,samtools,snpeffgrch37,snpeffgrch38,strelka,vepgrch37,vepgrch38 --singularityPublishDir containers/

# Build the smallGRCh37 reference
nextflow run SciLifeLab/CAW/buildReferences.nf --download --genome smallGRCh37
  -profile singularityTest

# Test the workflow on a test tiny set
nextflow run SciLifeLab/CAW --test --genome smallGRCh37 -profile singularityTest
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
