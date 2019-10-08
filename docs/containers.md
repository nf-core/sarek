# Containers

Our main container is designed using [Conda](https://conda.io/) to install all tools used in Sarek:

- [sarek](#sarek-)

For annotation, the main container can be used, but the cache has to be downloaded, or additional containers are available with cache (see [extra annotation documentation](annotation.md)):

- [sareksnpeff](#sareksnpeff-)
- [sarekvep](#sarekvep-)

## What is actually inside the containers

### sarek [![sarek-docker status](https://img.shields.io/docker/automated/nfcore/sarek.svg)](https://hub.docker.com/r/nfcore/sarek)

- Based on `nfcore/base:latest`
- Contain **[ASCAT](https://github.com/Crick-CancerGenomics/ascat)** 2.5.2
- Contain **[AlleleCount](https://github.com/cancerit/alleleCount)** 4.0.2
- Contain **[BCFTools](https://github.com/samtools/bcftools)** 1.9
- Contain **[BWA](https://github.com/lh3/bwa)** 0.7.17
- Contain **[Control-FREEC](https://github.com/BoevaLab/FREEC)** 11.4
- Contain **[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** 0.11.8
- Contain **[FreeBayes](https://github.com/ekg/freebayes)** 1.2.0
- Contain **[GATK4](https://github.com/broadinstitute/gatk)** 4.1.2.0
- Contain **[GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)** 1.0
- Contain **[HTSlib](https://github.com/samtools/htslib)** 1.9
- Contain **[Manta](https://github.com/Illumina/manta)** 1.5.0
- Contain **[MultiQC](https://github.com/ewels/MultiQC/)** 1.7
- Contain **[Qualimap](http://qualimap.bioinfo.cipf.es)** 2.2.2b
- Contain **[samtools](https://github.com/samtools/samtools)** 1.9
- Contain **[snpEff](http://snpeff.sourceforge.net/)** 4.3.1t
- Contain **[Strelka2](https://github.com/Illumina/strelka)** 2.9.10
- Contain **[TIDDIT](https://github.com/SciLifeLab/TIDDIT)** 2.7.1
- Contain **[VCFanno](https://github.com/brentp/vcfanno)** 0.3.1
- Contain **[VCFtools](https://vcftools.github.io/index.html)** 0.1.16
- Contain **[VEP](https://github.com/Ensembl/ensembl-vep)** 95.2

### sareksnpeff [![sareksnpeff-docker status](https://img.shields.io/docker/automated/nfcore/sareksnpeff.svg)](https://hub.docker.com/r/nfcore/sareksnpeff)

- Based on `nfcore/base:latest`
- Contain **[snpEff](http://snpeff.sourceforge.net/)** 4.3.1t
- Contains cache for `GRCh37`, `GRCh38`, `GRCm38` or `CanFam3.1`

### sarekvep [![sarekvep-docker status](https://img.shields.io/docker/automated/nfcore/sarekvep.svg)](https://hub.docker.com/r/nfcore/sarekvep)

- Based on `nfcore/base:latest`
- Contain **[GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)** 1.0
- Contain **[VEP](https://github.com/Ensembl/ensembl-vep)** 95.2
- Contain cache for `GRCh37`, `GRCh38`, `GRCm38` or `CanFam3.1`

## Using helper script

A helper script, used for testing can also be used to help with pulling docker containers, or building singularity images.
The following parameters can be used:

### Engine: -n

Specify which container engine to use: `docker` or `singularity`.
Default:`docker`

### Containers: -c

Specify which containers to build: `SNPEFF`, `VEP` or `ALL`.
Default:`ALL`

### Version: -T

Specify which release to pull or build: any tagged release, or `dev`.
Default:`dev`

### Genome: -g

Specify which release genome to use for annotation containers (`sareksnpeff`, `sarekvep`): `smallGRCh37`, `GRCh37`, `GRCh38`, `GRCm38` or `CanFam3.1`.
Default:`smallGRCh37`

### Singularity

To specify where to build singularity image, use the Nextflow ENV variable `NXF_SINGULARITY_CACHEDIR`, ie:

```bash
NXF_SINGULARITY_CACHEDIR=/data/singularity ./scripts/download_image.sh -n singularity -t ALL -T dev -g GRCh38
```

That will build the main container, plus the annotation containers (`sareksnpeff`, `sarekvep`) for `GRCh38`, in the `/data/singularity` folder.

## Building your own

Our containers are designed using [Conda](https://conda.io/).
The `environment.yml` file can easilly be modified if particular versions of tools are more suited to your needs.
