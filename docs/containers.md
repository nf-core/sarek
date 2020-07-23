# Containers

Our main container is designed using [Conda](https://conda.io/) to install all tools used in Sarek:

- [sarek](#sarek-)

For annotation, the main container can be used, but the cache has to be downloaded, or additional containers are available with cache (see [extra annotation documentation](annotation.md)):

- [sareksnpeff](#sareksnpeff-)
- [sarekvep](#sarekvep-)

## What is actually inside the containers

### sarek [![sarek-docker status](https://img.shields.io/docker/automated/nfcore/sarek.svg)](https://hub.docker.com/r/nfcore/sarek)

- Based on `nfcore/base:1.9`
- Contain **[ASCAT](https://github.com/Crick-CancerGenomics/ascat)** 2.5.2
- Contain **[AlleleCount](https://github.com/cancerit/alleleCount)** 4.0.2
- Contain **[BCFTools](https://github.com/samtools/bcftools)** 1.9
- Contain **[BWA](https://github.com/lh3/bwa)** 0.7.17
- Contain **[Control-FREEC](https://github.com/BoevaLab/FREEC)** 11.5
- Contain **[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** 0.11.9
- Contain **[FreeBayes](https://github.com/ekg/freebayes)** 1.3.2
- Contain **[GATK4-spark](https://github.com/broadinstitute/gatk)** 4.1.6.0
- Contain **[GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)** 1.0
- Contain **[ggplot2](https://github.com/tidyverse/ggplot2)** 3.3.0
- Contain **[HTSlib](https://github.com/samtools/htslib)** 1.9
- Contain **[Manta](https://github.com/Illumina/manta)** 1.6.0
- Contain **[msisensor](https://github.com/ding-lab/msisensor)** 0.5
- Contain **[MultiQC](https://github.com/ewels/MultiQC/)** 1.8
- Contain **[Qualimap](http://qualimap.bioinfo.cipf.es)** 2.2.2d
- Contain **[samtools](https://github.com/samtools/samtools)** 1.9
- Contain **[snpEff](http://snpeff.sourceforge.net/)** 4.3.1t
- Contain **[Strelka2](https://github.com/Illumina/strelka)** 2.9.10
- Contain **[TIDDIT](https://github.com/SciLifeLab/TIDDIT)** 2.7.1
- Contain **[pigz](https://zlib.net/pigz/)** 2.3.4
- Contain **[Trim Galore](https://github.com/FelixKrueger/TrimGalore)** 0.6.5
- Contain **[VCFanno](https://github.com/brentp/vcfanno)** 0.3.2
- Contain **[VCFtools](https://vcftools.github.io/index.html)** 0.1.16
- Contain **[VEP](https://github.com/Ensembl/ensembl-vep)** 99.2

### sareksnpeff [![sareksnpeff-docker status](https://img.shields.io/docker/automated/nfcore/sareksnpeff.svg)](https://hub.docker.com/r/nfcore/sareksnpeff)

- Based on `nfcore/base:1.9`
- Contain **[snpEff](http://snpeff.sourceforge.net/)** 4.3.1t
- Contains cache for `GRCh37`, `GRCh38`, `GRCm38`, `CanFam3.1` or `WBcel235`

### sarekvep [![sarekvep-docker status](https://img.shields.io/docker/automated/nfcore/sarekvep.svg)](https://hub.docker.com/r/nfcore/sarekvep)

- Based on `nfcore/base:1.9`
- Contain **[GeneSplicer](https://ccb.jhu.edu/software/genesplicer/)** 1.0
- Contain **[VEP](https://github.com/Ensembl/ensembl-vep)** 99.2
- Contain cache for `GRCh37`, `GRCh38`, `GRCm38`, `CanFam3.1` or `WBcel235`

## Building your own

Our containers are designed using [Conda](https://conda.io/).
The [`environment.yml`](../environment.yml) file can be modified if particular versions of tools are more suited to your needs.

The following commands can be used to build/download containers on your own system:

- Adjust `VERSION` for sarek version (typically a release or `dev`).

### Build with Conda

```Bash
conda env create -f environment.yml
```

### Build with Docker

- `sarek`

```Bash
docker build -t nfcore/sarek:<VERSION> .
```

- `sareksnpeff`

Adjust arguments for `GENOME` version and snpEff `CACHE_VERSION`

```Bash
docker build -t nfcore/sareksnpeff:<VERSION>.<GENOME> containers/snpeff/. --build-arg GENOME=<GENOME> --build-arg CACHE_VERSION=<CACHE_VERSION>
```

- `sarekvep`

Adjust arguments for `GENOME` version, `SPECIES` name and VEP `VEP_VERSION`

```Bash
docker build -t nfcore/sarekvep:<VERSION>.<GENOME> containers/vep/. --build-arg GENOME=<GENOME> --build-arg SPECIES=<SPECIES> --build-arg VEP_VERSION=<VEP_VERSION>
```

### Pull with Docker

- `sarek`

```Bash
docker pull nfcore/sarek:<VERSION>
```

- `sareksnpeff`

Adjust arguments for `GENOME` version

```Bash
docker pull nfcore/sareksnpeff:<VERSION>.<GENOME>
```

- `sarekvep`

Adjust arguments for `GENOME` version

```Bash
docker pull nfcore/sarekvep:<VERSION>.<GENOME>
```

### Pull with Singularity

You can directly pull singularity image, in the path used by the Nextflow ENV variable `NXF_SINGULARITY_CACHEDIR`, ie:

```Bash
cd $NXF_SINGULARITY_CACHEDIR
singularity build ...
```

- `sarek`

```Bash
singularity build nfcore-sarek-<VERSION>.img docker://nfcore/sarek:<VERSION>
```

- `sareksnpeff`

Adjust arguments for `GENOME` version

```Bash
singularity build nfcore-sareksnpeff-<VERSION>.<GENOME>.img docker://nfcore/sareksnpeff:<VERSION>.<GENOME>
```

- `sarekvep`

Adjust arguments for `GENOME` version

```Bash
singularity build nfcore-sarekvep-<VERSION>.<GENOME>.img docker://nfcore/sarekvep:<VERSION>.<GENOME>
```
