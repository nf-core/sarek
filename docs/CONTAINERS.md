# Containers

Subsets of all containers can be downloaded:

- For processing, germline variant calling and Reports:
  - [sarek](#sarek-)
- For somatic variant calling, add:
  - [r-base](#r-base-)
  - [runallelecount](#runallelecount-)
- For annotation for GRCh37, you will need:
  - [snpeffgrch37](#snpeffgrch37-)
  - [vepgrch37](#vepgrch37-)
- For annotation for GRCh38, you will need:
  - [snpeffgrch38](#snpeffgrch38-)
  - [vepgrch38](#vepgrch38-)
- For annotation with cache use:
  - [sarek](#sarek-)

## Building

Use the Nextflow script to build and/or push containers from Docker and/or Singularity.

All the containers have built in UPPMAX directories, so there is no need to add them for use on UPPMAX clusters.
- See the [Singularity UPPMAX guide](https://www.uppmax.uu.se/support-sv/user-guides/singularity-user-guide/)

### Usage

```bash
nextflow run buildContainers.nf [--docker] [--singularity] /
[--containerPath <path>] [--push] [--containers <container1[,container2..]>] /
[--repository <repository>] [--tag tag]
```

- `--containers`: Choose which containers to build.
Default: `all`.
Possible values (to separate by commas):
  - `all` -  all available containers.
  - `r-base` - the [r-base](#r-base-) container.
  - `runallelecount` - the [runallelecount](#runallelecount-) container.
  - `sarek` - the [sarek](#sarek-) container.
  - `snpeffgrch37` - the [snpeffgrch37](#snpeffgrch37-) container.
  - `snpeffgrch38` - the [snpeffgrch38](#snpeffgrch38-) container.
  - `vepgrch37` - the [vepgrch37](#vepgrch37-) container.
  - `vepgrch38` - the [vepgrch38](#vepgrch38-) container.

- `--docker`: Build containers using `Docker`
- `--push`: Push containers to `DockerHub`
- `--repository`: Build containers under given repository.
Default: `maxulysse`
- `--singularity`: Build containers using `Singularity`.
- `--containerPath`: Select where to download containers.
Default: `$PWD`
- `--tag`: Build containers using given tag.
Default is version number.

### Example

```bash
nextflow run buildContainers.nf --docker --singularity --push --containers sarek
```

### For lazy users
We provide script to build/push or pull all containers
```bash
./scripts/do_all.sh        # Build all docker containers
./scripts/do_all.sh --push # Build and push all Docker containers into DockerHub
./scripts/do_all.sh --pull # Pull all containers from DockerHub with Singularity
```

## What is actually inside the containers

### r-base [![r-base-docker status][r-base-docker-badge]][r-base-docker-link]

 - Based on `r-base:3.3.2`
 - Contain **RColorBrewer**

### runallelecount [![runallelecount-docker status][runallelecount-docker-badge]][runallelecount-docker-link]

- Based on `debian:8.9`
- Contain **[AlleleCount][allelecount-link]** 2.2.0

### sarek [![sarek-docker status][sarek-docker-badge]][sarek-docker-link]

- Based on `nfcore/base:latest`
- Contain **[BCFTools][bcftools-link]** 1.8
- Contain **[BWA][bwa-link]** 0.7.17
- Contain **[FastQC][fastqc-link]** 0.11.8
- Contain **[FreeBayes][freebayes-link]** 1.2.0
- Contain **[GATK4][gatk4-link]** 4.0.9.0
- Contain **[HTSlib][htslib-link]** 1.9
- Contain **[IGVtools][igvtools-link]** 2.3.93
- Contain **[Manta][manta-link]** 1.4.0
- Contain **[MultiQC][multiqc-link]** 1.6
- Contain **[Qualimap][qualimap-link]** 2.2.2b
- Contain **[samtools][samtools-link]** 1.8
- Contain **[snpEff][snpeff-link]** 4.3.1t
- Contain **[Strelka2][strelka-link]** 2.9.3
- Contain **[VCFanno][vcfanno-link]** 0.3.0
- Contain **[VCFtools][vcftools-link]** 0.1.16
- Contain **[VEP][vep-link]** 95.1

### snpeffgrch37 [![snpeffgrch37-docker status][snpeffgrch37-docker-badge]][snpeffgrch37-docker-link]

- Based on `nfcore/base:latest`
- Contain **[snpEff][snpeff-link]** 4.3.1t
- Contain GRCh37.75

### snpeffgrch38 [![snpeffgrch38-docker status][snpeffgrch38-docker-badge]][snpeffgrch38-docker-link]

- Based on `nfcore/base:latest`
- Contain **[snpEff][snpeff-link]** 4.3.1t
- Contain GRCh38.86

### vepgrch37 [![vepgrch37-docker status][vepgrch37-docker-badge]][vepgrch37-docker-link]

- Based on `nfcore/base:latest`
- Contain **[VEP][vep-link]** 95.1
- Contain GRCh37

### vepgrch38 [![vepgrch38-docker status][vepgrch38-docker-badge]][vepgrch38-docker-link]

- Based on `nfcore/base:latest`
- Contain **[VEP][vep-link]** 95.1
- Contain GRCh38

[allelecount-link]: https://github.com/cancerit/alleleCount
[bcftools-link]: https://github.com/samtools/bcftools
[bwa-link]: https://github.com/lh3/bwa
[fastqc-link]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[freebayes-link]: https://github.com/ekg/freebayes
[gatk4-link]: https://github.com/broadinstitute/gatk
[htslib-link]: https://github.com/samtools/htslib
[igvtools-link]: http://software.broadinstitute.org/software/igv/
[manta-link]: https://github.com/Illumina/manta
[multiqc-link]: https://github.com/ewels/MultiQC/
[qualimap-link]: http://qualimap.bioinfo.cipf.es
[r-base-docker-badge]: https://img.shields.io/docker/automated/maxulysse/r-base.svg
[r-base-docker-link]: https://hub.docker.com/r/maxulysse/r-base
[rcolorbrewer-link]: https://CRAN.R-project.org/package=RColorBrewer
[runallelecount-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runallelecount.svg
[runallelecount-docker-link]: https://hub.docker.com/r/maxulysse/runallelecount
[samtools-link]: https://github.com/samtools/samtools
[sarek-docker-badge]: https://img.shields.io/docker/automated/maxulysse/sarek.svg
[sarek-docker-link]: https://hub.docker.com/r/maxulysse/sarek
[snpeff-link]: http://snpeff.sourceforge.net/
[snpeffgrch37-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeffgrch37.svg
[snpeffgrch37-docker-link]: https://hub.docker.com/r/maxulysse/snpeffgrch37
[snpeffgrch38-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeffgrch38.svg
[snpeffgrch38-docker-link]: https://hub.docker.com/r/maxulysse/snpeffgrch38
[strelka-link]: https://github.com/Illumina/strelka
[vcfanno-link]: https://github.com/brentp/vcfanno
[vcftools-link]: https://vcftools.github.io/index.html
[vep-link]: https://github.com/Ensembl/ensembl-vep
[vepgrch37-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch37.svg
[vepgrch37-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch37
[vepgrch38-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch38.svg
[vepgrch38-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch38
