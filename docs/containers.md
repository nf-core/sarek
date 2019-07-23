# Containers

Our main container is designed using [Conda](https://conda.io/) to install all tools used in Sarek:

- [sarek](#sarek-)

For annotation, the main container can be used, but the cache has to be downloaded, or additional containers are available with cache (see [extra annotation documentation](annotation.md)):

- [sareksnpeff](#sareksnpeff-)
- [sarekvep](#sarekvep-)

## What is actually inside the containers

### sarek [![sarek-docker status][sarek-docker-badge]][sarek-docker-link]

- Based on `nfcore/base:latest`
- Contain **[AlleleCount][allelecount-link]** 4.0.2
- Contain **[BCFTools][bcftools-link]** 1.9
- Contain **[BWA][bwa-link]** 0.7.17
- Contain **[FastQC][fastqc-link]** 0.11.8
- Contain **[FreeBayes][freebayes-link]** 1.2.0
- Contain **[GATK4][gatk4-link]** 4.1.2.0
- Contain **[GeneSplicer][genesplicer-link]** 1.0
- Contain **[HTSlib][htslib-link]** 1.9
- Contain **[IGVtools][igvtools-link]** 2.3.93
- Contain **[Manta][manta-link]** 1.5.0
- Contain **[MultiQC][multiqc-link]** 1.7
- Contain **[Qualimap][qualimap-link]** 2.2.2b
- Contain **[R][r-link]** 3.5.1
- Contain **[RColorBrewer][rcolorbrewer-link]** 1.1
- Contain **[Rtracklayer][rtracklayer-link]** 1.42.1
- Contain **[samtools][samtools-link]** 1.9
- Contain **[snpEff][snpeff-link]** 4.3.1t
- Contain **[Strelka2][strelka-link]** 2.9.10
- Contain **[VCFanno][vcfanno-link]** 0.3.1
- Contain **[VCFtools][vcftools-link]** 0.1.16
- Contain **[VEP][vep-link]** 96.0

### sareksnpeff [![sareksnpeff-docker status][sareksnpeff-docker-badge]][sareksnpeff-docker-link]

- Based on `nfcore/base:latest`
- Contain **[snpEff][snpeff-link]** 4.3.1t
- Contain cache for `GRCh37`, `GRCh38`, or `GRCm38`

### sarekvep [![sarekvep-docker status][sarekvep-docker-badge]][sarekvep-docker-link]

- Based on `nfcore/base:latest`
- Contain **[GeneSplicer][genesplicer-link]** 1.0
- Contain **[VEP][vep-link]** 96.0
- Contain cache for `GRCh37`, `GRCh38`, or `GRCm38`

## Using helper script

An helper script, used for testing can also be used to help with pulling docker containers, or building singularity images.
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

Specify which release genome to use for annotation containers (`sareksnpeff`, `sarekvep`): `GRCh37`, `GRCh38`, `smallGRCh37`, `CanFan3.1`, `GRCm38`.
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

[allelecount-link]: https://github.com/cancerit/alleleCount
[bcftools-link]: https://github.com/samtools/bcftools
[bwa-link]: https://github.com/lh3/bwa
[fastqc-link]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[freebayes-link]: https://github.com/ekg/freebayes
[gatk4-link]: https://github.com/broadinstitute/gatk
[genesplicer-link]: https://ccb.jhu.edu/software/genesplicer/
[htslib-link]: https://github.com/samtools/htslib
[igvtools-link]: http://software.broadinstitute.org/software/igv/
[manta-link]: https://github.com/Illumina/manta
[multiqc-link]: https://github.com/ewels/MultiQC/
[qualimap-link]: http://qualimap.bioinfo.cipf.es
[r-link]: https://www.r-project.org/
[rcolorbrewer-link]: https://CRAN.R-project.org/package=RColorBrewer
[rtracklayer-link]: https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html
[samtools-link]: https://github.com/samtools/samtools
[sarek-docker-badge]: https://img.shields.io/docker/automated/nfcore/sarek.svg
[sarek-docker-link]: https://hub.docker.com/r/nfcore/sarek
[snpeff-link]: http://snpeff.sourceforge.net/
[sareksnpeff-docker-badge]: https://img.shields.io/docker/automated/nfcore/sareksnpeff.svg
[sareksnpeff-docker-link]: https://hub.docker.com/r/nfcore/sareksnpeff
[strelka-link]: https://github.com/Illumina/strelka
[vcfanno-link]: https://github.com/brentp/vcfanno
[vcftools-link]: https://vcftools.github.io/index.html
[vep-link]: https://github.com/Ensembl/ensembl-vep
[sarekvep-docker-badge]: https://img.shields.io/docker/automated/nfcore/sarekvep.svg
[sarekvep-docker-link]: https://hub.docker.com/r/nfcore/sarekvep
