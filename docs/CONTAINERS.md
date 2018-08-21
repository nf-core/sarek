# Containers

Subsets of all containers can be dowloaded:

For processing, germline and somatic variant calling and Reports:
 - [sarek](#sarek-)

For annotation for GRCh37, you will need:
 - [snpeffgrch37](#snpeffgrch37-)
 - [vepgrch37](#vepgrch37-)

For annotation for GRCh38, you will need:
 - [snpeffgrch38](#snpeffgrch38-)
 - [vepgrch38](#vepgrch38-)

## r-base [![r-base-docker status][r-base-docker-badge]][r-base-docker-link]

 - Based on `debian:8.9`
 - Contain **[AlleleCount][allelecount-link]** 2.2.0

## runallelecount [![runallelecount-docker status][runallelecount-docker-badge]][runallelecount-docker-link]

- Based on `debian:8.9`
- Contain **[AlleleCount][allelecount-link]** 2.2.0

## sarek [![sarek-docker status][sarek-docker-badge]][sarek-docker-link]

- Based on `debian:8.9`
- Contain **[BCFTools][bcftools-link]** 1.5
- Contain **[BWA][bwa-link]** 0.7.16
- Contain **[HTSlib][htslib-link]** 1.5
- Contain **[Manta][manta-link]** 1.1.1
- Contain **[samtools][samtools-link]** 1.5
- Contain **[Strelka][strelka-link]** 2.8.2

## snpeff [![snpeff-docker status][snpeff-docker-badge]][snpeff-docker-link]

- Based on `openjdk:8-slim`
- Contain **[snpEff][snpeff-link]** 4.3i

## snpeffgrch37 [![snpeffgrch37-docker status][snpeffgrch37-docker-badge]][snpeffgrch37-docker-link]

- Based on `maxulysse/snpeff`
- Contain **[snpEff][snpeff-link]** 4.3i
- Contain GRCh37.75

## snpeffgrch38 [![snpeffgrch38-docker status][snpeffgrch38-docker-badge]][snpeffgrch38-docker-link]

- Based on `maxulysse/snpeff`
- Contain **[snpEff][snpeff-link]** 4.3i
- Contain GRCh38.86

## vepgrch37 [![vepgrch37-docker status][vepgrch37-docker-badge]][vepgrch37-docker-link]

- Based on `willmclaren/ensembl-vep:release_90.6`
- Contain **[VEP][vep-link]** 90.5
- Contain GRCh37

## vepgrch38 [![vepgrch38-docker status][vepgrch38-docker-badge]][vepgrch38-docker-link]

- Based on `willmclaren/ensembl-vep:release_90.6`
- Contain **[VEP][vep-link]** 90.5
- Contain GRCh38

---
[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[allelecount-link]: https://github.com/cancerit/alleleCount
[bcftools-link]: https://github.com/samtools/bcftools
[bwa-link]: https://github.com/lh3/bwa
[fastqc-link]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[freebayes-link]: https://github.com/ekg/freebayes
[gatk-link]: https://github.com/broadgsa/gatk-protected
[htslib-link]: https://github.com/samtools/htslib
[igvtools-link]: http://software.broadinstitute.org/software/igv/
[manta-link]: https://github.com/Illumina/manta
[multiqc-link]: https://github.com/ewels/MultiQC/
[mutect1-link]: https://github.com/broadinstitute/mutect
[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[picard-link]: https://github.com/broadinstitute/picard
[qualimap-link]: http://qualimap.bioinfo.cipf.es
[rcolorbrewer-link]: https://CRAN.R-project.org/package=RColorBrewer
[runallelecount-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runallelecount.svg
[runallelecount-docker-link]: https://hub.docker.com/r/maxulysse/runallelecount
[r-base-docker-badge]: https://img.shields.io/docker/automated/maxulysse/r-base.svg
[r-base-docker-link]: https://hub.docker.com/r/maxulysse/r-base
[samtools-link]: https://github.com/samtools/samtools
[sarek-docker-badge]: https://img.shields.io/docker/automated/maxulysse/sarek.svg
[sarek-docker-link]: https://hub.docker.com/r/maxulysse/sarek
[scilifelab-link]: https://www.scilifelab.se/
[snpeff-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeff.svg
[snpeff-docker-link]: https://hub.docker.com/r/maxulysse/snpeff
[snpeff-link]: http://snpeff.sourceforge.net/
[snpeffgrch37-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeffgrch37.svg
[snpeffgrch37-docker-link]: https://hub.docker.com/r/maxulysse/snpeffgrch37
[snpeffgrch38-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeffgrch38.svg
[snpeffgrch38-docker-link]: https://hub.docker.com/r/maxulysse/snpeffgrch38
[strelka-link]: https://github.com/Illumina/strelka
[vcftools-link]: https://vcftools.github.io/index.html
[vep-link]: https://github.com/Ensembl/ensembl-vep
[vepgrch37-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch37.svg
[vepgrch37-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch37
[vepgrch38-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch38.svg
[vepgrch38-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch38
