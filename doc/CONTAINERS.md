# Containers

Subsets of all containers can be dowloaded:

For processing + germline variant calling + Reports:
 - [fastqc](#fastqc-)
 - [gatk](#gatk-)
 - [multiqc](#multiqc-)
 - [picard](#picard-)
 - [qualimap](#qualimap-)
 - [sarek](#sarek-)

For processing + somatic variant calling + Reports:
 - [fastqc](#fastqc-)
 - [freebayes](#freebayes-)
 - [gatk](#gatk-)
 - [multiqc](#multiqc-)
 - [mutect1](#mutect1-)
 - [picard](#picard-)
 - [qualimap](#qualimap-)
 - [r-base](#r-base-)
 - [runallelecount](#runallelecount-)
 - [sarek](#sarek-)

For annotation for GRCh37, you will need:
 - [snpeffgrch37](#snpeffgrch37-)
 - [vepgrch37](#vepgrch37-)

For annotation for GRCh38, you will need:
 - [snpeffgrch38](#snpeffgrch38-)
 - [vepgrch38](#vepgrch38-)

A container named after the process is made for each process. If a container can be reused, it will be named after the tool used.

## fastqc [![fastqc-docker status][fastqc-docker-badge]][fastqc-docker-link]

- Based on `openjdk:8`
- Contain **[FastQC][fastqc-link]** 0.11.5

## freebayes [![freebayes-docker status][freebayes-docker-badge]][freebayes-docker-link]

- Based on `debian:8.6`
- Contain **[FreeBayes][freebayes-link]** 1.1.0

## gatk [![gatk-docker status][gatk-docker-badge]][gatk-docker-link]

- Based on `broadinstitute/gatk3:3.8-0`
- Contain **[GATK][gatk-link]** 3.8

## igvtools [![igvtools-docker status][igvtools-docker-badge]][igvtools-docker-link]

- Based on `openjdk:8-slim`
- Contain **[IGVTools][igvtools-link]** 2.3.98

## multiqc [![multiqc-docker status][multiqc-docker-badge]][multiqc-docker-link]

- Based on `ewels/multiqc:v1.4`
- Contain **[MultiQC][multiqc-link]** 143

## mutect1 [![mutect1-docker status][mutect1-docker-badge]][mutect1-docker-link]

- Based on `openjdk:7-slim`
- Contain **[MuTect1][mutect1-link]** 1.5

## picard [![picard-docker status][picard-docker-badge]][picard-docker-link]

- Based on `openjdk:8-slim`
- Contain **[Picard][picard-link]** 2.0.1

## qualimap [![qualimap-docker status][qualimap-docker-badge]][qualimap-docker-link]

- Based on `openjdk:8`
- Contain **[qualimap][qualimap-link]** 2.2.1

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
[fastqc-docker-badge]: https://img.shields.io/docker/automated/maxulysse/fastqc.svg
[fastqc-docker-link]: https://hub.docker.com/r/maxulysse/fastqc
[fastqc-link]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[freebayes-docker-badge]: https://img.shields.io/docker/automated/maxulysse/freebayes.svg
[freebayes-docker-link]: https://hub.docker.com/r/maxulysse/freebayes
[freebayes-link]: https://github.com/ekg/freebayes
[gatk-docker-badge]: https://img.shields.io/docker/automated/maxulysse/gatk.svg
[gatk-docker-link]: https://hub.docker.com/r/maxulysse/gatk
[gatk-link]: https://github.com/broadgsa/gatk-protected
[htslib-link]: https://github.com/samtools/htslib
[igvtools-docker-badge]: https://img.shields.io/docker/automated/maxulysse/igvtools.svg
[igvtools-docker-link]: https://hub.docker.com/r/maxulysse/igvtools
[igvtools-link]: http://software.broadinstitute.org/software/igv/
[manta-link]: https://github.com/Illumina/manta
[multiqc-docker-badge]: https://img.shields.io/docker/automated/maxulysse/multiqc.svg
[multiqc-docker-link]: https://hub.docker.com/r/maxulysse/multiqc
[multiqc-link]: https://github.com/ewels/MultiQC/
[mutect1-docker-badge]: https://img.shields.io/docker/automated/maxulysse/mutect1.svg
[mutect1-docker-link]: https://hub.docker.com/r/maxulysse/mutect1
[mutect1-link]: https://github.com/broadinstitute/mutect
[nbis-link]: https://www.nbis.se/
[ngi-link]: https://ngisweden.scilifelab.se/
[picard-docker-badge]: https://img.shields.io/docker/automated/maxulysse/picard.svg
[picard-docker-link]: https://hub.docker.com/r/maxulysse/picard
[picard-link]: https://github.com/broadinstitute/picard
[qualimap-docker-badge]: https://img.shields.io/docker/automated/maxulysse/qualimap.svg
[qualimap-docker-link]: https://hub.docker.com/r/maxulysse/qualimap
[qualimap-link]: http://qualimap.bioinfo.cipf.es
[rcolorbrewer-link]: https://CRAN.R-project.org/package=RColorBrewer
[runallelecount-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runallelecount.svg
[runallelecount-docker-link]: https://hub.docker.com/r/maxulysse/runallelecount
[runascat-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runascat.svg
[runascat-docker-link]: https://hub.docker.com/r/maxulysse/runascat
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
[vep-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vep.svg
[vep-docker-link]: https://hub.docker.com/r/maxulysse/vep
[vep-link]: https://github.com/Ensembl/ensembl-vep
[vepgrch37-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch37.svg
[vepgrch37-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch37
[vepgrch38-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch38.svg
[vepgrch38-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch38
