# Containers

A container named after the process is made for each process. If a container can be reused, it will be named after the tool used.

## bcftools [![bcftools-docker status][bcftools-docker-badge]][bcftools-docker-link]

- Based on `debian:8.6`
- Contain **[BCFTools][bcftools-link]** 1.4

## concatvcf [![concatvcf-docker status][concatvcf-docker-badge]][concatvcf-docker-link]

- Based on `debian:8.6`
- Contain **[pigz][pigz-link]** 2.3.4

## fastqc [![fastqc-docker status][fastqc-docker-badge]][fastqc-docker-link]

- Based on `openjdk:8`
- Contain **[FastQC][fastqc-link]** 0.11.5

## freebayes [![freebayes-docker status][freebayes-docker-badge]][freebayes-docker-link]

- Based on `debian:8.6`
- Contain **[FreeBayes][freebayes-link]** 1.1.0

## gatk [![gatk-docker status][gatk-docker-badge]][gatk-docker-link]

- Based on `openjdk:8-slim`
- Contain **[GATK][gatk-link]** 3.7

## htslib [![htslib-docker status][htslib-docker-badge]][htslib-docker-link]

- Based on `debian:8.6`
- Contain **[HTSlib][htslib-link]** 1.4

## igvtools [![igvtools-docker status][igvtools-docker-badge]][igvtools-docker-link]

- Based on `openjdk:8-slim`
- Contain **[IGVTools][igvtools-link]** 2.3.91

## mapreads [![mapreads-docker status][mapreads-docker-badge]][mapreads-docker-link]

- Based on `maxulysse/samtools:1.1`
- Contain **[BWA][bwa-link]** 0.7.8

## multiqc [![multiqc-docker status][multiqc-docker-badge]][multiqc-docker-link]

- Based on `openjdk:8-slim`
- Contain **[MultiQC][multiqc-link]** 1.0

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

- Based on `maxulysse/samtools:1.1`
- Contain **[AlleleCount][allelecount-link]** 2.2.0

## runascat [![runascat-docker status][runascat-docker-badge]][runascat-docker-link]

- Based on `maxulysse/samtools:1.1`
- Contain **[RColorBrewer][rcolorbrewer-link]**

## runconvertallelecounts [![runconvertallelecounts-docker status][runconvertallelecounts-docker-badge]][runconvertallelecounts-docker-link]

- Based on `r-base:3.3.2`

## runmanta [![runmanta-docker status][runmanta-docker-badge]][runmanta-docker-link]

- Based on `maxulysse/samtools:1.1`
- Contain **[Manta][manta-link]** 1.0.3

## samtools [![samtools-docker status][samtools-docker-badge]][samtools-docker-link]

- Based on `debian:8.6`
- Contain **[samtools][samtools-link]** 1.4

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

## strelka [![strelka-docker status][strelka-docker-badge]][strelka-docker-link]

- Based on `debian:8.6`
- Contain **[Strelka][strelka-link]** 1.0.15

## vep [![vep-docker status][vep-docker-badge]][vep-docker-link]

- Based on `ubuntu:16.04`
- Contain **[VEP][vep-link]** 90.1

## vepgrch37 [![vepgrch37-docker status][vepgrch37-docker-badge]][vepgrch37-docker-link]

- Based on `maxulysse/vep`
- Contain **[VEP][vep-link]** 90.1
- Contain GRCh37

## vepgrch38 [![vepgrch38-docker status][vepgrch38-docker-badge]][vepgrch38-docker-link]

- Based on `maxulysse/vep`
- Contain **[VEP][vep-link]** 90.1
- Contain GRCh38

---
[![](images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](images/NGI_logo.png "NGI")][ngi-link]
[![](images/NBIS_logo.png "NBIS")][nbis-link]

[allelecount-link]: https://github.com/cancerit/alleleCount
[bcftools-docker-badge]: https://img.shields.io/docker/automated/maxulysse/bcftools.svg
[bcftools-docker-link]: https://hub.docker.com/r/maxulysse/bcftools
[bcftools-link]: https://github.com/samtools/bcftools
[bwa-link]: https://github.com/lh3/bwa
[concatvcf-docker-badge]: https://img.shields.io/docker/automated/maxulysse/concatvcf.svg
[concatvcf-docker-link]: https://hub.docker.com/r/maxulysse/concatvcf
[fastqc-docker-badge]: https://img.shields.io/docker/automated/maxulysse/fastqc.svg
[fastqc-docker-link]: https://hub.docker.com/r/maxulysse/fastqc
[fastqc-link]: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[freebayes-docker-badge]: https://img.shields.io/docker/automated/maxulysse/freebayes.svg
[freebayes-docker-link]: https://hub.docker.com/r/maxulysse/freebayes
[freebayes-link]: https://github.com/ekg/freebayes
[gatk-docker-badge]: https://img.shields.io/docker/automated/maxulysse/gatk.svg
[gatk-docker-link]: https://hub.docker.com/r/maxulysse/gatk
[gatk-link]: https://github.com/broadgsa/gatk-protected
[htslib-docker-badge]: https://img.shields.io/docker/automated/maxulysse/htslib.svg
[htslib-docker-link]: https://hub.docker.com/r/maxulysse/htslib
[htslib-link]: https://github.com/samtools/htslib
[igvtools-docker-badge]: https://img.shields.io/docker/automated/maxulysse/igvtools.svg
[igvtools-docker-link]: https://hub.docker.com/r/maxulysse/igvtools
[igvtools-link]: http://software.broadinstitute.org/software/igv/
[manta-link]: https://github.com/Illumina/manta
[mapreads-docker-badge]: https://img.shields.io/docker/automated/maxulysse/mapreads.svg
[mapreads-docker-link]: https://hub.docker.com/r/maxulysse/mapreads
[multiqc-docker-badge]: https://img.shields.io/docker/automated/maxulysse/multiqc.svg
[multiqc-docker-link]: https://hub.docker.com/r/maxulysse/multiqc
[multiqc-link]: https://github.com/ewels/MultiQC/
[mutect1-docker-badge]: https://img.shields.io/docker/automated/maxulysse/mutect1.svg
[mutect1-docker-link]: https://hub.docker.com/r/maxulysse/mutect1
[mutect1-link]: https://github.com/broadinstitute/mutect
[nbis-link]: https://www.nbis.se/
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A50.22.2-brightgreen.svg
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[picard-docker-badge]: https://img.shields.io/docker/automated/maxulysse/picard.svg
[picard-docker-link]: https://hub.docker.com/r/maxulysse/picard
[picard-link]: https://github.com/broadinstitute/picard
[pigz-link]: https://zlib.net/pigz/
[qualimap-docker-badge]: https://img.shields.io/docker/automated/maxulysse/qualimap.svg
[qualimap-docker-link]: https://hub.docker.com/r/maxulysse/qualimap
[qualimap-link]: http://qualimap.bioinfo.cipf.es
[rcolorbrewer-link]: https://CRAN.R-project.org/package=RColorBrewer
[runallelecount-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runallelecount.svg
[runallelecount-docker-link]: https://hub.docker.com/r/maxulysse/runallelecount
[runascat-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runascat.svg
[runascat-docker-link]: https://hub.docker.com/r/maxulysse/runascat
[runconvertallelecounts-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runconvertallelecounts.svg
[runconvertallelecounts-docker-link]: https://hub.docker.com/r/maxulysse/runconvertallelecounts
[runmanta-docker-badge]: https://img.shields.io/docker/automated/maxulysse/runmanta.svg
[runmanta-docker-link]: https://hub.docker.com/r/maxulysse/runmanta
[samtools-docker-badge]: https://img.shields.io/docker/automated/maxulysse/samtools.svg
[samtools-docker-link]: https://hub.docker.com/r/maxulysse/samtools
[samtools-link]: https://github.com/samtools/samtools
[scilifelab-link]: https://www.scilifelab.se/
[scilifelab-stockholm-link]: https://www.scilifelab.se/platforms/ngi/
[snpeff-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeff.svg
[snpeff-docker-link]: https://hub.docker.com/r/maxulysse/snpeff
[snpeff-link]: http://snpeff.sourceforge.net/
[snpeffgrch37-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeffgrch37.svg
[snpeffgrch37-docker-link]: https://hub.docker.com/r/maxulysse/snpeffgrch37
[snpeffgrch38-docker-badge]: https://img.shields.io/docker/automated/maxulysse/snpeffgrch38.svg
[snpeffgrch38-docker-link]: https://hub.docker.com/r/maxulysse/snpeffgrch38
[strelka-docker-badge]: https://img.shields.io/docker/automated/maxulysse/strelka.svg
[strelka-docker-link]: https://hub.docker.com/r/maxulysse/strelka
[strelka-link]: https://github.com/Illumina/strelka
[vep-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vep.svg
[vep-docker-link]: https://hub.docker.com/r/maxulysse/vep
[vep-link]: https://github.com/Ensembl/ensembl-vep
[vepgrch37-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch37.svg
[vepgrch37-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch37
[vepgrch38-docker-badge]: https://img.shields.io/docker/automated/maxulysse/vepgrch38.svg
[vepgrch38-docker-link]: https://hub.docker.com/r/maxulysse/vepgrch38
[version-badge]: https://img.shields.io/github/release/MaxUlysse/CAW-containers.svg
[version-link]: https://github.com/MaxUlysse/CAW-containers/releases/latest
