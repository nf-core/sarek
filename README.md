# CAW
[![caw version][version-badge]][version-link] [![mit license][license-badge]][license-link] [![nextflow version][nextflow-badge]][nextflow-link]

Nextflow Cancer Analysis Workflow developed at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

## Authors
- Sebastian DiLorenzo (@Sebastian-D)
- Jesper Eisfeldt (@J35P312)
- Maxime Garcia (@MaxUlysse)
- Szilveszter Juhos (@szilvajuhos)
- Max Käller (@gulfshores)
- Malin Larsson (@malinlarsson)
- Björn Nystedt (@bjornnystedt)
- Pall Olason (@pallolason)
- Pelin Sahlén (@pelinakan)

## Installation and first execution
To use this pipeline, you need to have a working version of Nextflow installed.
- See the [Install Nextflow documentation](https://github.com/SciLifeLab/NGI-NextflowDocs/blob/master/docs/INSTALL.md)
- See the [first execution documentation](doc/FIRST_RUN.md)

## Usage
See the [usage documentation](doc/USAGE.md)

## Workflow processes
See the [processes documentation](doc/PROCESS.md)

## TSV files
See the [workflow TSV file documentation](doc/TSV.md)

## Use cases
See the [workflow use cases documentation](doc/USE_CASES.md)

## Tools and dependencies
- alleleCount 2.2.0
- bwa 0.7.8
- FastQC 0.11.5
- freebayes
- GATK 3.7
- gcc 4.9.2
- java jdk 7
- java jdk 8
- manta 1.0.0
- MultiQC 0.9
- MuTecT 1.1.5
- nextflow >= 0.22.1
- perl 5.18.4
- picard 2.0.1
- R 3.2.3
- samtools 1.3
- strelka 1.0.15

[license-badge]: https://img.shields.io/badge/license-MIT-blue.svg
[license-link]: https://github.com/SciLifeLab/CAW/blob/master/LICENSE
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A50.22.2-brightgreen.svg
[nextflow-link]: https://www.nextflow.io/
[version-badge]: https://img.shields.io/badge/version-1.0-green.svg
[version-link]: https://github.com/SciLifeLab/CAW/releases/tag/1.0
