# [![Sarek](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/Sarek_logo.png "Sarek")](http://opensource.scilifelab.se/projects/sarek/)

#### An open-source analysis pipeline to detect germline or somatic variants from whole genome sequencing

[![Nextflow version][nextflow-badge]][nextflow-link]
[![Travis build status][travis-badge]][travis-link]
[![Join the chat on https://gitter.im/SciLifeLab/Sarek][gitter-badge]][gitter-link]

[![MIT License][license-badge]][license-link]
[![Sarek version][version-badge]][version-link]
[![DOI][zenodo-badge]][zenodo-link]

[![Install with bioconda][bioconda-badge]][bioconda-link]
[![Docker Container available][docker-badge]][docker-link]

## Introduction

<img align="right" title="CAW" src="https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/CAW_logo.png">

Previously known as the Cancer Analysis Workflow (CAW),
Sarek is a workflow designed to run analyses on WGS data from regular samples or tumour / normal pairs, including relapse samples if required.

It's built using [Nextflow][nextflow-link], a domain specific language for workflow building.
Software dependencies are handled using [Docker](https://www.docker.com) or [Singularity](http://singularity.lbl.gov) - container technologies that provide excellent reproducibility and ease of use.
Singularity has been designed specifically for high-performance computing environments.
This means that although Sarek has been primarily designed for use with the Swedish [UPPMAX HPC systems](https://www.uppmax.uu.se), it should be able to run on any system that supports these two tools.

Sarek was developed at the [National Genomics Infastructure][ngi-link] and [National Bioinformatics Infastructure Sweden][nbis-link] which are both platforms at [SciLifeLab][scilifelab-link].
It is listed on the [Elixir - Tools and Data Services Registry](https://bio.tools/Sarek).

## Workflow steps

Sarek is built with several workflow scripts.
A wrapper script contained within the repository makes it easy to run the different workflow scripts as a single job.
To test your installation, follow the [tests documentation.](https://github.com/SciLifeLab/Sarek/blob/master/docs/TESTS.md)

Raw FastQ files or aligned BAM files (with or without realignment & recalibration) can be used as inputs.
You can choose which variant callers to use, plus the pipeline is capable of accommodating additional variant calling software or CNV callers if required.

The worflow steps and tools used are as follows:

1. **Preprocessing** - `main.nf` _(based on [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/))_
    * Map reads to Reference
        * [BWA](http://bio-bwa.sourceforge.net/)
    * Mark Duplicates
        * [GATK MarkDuplicates](https://github.com/broadinstitute/gatk)
    * Base (Quality Score) Recalibration
        * [GATK BaseRecalibrator](https://github.com/broadinstitute/gatk)
        * [GATK ApplyBQSR](https://github.com/broadinstitute/gatk)
2. **Germline variant calling** - `germlineVC.nf`
    * SNVs and small indels
        * [GATK HaplotyeCaller](https://github.com/broadinstitute/gatk)
        * [Strelka2](https://github.com/Illumina/strelka)
    * Structural variants
        * [Manta](https://github.com/Illumina/manta)
3. **Somatic variant calling** - `somaticVC.nf` _(optional)_
    * SNVs and small indels
        * [MuTect2](https://github.com/broadinstitute/gatk)
        * [Freebayes](https://github.com/ekg/freebayes)
        * [Strelka2](https://github.com/Illumina/strelka)
    * Structural variants
        * [Manta](https://github.com/Illumina/manta)
    * Sample heterogeneity, ploidy and CNVs
        * [ASCAT](https://github.com/Crick-CancerGenomics/ascat)
4. **Annotation** - `annotate.nf` _(optional)_
    * Variant annotation
        * [SnpEff](http://snpeff.sourceforge.net/)
        * [VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html)
5. **Reporting** - `runMultiQC.nf`
    * Reporting
        * [MultiQC](http://multiqc.info)

## Documentation

The Sarek pipeline comes with documentation in the `docs/` directory:

01. [Installation documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/INSTALL.md)
02. [Installation documentation specific for UPPMAX `rackham`](https://github.com/SciLifeLab/Sarek/blob/master/docs/INSTALL_RACKHAM.md)
03. [Installation documentation specific for UPPMAX `bianca`](https://github.com/SciLifeLab/Sarek/blob/master/docs/INSTALL_BIANCA.md)
04. [Tests documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/TESTS.md)
05. [Reference files documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/REFERENCES.md)
06. [Configuration and profiles documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/CONFIG.md)
07. [Intervals documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/INTERVALS.md)
08. [Running the pipeline](https://github.com/SciLifeLab/Sarek/blob/master/docs/USAGE.md)
09. [Examples](https://github.com/SciLifeLab/Sarek/blob/master/docs/USE_CASES.md)
10. [TSV file documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/TSV.md)
11. [Processes documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/PROCESS.md)
12. [Documentation about containers](https://github.com/SciLifeLab/Sarek/blob/master/docs/CONTAINERS.md)
13. [Documentation about building](https://github.com/SciLifeLab/Sarek/blob/master/docs/BUILD.md)
14. [More information about ASCAT](https://github.com/SciLifeLab/Sarek/blob/master/docs/ASCAT.md)
15. [Folder structure](https://github.com/SciLifeLab/Sarek/blob/master/docs/FOLDER.md)

## Contributions & Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](https://github.com/SciLifeLab/Sarek/blob/master/.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Gitter][gitter-link] or contact us: maxime.garcia@scilifelab.se, szilveszter.juhos@scilifelab.se

## CHANGELOG

- [CHANGELOG](https://github.com/SciLifeLab/Sarek/blob/master/CHANGELOG.md)

## Credits

Main authors:
* [Maxime Garcia](https://github.com/MaxUlysse)
* [Szilveszter Juhos](https://github.com/szilvajuhos)

Helpful contributors:
* [Sebastian DiLorenzo](https://github.com/Sebastian-D)
* [Jesper Eisfeldt](https://github.com/J35P312)
* [Phil Ewels](https://github.com/ewels)
* [Max Käller](https://github.com/gulfshores)
* [Malin Larsson](https://github.com/malinlarsson)
* [Marcel Martin](https://github.com/marcelm)
* [Björn Nystedt](https://github.com/bjornnystedt)
* [Pall Olason](https://github.com/pallolason)

--------------------------------------------------------------------------------

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![NGI](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/NGI_logo.png "NGI")][ngi-link]
[![NBIS](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/NBIS_logo.png "NBIS")][nbis-link]

[bioconda-badge]:https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAFAAAABQCAMAAAC5zwKfAAAC7lBMVEUAAAD/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////jtXoAAAA+XRSTlMAAQIDBAUGBwgJCgsMDQ4PEBESExQVFhcYGRobHB0eICEiIyQlJicoKissLS4wMTIzNDU2Nzg5Ojs8PT4/QEFCQ0RFRkdISUpLTE1OT1BRUlNUVVZXWFlaW1xdXl9gYWJj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
[bioconda-link]:http://bioconda.github.io/
[docker-badge]: https://img.shields.io/docker/automated/maxulysse/sarek.svg?logo=docker
[docker-link]: https://hub.docker.com/r/maxulysse/sarek
[gitter-badge]: https://img.shields.io/gitter/room/SciLifeLab/Sarek.svg?logo=gitter&logoColor=white&colorB=4fb99a
[gitter-link]: https://gitter.im/SciLifeLab/Sarek
[license-badge]: https://img.shields.io/github/license/SciLifeLab/Sarek.svg
[license-link]: https://github.com/SciLifeLab/Sarek/blob/master/LICENSE
[nbis-link]: https://www.nbis.se/
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A50.31.0-brightgreen.svg?logo=data:image/png;base64,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
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
[travis-badge]: https://img.shields.io/travis/SciLifeLab/Sarek.svg?logo=travis
[travis-link]: https://travis-ci.org/SciLifeLab/Sarek
[version-badge]: https://img.shields.io/github/release/SciLifeLab/Sarek.svg?logo=github&logoColor=white
[version-link]: https://github.com/SciLifeLab/Sarek/releases/latest
[zenodo-badge]: https://zenodo.org/badge/54024046.svg
[zenodo-link]: https://zenodo.org/badge/latestdoi/54024046
