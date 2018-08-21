# [![Sarek](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/Sarek_logo.png "Sarek")](http://opensource.scilifelab.se/projects/sarek/)

#### An open-source analysis pipeline to detect germline or somatic variants from whole genome sequencing

[![Nextflow version][nextflow-badge]][nextflow-link]
[![Travis build status][travis-badge]][travis-link]
[![Join the chat at [gitter](gitter-link)][gitter-badge]][gitter-link]

[![MIT License][license-badge]][license-link]
[![Sarek version][version-badge]][version-link]
[![DOI][zenodo-badge]][zenodo-link]

[![Install with bioconda][bioconda-badge]][bioconda-link]
[![Docker Container available][docker-badge]][docker-link]

## Introduction

<img align="right" title="CAW" src="https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/CAW_logo.png">

Previously known as the Cancer Analysis Workflow (CAW),
Sarek is a workflow designed to run analyses on WGS data from regular samples or tumour / normal pairs, including relapse samples if required.

It's built using [Nextflow][nextflow-link], a bioinformatics domain specific language for workflow building. Software dependencies are handled using [Docker](https://www.docker.com) or [Singularity](http://singularity.lbl.gov) - container technologies that provide excellent reproducibility and ease of use.
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
    * Read alignment
        * [BWA](http://bio-bwa.sourceforge.net/)
    * Read realignment and recalibration of short-read data
        * [GATK](https://github.com/broadgsa/gatk-protected)
2. **Germline variant calling** - `germlineVC.nf`
    * SNVs and small indels
        * [GATK HaplotyeCaller](https://github.com/broadgsa/gatk-protected)
        * [Strelka](https://github.com/Illumina/strelka)
    * Structural variants
        * [Manta](https://github.com/Illumina/manta)
3. **Somatic variant calling** - `somaticVC.nf` _(optional)_
    * SNVs and small indels
        * [MuTect2](https://github.com/broadgsa/gatk-protected)
        * [Freebayes](https://github.com/ekg/freebayes)
        * [Strelka](https://github.com/Illumina/strelka)
    * Structural variants
        * [Manta](https://github.com/Illumina/manta)
    * Sample heterogeneity, ploidy and CNVs
        * [ASCAT](https://github.com/Crick-CancerGenomics/ascat)
4. **Annotation** - `annotate.nf` _(optional)_
    * Variant annotation
        * [SnpEff](http://snpeff.sourceforge.net/)
        * [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) (Variant Effect Predictor)
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

[bioconda-badge]:https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg
[bioconda-link]:http://bioconda.github.io/
[docker-badge]: https://img.shields.io/docker/automated/maxulysse/sarek.svg
[docker-link]: https://hub.docker.com/r/maxulysse/sarek
[gitter-badge]: https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg
[gitter-link]: https://gitter.im/SciLifeLab/Sarek
[license-badge]: https://img.shields.io/github/license/SciLifeLab/Sarek.svg
[license-link]: https://github.com/SciLifeLab/Sarek/blob/master/LICENSE
[nbis-link]: https://www.nbis.se/
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A50.31.0-brightgreen.svg
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
[travis-badge]: https://api.travis-ci.org/SciLifeLab/Sarek.svg
[travis-link]: https://travis-ci.org/SciLifeLab/Sarek
[version-badge]: https://img.shields.io/github/release/SciLifeLab/Sarek.svg
[version-link]: https://github.com/SciLifeLab/Sarek/releases/latest
[zenodo-badge]: https://zenodo.org/badge/54024046.svg
[zenodo-link]: https://zenodo.org/badge/latestdoi/54024046
