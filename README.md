# [![nf-core/sarek](docs/images/nf-core-sarek_logo.png "nf-core/sarek")](https://nf-co.re/sarek)

> **An open-source analysis pipeline to detect germline or somatic variants from whole genome or targeted sequencing**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/nf--core-pipeline-brightgreen.svg)](https://nf-co.re/)
[![DOI](https://zenodo.org/badge/184289291.svg)](https://zenodo.org/badge/latestdoi/184289291)

[![GitHub Actions CI status](https://github.com/nf-core/sarek/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/sarek/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting status](https://github.com/nf-core/sarek/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/sarek/actions?query=workflow%3A%22nf-core+linting%22)
[![CircleCi build status](https://img.shields.io/circleci/project/github/nf-core/sarek?logo=circleci)](https://circleci.com/gh/nf-core/sarek/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/sarek.svg)](https://hub.docker.com/r/nfcore/sarek)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23sarek-4A154B?logo=slack)](https://nfcore.slack.com/channels/sarek)

## Introduction

Sarek is a workflow designed to detect variants on whole genome or targeted sequencing data.
Initially designed for Human, and Mouse, it can work on any species with a reference genome.
Sarek can also handle tumour / normal pairs and could include additional relapses.

The pipeline is built using [`Nextflow`](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with `Docker` containers making installation trivial and results highly reproducible.

<p align="center">
    <img title="Sarek Workflow" src="docs/images/sarek_workflow.png" width=40%>
</p>

It's listed on [Elixir - Tools and Data Services Registry](https://bio.tools/nf-core-sarek) and [Dockstore](https://dockstore.org/workflows/github.com/nf-core/sarek).

## Quick Start

1. Install [`Nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/sarek -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute.
    > If so, you can simply use `-profile <institute>` in your command.
    > This will enable either `Docker` or `Singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    ```bash
    nextflow run nf-core/sarek -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input '*.tsv' --genome GRCh38
    ```

See [usage docs](https://nf-co.re/sarek/usage) for all of the available options when running the pipeline.

## Pipeline Summary

By default, the pipeline currently performs the following:

* Sequencing quality control (`FastQC`)
* Map Reads to Reference (`BWA mem`)
* Mark Duplicates (`GATK MarkDuplicatesSpark`)
* Base (Quality Score) Recalibration (`GATK BaseRecalibrator`, `GATK ApplyBQSR`)
* Preprocessing quality control (`samtools stats`)
* Preprocessing quality control (`Qualimap bamqc`)
* Overall pipeline run summaries (`MultiQC`)

## Documentation

The nf-core/sarek pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/sarek/usage) and [output](https://nf-co.re/sarek/output).

## Credits

Sarek was developed at the [National Genomics Infastructure](https://ngisweden.scilifelab.se) and [National Bioinformatics Infastructure Sweden](https://nbis.se) which are both platforms at [SciLifeLab](https://scilifelab.se), with the support of [The Swedish Childhood Tumor Biobank (Barntumörbanken)](https://ki.se/forskning/barntumorbanken).
[QBiC](https://www.qbic.uni-tuebingen.de/) later joined and helped with further development.

Main authors:

* [Gisela Gabernet](https://github.com/ggabernet)
* [Maxime Garcia](https://github.com/maxulysse)
* [Friederike Hanssen](https://github.com/FriederikeHanssen)
* [Szilveszter Juhos](https://github.com/szilvajuhos)

Helpful contributors:

* [Adrian Lärkeryd](https://github.com/adrlar)
* [Alexander Peltzer](https://github.com/apeltzer)
* [Chela James](https://github.com/chelauk)
* [David Mas-Ponte](https://github.com/davidmasp)
* [Francesco L](https://github.com/nibscles)
* [Harshil Patel](https://github.com/drpatelh)
* [James A. Fellows Yates](https://github.com/jfy133)
* [Jesper Eisfeldt](https://github.com/J35P312)
* [Johannes Alneberg](https://github.com/alneberg)
* [José Fernández Navarro](https://github.com/jfnavarro)
* [Lucia Conde](https://github.com/lconde-ucl)
* [Malin Larsson](https://github.com/malinlarsson)
* [Marcel Martin](https://github.com/marcelm)
* [Nilesh Tawari](https://github.com/nilesh-tawari)
* [Olga Botvinnik](https://github.com/olgabot)
* [Paul Cantalupo](https://github.com/pcantalupo)
* [Phil Ewels](https://github.com/ewels)
* [Sabrina Krakau](https://github.com/skrakau)
* [Sebastian-D](https://github.com/Sebastian-D)
* [Tobias Koch](https://github.com/KochTobi)
* [Winni Kretzschmar](https://github.com/winni2k)
* [arontommi](https://github.com/arontommi)
* [bjornnystedt](https://github.com/bjornnystedt)
* [cgpu](https://github.com/cgpu)
* [gulfshores](https://github.com/gulfshores)
* [pallolason](https://github.com/pallolason)
* [silviamorins](https://github.com/silviamorins)

## Contributions & Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#sarek` channel](https://nfcore.slack.com/channels/sarek) (you can join with [this invite](https://nf-co.re/join/slack)), or contact us: [Maxime Garcia](mailto:maxime.garcia@scilifelab.se?subject=[GitHub]%20nf-core/sarek), [Szilvester Juhos](mailto:szilveszter.juhos@scilifelab.se?subject=[GitHub]%20nf-core/sarek)

## CHANGELOG

* [CHANGELOG](CHANGELOG.md)

## Acknowledgements

[![Barntumörbanken](docs/images/BTB_logo.png)](https://ki.se/forskning/barntumorbanken) | [![SciLifeLab](docs/images/SciLifeLab_logo.png)](https://scilifelab.se)
:-:|:-:
[![National Genomics Infrastructure](docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/) | [![National Bioinformatics Infrastructure Sweden](docs/images/NBIS_logo.png)](https://nbis.se)
[![QBiC](docs/images/QBiC_logo.png)](hhttps://www.qbic.uni-tuebingen.de) |

## Citations

If you use `nf-core/sarek` for your analysis, please cite the `Sarek` article as follows:
> Garcia M, Juhos S, Larsson M et al. **Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants [version 2; peer review: 2 approved]** *F1000Research* 2020, 9:63 [doi: 10.12688/f1000research.16665.2](http://dx.doi.org/10.12688/f1000research.16665.2).

You can cite the sarek zenodo record for a specific version using the following [doi: 10.5281/zenodo.3476426](https://zenodo.org/badge/latestdoi/184289291)

In addition, references of tools and data used in this pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
