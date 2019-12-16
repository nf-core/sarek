# [![Sarek](docs/images/nf-core_sarek_logo.png "Sarek")](https://nf-co.re/sarek)

> **An open-source analysis pipeline to detect germline or somatic variants from whole genome or targeted sequencing**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![nf-core](https://img.shields.io/badge/nf--core-pipeline-brightgreen.svg)](https://nf-co.re/)
[![DOI](https://zenodo.org/badge/184289291.svg)](https://zenodo.org/badge/latestdoi/184289291)

[![GitHub Actions CI status](https://github.com/nf-core/sarek/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/sarek/actions?query=workflow%3A%22sarek+CI%22)
[![GitHub Actions Linting status](https://github.com/nf-core/sarek/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/sarek/actions?query=workflow%3A%22sarek+linting%22)
[![CircleCi build status](https://img.shields.io/circleci/project/github/nf-core/sarek?logo=circleci)](https://circleci.com/gh/nf-core/sarek/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/sarek.svg)](https://hub.docker.com/r/nfcore/sarek/)
[![Install with Singularity](https://img.shields.io/badge/use%20with-singularity-purple.svg)](https://www.sylabs.io/docs/)

[![Join us on Slack](https://img.shields.io/badge/slack-nfcore/sarek-blue.svg)](https://nfcore.slack.com/channels/sarek)

## Introduction

Sarek is a workflow designed to run analyses on whole genome or targeted sequencing data from regular samples or tumour / normal pairs and could include additional relapses.

It's built using [Nextflow](https://www.nextflow.io),
a domain specific language for workflow building,
across multiple compute infrastructures in a very portable manner.
Software dependencies are handled using [Conda](https://conda.io/), [Docker](https://www.docker.com) or [Singularity](https://www.sylabs.io/singularity/) - environment/container technologies that provide excellent reproducibility and ease of use.
Thus making installation trivial and results highly reproducible.

<p align="center">
    <a href="./docs/images/sarek_workflow.png"><img title="Sarek Workflow" src="./docs/images/sarek_workflow.png" width=40%></a>
</p>

It's listed on the [Elixir - Tools and Data Services Registry](https://bio.tools/Sarek), [Dockstore](https://dockstore.org/workflows/github.com/nf-core/sarek) and [omicX - Bioinformatics tools](https://omictools.com/sarek-tool).

## Documentation

The nf-core/sarek pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Install on a secure cluster](docs/install_bianca.md)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
    * [Extra documentation on reference](docs/reference.md)
3. [Running the pipeline](docs/usage.md)
    * [Examples](docs/use_cases.md)
    * [Input files documentation](docs/input.md)
    * [Documentation about containers](docs/containers.md)
4. [Output and how to interpret the results](docs/output.md)
    * [Complementary information about ASCAT](docs/ascat.md)
    * [Extra documentation on annotation](docs/annotation.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Credits

Sarek was developed at the [National Genomics Infastructure](https://ngisweden.scilifelab.se) and [National Bioinformatics Infastructure Sweden](https://nbis.se) which are both platforms at [SciLifeLab](https://scilifelab.se), with the support of [The Swedish Childhood Tumor Biobank (Barntumörbanken)](ttps://ki.se/forskning/barntumorbanken-0).

Main authors:

* [Maxime Garcia](https://github.com/MaxUlysse)
* [Szilveszter Juhos](https://github.com/szilvajuhos)

Helpful contributors:

* [Johannes Alneberg](https://github.com/alneberg)
* [Phil Ewels](https://github.com/ewels)
* [Jesper Eisfeldt](https://github.com/J35P312)
* [Malin Larsson](https://github.com/malinlarsson)
* [Marcel Martin](https://github.com/marcelm)
* [Alexander Peltzer](https://github.com/apeltzer)
* [Nilesh Tawari](https://github.com/nilesh-tawari)
* [arontommi](https://github.com/arontommi)
* [bjornnystedt](https://github.com/bjornnystedt)
* [gulfshores](https://github.com/gulfshores)
* [KochTobi](https://github.com/KochTobi)
* [pallolason](https://github.com/pallolason)
* [Sebastian-D](https://github.com/Sebastian-D)
* [silviamorins](https://github.com/silviamorins)

## Contributions & Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/channels/sarek) (you can join with [this invite](https://nf-co.re/join/slack)) or contact us: maxime.garcia@scilifelab.se, szilveszter.juhos@scilifelab.se

## CHANGELOG

* [CHANGELOG](CHANGELOG.md)

## Aknowledgements

[![Barntumörbanken](docs/images/BTB_logo.png)](https://ki.se/forskning/barntumorbanken) | [![SciLifeLab](docs/images/SciLifeLab_logo.png)](https://scilifelab.se)
:-:|:-:
[![National Genomics Infrastructure](docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/) | [![National Bioinformatics Infrastructure Sweden](docs/images/NBIS_logo.png)](https://nbis.se)

## Citation

If you use `nf-core/sarek` for your analysis, please cite the `Sarek` pre-print as follows:
> Garcia MU, Juhos S, Larsson M, Olason PI, Martin M, Eisfeldt J, DiLorenzo S, Sandgren J, de Ståhl TD, Wirta V, Nistér M, Nystedt B, Käller M. **Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants**. *bioRxiv*. 2018. p. 316976. [doi: 10.1101/316976](https://www.biorxiv.org/content/10.1101/316976v1).

You can cite the sarek zenodo record for a specific version using the following [doi: 10.5281/zenodo.3476426](https://zenodo.org/badge/latestdoi/184289291)

You can cite the `nf-core` pre-print as follows:
> Ewels PA, Peltzer A, Fillinger S, Alneberg JA, Patel H, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. **nf-core: Community curated bioinformatics pipelines**. *bioRxiv*. 2019. p. 610741. [doi: 10.1101/610741](https://www.biorxiv.org/content/10.1101/610741v3).
