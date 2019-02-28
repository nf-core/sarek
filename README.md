# [![Sarek](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/Sarek_logo.png "Sarek")](http://sarek.scilifelab.se/)

#### An open-source analysis pipeline to detect germline or somatic variants from whole genome or targeted sequencing

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
Software dependencies are handled using [Docker](https://www.docker.com) or [Singularity](https://www.sylabs.io/singularity/) - container technologies that provide excellent reproducibility and ease of use.
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
09. [Running the pipeline using Conda](https://github.com/SciLifeLab/Sarek/blob/master/docs/CONDA.md)
10. [Command line parameters](https://github.com/SciLifeLab/Sarek/blob/master/docs/PARAMETERS.md)
11. [Examples](https://github.com/SciLifeLab/Sarek/blob/master/docs/USE_CASES.md)
12. [Input files documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/INPUT.md)
13. [Processes documentation](https://github.com/SciLifeLab/Sarek/blob/master/docs/PROCESS.md)
14. [Documentation about containers](https://github.com/SciLifeLab/Sarek/blob/master/docs/CONTAINERS.md)
15. [Complementary information about ASCAT](https://github.com/SciLifeLab/Sarek/blob/master/docs/ASCAT.md)
16. [Complementary information about annotations](https://github.com/SciLifeLab/Sarek/blob/master/docs/ANNOTATION.md)
17. [Output documentation structure](https://github.com/SciLifeLab/Sarek/blob/master/docs/OUTPUT.md)

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
* [Johannes Alneberg](https://github.com/alneberg)
* [Sebastian DiLorenzo](https://github.com/Sebastian-D)
* [Jesper Eisfeldt](https://github.com/J35P312)
* [Phil Ewels](https://github.com/ewels)
* [Max Käller](https://github.com/gulfshores)
* [Malin Larsson](https://github.com/malinlarsson)
* [Marcel Martin](https://github.com/marcelm)
* [Björn Nystedt](https://github.com/bjornnystedt)
* [Pall Olason](https://github.com/pallolason)
* [Aron Skaftason](https://github.com/arontommi)
* [Nilesh Tawari](https://github.com/nilesh-tawari)

--------------------------------------------------------------------------------

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![NGI](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/NGI_logo.png "NGI")][ngi-link]
[![NBIS](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/docs/images/NBIS_logo.png "NBIS")][nbis-link]

[bioconda-badge]:https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADEAAAAyCAYAAAD1CDOyAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAN1wAADdcBQiibeAAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAAa2SURBVGiBxZprjFVXFcd/e2aA6UAoYGrk0aZYKvXdwWAyoDyswRqhxmpNjQFrNOIHTR+aJhoxrWBoAz4aGz80bdWCsW1qq5IGSlvDQA1aUGw7BEuR4dFCrVSY0qFYYH5+WHt674zzOHdm7sw/uTn7nLP2Put/z9prr7X2SQwh1InASqAJmAyMBcYDbUA7cAR4HngOaAZ2ppTODqUOA4Jar16mTsjnU9THLIYT6j3qPDWNlPI/V29X31T3qV9Ux6tJ/WlBIp14Vl2m1lZb8Tnqwtz+XH54i7olt9eoreqMTOSOComo/kVtrIbyo9Ufqe3qWLVR3azuzg++LR9vzcfvq+/NRO4bAJEz6koLvpWaAgQmAVuAm4DtKaV2YBlwBfBIFuucnOOADmAKsCalJPDriv6xQB3wPeBx9YL+hPskoU4hvEhTvvRCPp7IfccBp4HZ+V4jsBeYASxXa4AVlXN4CwuBreqFfQn1SkJtAL4N7AG2AvuBV/LtscBh4FribSwANgMfBp4G/pRSOgzcCMwdBAmAy4Bt6rRBjtMV6i3qDdl+V+TjLfn4NUtu99QA5kNv2G2sQ/+HHn2zegmwBJgEzAcOAuuB4ymlHVmmFvgK8BFgFvBX4HJgaUrpWfVtwCjgVD5OA94DzMtjTx3A//cosCTPtd6hvl99PbPfpD6S283q17PMSnV2bjeoi8yutwjUWvXThnuuFDcWGXyz4Sr/mzvtVNfl9t1Z7ol8fldRxft43nL13xWQeMOwlF4H/WAWbM9E9ufz/cZCtifL3aduVScPhkTZc6dbWnOK4A99DTY/K38gC/9G/V1uH1NXZLkr1fOGgkDZsyeoT1ZAZF5Pg0xVP5oFHlbvVM+qe9QfG6vovqFUvAcdxqnPFSTxaPfO09WfGK7xP1nouLpK3WG4ytvsb1INDZFLy3ToCx3qzPKOt2alG9Ql6sYspGH7q9TvWu0Is6TPsoJv4wflnf6ZL35LPV+9X12oXmX4+2GFWmOE5v1hb2eHi/KFM+qasoHOM5KV76gb1DnDTGRJwbdxMeoX1O1G6FyrfsaYGzeUCR4wgrnhJJEsufi+cF0N8C8iWhwD3A6sBe4G7gDuyWM+kFLqGE4SObR4qIDoLOCtgK4j/14wXOxydZQReiyuqsa9QP1EgTexKakfB64DJgIX5t+EPM43iaTGlNKJESDxdsJS+sK+pL5KRKsALwOHgKNEmeUUsDqldKhqmvYD9SSRfPWGYxiVip5w1lh0BpOZDRrq4X7M6XQdkSfUAqOJ3HYUUJ+vTQSOjRiDQH8OJdUB19D1db1BVOqOAgeAjVVRrTjO7+f+63XA9UQhYAxB5gKiBNkIfAmYpLallI5XU9OeYKSj/ZFoQ61Tf9bNzl4zQpCp2SavHA6lu0NdUMDFPlkHfBZYRZjNHOBiYDuwDthG5MZNwKYR4FEk5d2LulQ9alQpGtSrjSrf/WVs9zgCBV+LZXvLO3OJThw0MqxLM5GPqavVv6vzh5lAEVNSnVmXUmpVXyJKKE8R5vM34DHgGeBVYCml6t9wEEjA6gKiL6aUnu/stCaz+oD6DXW9USzQiKXWGZHu+6qqfUY26SJYW95pprG/ME09lwVeU39hKRx+ybJ8o4oEphlztAgau3depl6bb/7RrpWHjca+wYtG5je6SgTq83OKoLmnAWoykXvV01mwLZ+fVA+pDxrZ3ga1fogJjFV/X5CA9rZ2GRWPTmyztPfWalT9Dlh6W09YYO+gIIEpRlWlKLbam8tXZxt12HvVI7nDP9SncnujelPZYK+onx8kgWssPgc0agFdHEyXvDlXvK8HvkzET7uIvGIu0EJsoHTmHmeAPwMz1B+qCypQvFb9pLoNeBB4RwW8V6WUWrro3cMDRhHbW4kICmcBuzMZgV8SIfpB4GYikfoUsRFzCbG+PA60EtFwGxHmTyVK+/OBxQystN8MXJFSOtcniUykAfgQEbvUE3sPY4hUcTxwF7EgLiJ2iBYBDwNXD0CxotgPzEkp9ZeulqBOVH9leIynjZJ6u/pVY8+iQ91leLI31WcqsOtK8bI6Y0DUjVrUkW4DXmUpMPttPm6xemhV39WXnn0WxFJKu4md0R1llycD7yZs/fJ8rVop7HZgbkpp76BHMkL0Ow0TWm9EtRvyP1UNUzqnrjWczNDCCM13qjdbCkuah5jALrWpf20GR6RWfadRJdTSvBgsWoywp66qBHogs9j45qNtgIqfMCLlhQ6iYD0kKac6hsjDm4gqyXTgIqCBqKC0AScpfbTVQumjrXM9jVkJ/gfEGHquO3j8DQAAAABJRU5ErkJggg==
[bioconda-link]:http://bioconda.github.io/
[docker-badge]: https://img.shields.io/docker/automated/maxulysse/sarek.svg?logo=docker
[docker-link]: https://hub.docker.com/r/maxulysse/sarek
[gitter-badge]: https://img.shields.io/gitter/room/SciLifeLab/Sarek.svg?logo=gitter&logoColor=white&colorB=4fb99a
[gitter-link]: https://gitter.im/SciLifeLab/Sarek
[license-badge]: https://img.shields.io/github/license/SciLifeLab/Sarek.svg
[license-link]: https://github.com/SciLifeLab/Sarek/blob/master/LICENSE
[nbis-link]: https://www.nbis.se/
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg?logo=data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+PHN2ZyAgIHhtbG5zOmRjPSJodHRwOi8vcHVybC5vcmcvZGMvZWxlbWVudHMvMS4xLyIgICB4bWxuczpjYz0iaHR0cDovL2NyZWF0aXZlY29tbW9ucy5vcmcvbnMjIiAgIHhtbG5zOnJkZj0iaHR0cDovL3d3dy53My5vcmcvMTk5OS8wMi8yMi1yZGYtc3ludGF4LW5zIyIgICB4bWxuczpzdmc9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIiAgIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgICB4bWxuczpzb2RpcG9kaT0iaHR0cDovL3NvZGlwb2RpLnNvdXJjZWZvcmdlLm5ldC9EVEQvc29kaXBvZGktMC5kdGQiICAgeG1sbnM6aW5rc2NhcGU9Imh0dHA6Ly93d3cuaW5rc2NhcGUub3JnL25hbWVzcGFjZXMvaW5rc2NhcGUiICAgd2lkdGg9IjEyLjc5OTIyOG1tIiAgIGhlaWdodD0iMTIuODA0ODA5bW0iICAgdmlld0JveD0iMCAwIDQ1LjM1MTU5NCA0NS4zNzEzNjkiICAgaWQ9InN2Zzc2NTIiICAgdmVyc2lvbj0iMS4xIiAgIGlua3NjYXBlOnZlcnNpb249IjAuOTEgcjEzNzI1IiAgIHNvZGlwb2RpOmRvY25hbWU9Im5leHRmbG93LWZhdmljb24td2hpdGUuc3ZnIj4gIDxkZWZzICAgICBpZD0iZGVmczc2NTQiIC8+ICA8c29kaXBvZGk6bmFtZWR2aWV3ICAgICBpZD0iYmFzZSIgICAgIHBhZ2Vjb2xvcj0iI2ZmZmZmZiIgICAgIGJvcmRlcmNvbG9yPSIjNjY2NjY2IiAgICAgYm9yZGVyb3BhY2l0eT0iMS4wIiAgICAgaW5rc2NhcGU6cGFnZW9wYWNpdHk9IjAuMCIgICAgIGlua3NjYXBlOnBhZ2VzaGFkb3c9IjIiICAgICBpbmtzY2FwZTp6b29tPSI3LjkxOTU5NTkiICAgICBpbmtzY2FwZTpjeD0iMjAuMTEzMjM1IiAgICAgaW5rc2NhcGU6Y3k9IjIzLjE2MzkwOCIgICAgIGlua3NjYXBlOmRvY3VtZW50LXVuaXRzPSJweCIgICAgIGlua3NjYXBlOmN1cnJlbnQtbGF5ZXI9ImxheWVyMSIgICAgIHNob3dncmlkPSJmYWxzZSIgICAgIGZpdC1tYXJnaW4tdG9wPSIwIiAgICAgZml0LW1hcmdpbi1sZWZ0PSIwIiAgICAgZml0LW1hcmdpbi1yaWdodD0iMCIgICAgIGZpdC1tYXJnaW4tYm90dG9tPSIwIiAgICAgaW5rc2NhcGU6d2luZG93LXdpZHRoPSIxOTIwIiAgICAgaW5rc2NhcGU6d2luZG93LWhlaWdodD0iMTAxNSIgICAgIGlua3NjYXBlOndpbmRvdy14PSIwIiAgICAgaW5rc2NhcGU6d2luZG93LXk9IjAiICAgICBpbmtzY2FwZTp3aW5kb3ctbWF4aW1pemVkPSIxIiAvPiAgPG1ldGFkYXRhICAgICBpZD0ibWV0YWRhdGE3NjU3Ij4gICAgPHJkZjpSREY+ICAgICAgPGNjOldvcmsgICAgICAgICByZGY6YWJvdXQ9IiI+ICAgICAgICA8ZGM6Zm9ybWF0PmltYWdlL3N2Zyt4bWw8L2RjOmZvcm1hdD4gICAgICAgIDxkYzp0eXBlICAgICAgICAgICByZGY6cmVzb3VyY2U9Imh0dHA6Ly9wdXJsLm9yZy9kYy9kY21pdHlwZS9TdGlsbEltYWdlIiAvPiAgICAgICAgPGRjOnRpdGxlPjwvZGM6dGl0bGU+ICAgICAgPC9jYzpXb3JrPiAgICA8L3JkZjpSREY+ICA8L21ldGFkYXRhPiAgPGcgICAgIGlua3NjYXBlOmxhYmVsPSJMYXllciAxIiAgICAgaW5rc2NhcGU6Z3JvdXBtb2RlPSJsYXllciIgICAgIGlkPSJsYXllcjEiICAgICB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxMTQuMTA0MzcsLTQ1Mi41MzM2NikiPiAgICA8cGF0aCAgICAgICBzdHlsZT0iZmlsbDojZmZmZmZmIiAgICAgICBkPSJtIC0xMTQuMTA0MzcsNDU1LjQ2NTYyIDAsOC44NjEzMyAwLjIwMzEzLDAuMDYwNSBjIDMuODcyOTMsMS4xMzk0MyA4LjY1MjUxLDQuMzgzMiAxMi44MDA3OCw4LjY4NzUgMC45MTM2MywwLjk0ODAxIDEuOTcyNTY0LDIuMTA2ODQgMi4zNTM1MjQsMi41NzYxOCBsIDAuNjkxNCwwLjg1MzUxIC0wLjg2OTE0LDAuNzc1MzkgYyAtNC4xOTk5MDQsMy43NDE5MyAtOC45NzE5MDQsNi43NjYzNyAtMTQuMTA1NDc0LDguOTQxNDEgLTAuMzA5NzUsMC4xMzEyNCAtMC42OTcyMiwwLjI4MTIzIC0xLjA3NDIyLDAuNDI3NzMgbCAwLDkuMzA0NjkgYyAyLjY1OTkzLC0wLjg3NzkyIDUuMzA2MzksLTEuOTc1IDguMDYwNTUsLTMuMzUxNTYgNC4yNTYyMywtMi4xMjczMiA3LjU0MzI1NCwtNC4yNTc2NCAxMS4wMzcxMTQsLTcuMTU2MjUgMC45MjU4MSwtMC43NjgwOCAxLjgyMTA5LC0xLjUwNzAyIDEuOTkwMjMsLTEuNjQyNTggMC4yNzkzMSwtMC4yMjM4NCAwLjQ5MzMyLC0wLjA1MTQgMi4zMjQyMiwxLjg3ODkxIDYuMjIyNjUsNi41NjA0MSAxMy43ODMzNywxMC43NDQ0MyAyMS45Mzk0NiwxMi4yMjI2NSBsIDAsLTguOTQxNCBjIC00Ljc5NTM3LC0xLjE5NTkgLTkuNDIwMzEsLTMuNjQ1MTEgLTEzLjI1NzgyLC03LjA2NDQ2IC0xLjY4MzUxLC0xLjUwMDA2IC00LjI4NjgxLC00LjM1MDA5IC00LjM5MjU4LC00LjgwODU5IC0wLjA2ODYsLTAuMjk3MyA1LjQ3NDgsLTUuNzA5NzcgNy4yOTQ5MywtNy4xMjMwNSAzLjQ4MjczLC0yLjcwNDI0IDYuNTg4MjUsLTQuMTIwNDIgMTAuMjc1MzksLTQuNjg1NTQgMC4wMjc1LC0wLjAwNCAwLjA1MjcsLTAuMDA4IDAuMDgwMSwtMC4wMTE3IGwgMCwtOC43NSBjIC03LjkzOTI3LDIuMDIxMTQgLTE0Ljg3MDAxLDUuODc3MzggLTIxLjUsMTEuOTQzMzYgbCAtMS42MzA4NiwxLjQ5MjE4IC0yLjk5NjEsLTMuMDA3ODEgYyAtMS42NDc1NiwtMS42NTQ3IC0zLjc0MDI1LC0zLjYwMTU3IC00LjY1MjM0LC00LjMyNjE3IC01LjAwODU1NCwtMy45Nzg5OSAtMTAuMTUyOTU0LC02LjQ5OTIzIC0xNC41NzIyNzQsLTcuMTU2MjUgeiIgICAgICAgaWQ9InBhdGg3NjIwIiAgICAgICBpbmtzY2FwZTpjb25uZWN0b3ItY3VydmF0dXJlPSIwIiAgICAgICBzb2RpcG9kaTpub2RldHlwZXM9ImNjY3NjY2NzY2Nzc2NzY2NzY3NjY2NjY2NzYyIgLz4gIDwvZz48L3N2Zz4=
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
[travis-badge]: https://img.shields.io/travis/SciLifeLab/Sarek.svg?logo=travis
[travis-link]: https://travis-ci.org/SciLifeLab/Sarek
[version-badge]: https://img.shields.io/github/release/SciLifeLab/Sarek.svg?logo=github&logoColor=white
[version-link]: https://github.com/SciLifeLab/Sarek/releases/latest
[zenodo-badge]: https://zenodo.org/badge/54024046.svg
[zenodo-link]: https://zenodo.org/badge/latestdoi/54024046
