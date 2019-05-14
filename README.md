# [![Sarek](docs/images/Sarek_logo.png "Sarek")](https://sarek.scilifelab.se/)
[![nf-core](docs/images/nf-core_logo.png "Sarek")](https://nf-co.re/)

**An open-source analysis pipeline to detect germline or somatic variants from whole genome or targeted sequencing**.

> :warning: This pipeline is a work in progress being ported to nf-core from [SciLifeLab/Sarek](https://github/SciLifeLab/Sarek/)

[![Nextflow version][nextflow-badge]](https://www.nextflow.io/)
[![nf-core][nf-core-badge]](https://nf-co.re/)

[![Travis build status][travis-badge]](https://travis-ci.com/nf-core/sarek/)
[![CircleCi build status][circleci-badge]](https://circleci.com/gh/nf-core/sarek/)

[![Install with bioconda][bioconda-badge]](https://bioconda.github.io/)
[![Docker Container available][docker-sarek-badge]](https://hub.docker.com/r/nfcore/sarek/)
[![Install with Singularity][singularity-badge]](https://www.sylabs.io/docs/)

[![Join us on Slack][slack-badge]](https://nfcore.slack.com/messages/CGFUX04HZ/)

## Introduction
<img align="right" title="CAW" src="/docs/images/CAW_logo.png">

Previously known as the Cancer Analysis Workflow (CAW),
Sarek is a workflow designed to run analyses on whole genome or targeted sequencing data from regular samples or tumour / normal pairs and could including additional relapses.

It's built using [Nextflow](https://www.nextflow.io),
a domain specific language for workflow building,
across multiple compute infrastructures in a very portable manner.
Software dependencies are handled using [Conda](https://conda.io/), [Docker](https://www.docker.com) or [Singularity](https://www.sylabs.io/singularity/) - environment/container technologies that provide excellent reproducibility and ease of use.
Thus making installation trivial and results highly reproducible.

It is listed on the [Elixir - Tools and Data Services Registry](https://bio.tools/Sarek), [Dockstore](https://dockstore.org/workflows/github.com/SciLifeLab/Sarek/) and [omicX - Bioinformatics tools](https://omictools.com/sarek-tool).

## Documentation
The nf-core/sarek pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
    * [More information on reference](docs/reference.md)
3. [Running the pipeline](docs/usage.md)
    * [Intervals documentation](docs/INTERVALS.md)
    * [Running the pipeline](docs/USAGE.md)
    * [Command line parameters](docs/PARAMETERS.md)
    * [Examples](docs/USE_CASES.md)
    * [Input files documentation](docs/INPUT.md)
    * [Processes documentation](docs/PROCESS.md)
    * [Documentation about containers](docs/CONTAINERS.md)
4. [Output and how to interpret the results](docs/output.md)
    * [Complementary information about ASCAT](docs/ASCAT.md)
    * [Complementary information about annotations](docs/ANNOTATION.md)
    * [Output documentation structure](docs/OUTPUT.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

## Workflow steps

Sarek is built with several workflow scripts.
A wrapper script contained within the repository makes it easy to run the different workflow scripts as a single job.
To test your installation, follow the [tests documentation.](docs/TESTS.md)

Raw FastQ files or BAM files (unmapped, aligned or recalibrated) can be used as inputs.
You can choose which variant callers to use, plus the pipeline is capable of accommodating additional variant calling software or CNV callers if required.

The worflow steps and tools used are as follows:

1. **Preprocessing** _(based on [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/))_
    * Map reads to Reference
        * [BWA](http://bio-bwa.sourceforge.net/)
    * Mark Duplicates
        * [GATK MarkDuplicates](https://github.com/broadinstitute/gatk)
    * Base (Quality Score) Recalibration
        * [GATK BaseRecalibrator](https://github.com/broadinstitute/gatk)
        * [GATK ApplyBQSR](https://github.com/broadinstitute/gatk)
2. **Germline variant calling**
    * SNVs and small indels
        * [GATK HaplotypeCaller](https://github.com/broadinstitute/gatk)
        * [Strelka2](https://github.com/Illumina/strelka)
    * Structural variants
        * [Manta](https://github.com/Illumina/manta)
3. **Somatic variant calling**
    * SNVs and small indels
        * [MuTect2](https://github.com/broadinstitute/gatk)
        * [Freebayes](https://github.com/ekg/freebayes)
        * [Strelka2](https://github.com/Illumina/strelka)
    * Structural variants
        * [Manta](https://github.com/Illumina/manta)
    * Sample heterogeneity, ploidy and CNVs
        * [ASCAT](https://github.com/Crick-CancerGenomics/ascat)
4. **Annotation**
    * Variant annotation
        * [SnpEff](http://snpeff.sourceforge.net/)
        * [VEP (Variant Effect Predictor)](https://www.ensembl.org/info/docs/tools/vep/index.html)
5. **QC and Reporting**
    * QC
      * [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
      * [Qualimap bamqc](http://qualimap.bioinfo.cipf.es/)
      * [samtools stats](https://www.htslib.org/doc/samtools.html)
      * [GATK MarkDuplicates](https://github.com/broadinstitute/gatk)
      * [bcftools stats](https://samtools.github.io/bcftools/)
      * [VCFtools](https://vcftools.github.io/)
      * [SnpEff](http://snpeff.sourceforge.net/)
    * Reporting
        * [MultiQC](https://multiqc.info/)

## Credits

Sarek was developed at the [National Genomics Infastructure][ngi-link] and [National Bioinformatics Infastructure Sweden][nbis-link] which are both platforms at [SciLifeLab][scilifelab-link], with the support of [The Swedish Childhood Tumor Biobank (Barntumörbanken)][btb-link].

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

For further information or help, don't hesitate to get in touch on [Slack](https://nfcore.slack.com/messages/CGFUX04HZ/) or contact us: maxime.garcia@scilifelab.se, szilveszter.juhos@scilifelab.se

## CHANGELOG

* [CHANGELOG](CHANGELOG.md)

## Aknowledgements
[![Barntumörbanken](docs/images/BTB_logo.png)](https://ki.se/forskning/barntumorbanken-0) | [![SciLifeLab](docs/images/SciLifeLab_logo.png)](https://scilifelab.se)
:-:|:-:
[![National Genomics Infrastructure](docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/) | [![National Bioinformatics Infrastructure Sweden](docs/images/NBIS_logo.png)](https://nbis.se)

[bioconda-badge]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACQAAAAkCAYAAADhAJiYAAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4wUJDSc29Eu32QAACDBJREFUWMO1mGtwVOUZx//P875nz16zSxIgASKxI4URbWurDrYzWlr90GIRQZM14aK0oxVn1A86fvKD44eOnVHH6ThgR4cSLiaUFhzsdFov2NJ4+4C3ip2iYkwg2JIEsmd3z57zvu/TD8o9ASLp82nP7tk9v32e//vcCOdpql3B9tjj1+kV6YY4jn+mlb48tvFUxUqY+FBkon0ABgF8xsz9DemG4NBzh9z5PofOCtGmYLfak6F+oki1OnE8s2Hmtv7D/deIyJMCaT3zh8mC8BaB3tBab4+2RL0AkF2RRbAxmDiQLmqYbgO6hRoKucK8KI6iSlT5KxP/2zk3Xym1W0TW+56/u1qr7hFI01n+W1WxetX22BsBQLdrmB4z5o18+hv1d9QDAEy3gSqqZcR0eGT9SK8VuxBAQESjAEIIktbZZxbMWXBEsbqBiGpnAUpZZxdRG8WqXS06BjPrl7PODtR0ZxOG1w/jigevSFAbrbXWbtNK/4bbOFWLa/clVOKP1tprFavXBdJIIOzZv+fKbCo7BKB2Ln2IiLbOvki30mMAMLBuAE13NZ07ZN5t3pLYxDcqVlMJ1CWQvSLSC+BjYjoCh5lKqfeNNVctnLdw/q6PdvU6uO9hAsbM21yPu/UM3Z5ydTXgX+dnxcmAc+7vWumeS2df+kEYhSUAfcaapcxcIaIREcnnM/niJ//95CoRWUBEB4gpEpEYQHIsOZzmrUv5Mk7Jh/KyLmq4f7oTHvI7fNS21JBflZ9Rqpb+wMTTnTgNoCGXyhWDMPgHga6bVT/rjYHhgZsUq6RA9kfPRztPf1DLXS1TDh45uFRErhRIu4hMGf9EEZRSPzXd5s/ZlVkEXcGpIVPtaqcipSMbLSSiQ0QUisi7GT9zf6VWOZDQiXvCLeG6ZEdSh1tCg0VAa2srPnv6MwBAbkUOpY0lAEDhjoKqRtVZ1tqnjDU3jQtFFMpWSZ2A/AGAXiB/ez5VCSsHjDXkaW9XbOLriajf056xsb0qkUisrmyurMssz6C8qXxOjdStqsPohtEvNVn0no5tvGbcfMdqq+2x7cdDll6RbnXWDRlrPnRwZQIFIpIFkNKsD0Xd0QJcgGV+keFaqbbWOHPneHlKs/6h6TFv87TV0ziKooejONICecU5N8+Jm5nQifcF0uBp78Fj30oUExOG0UWN8rNlZ3rMXZ7y/jJunhK7pOXuFuKRykjGOnsZGM9m/Mwa3/M3Q9AV2ej1DfdvKJT3lncfS2BRdzRhINN9IiNbsWuIaLy68f2Dwwez7JybTURzxUm/0grN+ea7AQwyuHn1U6vvLVxW0APrBnCh5hU9pBKpzwFsHTsP4GoRmcHJRNICWDanec5zQTXo7h/unw/Cj40zDwEIj3QdMZgEi7tjlDeWDYHeAeDO5JEUM7dwzdY+Fcht+wb3vWusmZfQiTtvvvzmJUy8J5/K78QkW1Oh6U9ENDLWZ05cM3Ebd4CwAcCQiBgmTilWr13SfEnn3if2hvg/GLVRn4hcNMbxf4i10jc467Szbro4mWmdnSKQpQcOH1gOAINHBycN5HDl8Jd1jDgYR0fNuiHXsHuoNJTxtV/1tIdStcQEuiiIghEAaM43TxpQY7rxGJCzsGPWEj14ZNAj0K2xjQOEcCAoxWpL2kvvKaF0vM5NhqWWp1DdVEVs49w4te0L1qQ3kNAmAFkQEkz8QcpP/ToIg93ZVdnGyYIBgOqmKgqrCzMApMeOmBxkrfS3mPnlXDL3o4yf+dWU9JTOoBLsEMjMMAyXTLagR4PRxQDqxym0gzo2cdXB/a5cK79GoLtFpJWYniOhFhCm1a2s06NHRw1euDCQdGcaSikVhMF3zujDvkpV1tp+ZsV9AN62ztZnk9mGmqk9y+AOEBpSXuq90WCU66fXX7BnKpsrCKNwFoD2cbzzJjMPcNJLlgH8K+kl7wvC4F5r7cVO3DcYPL8Ull70k/7Dw78dxrFGbqKW7EieSHxwj4tIYZxb38wms2UudZWsVnp77OKPBNLGzAdBiK3YFsXqUwBrACC9PJ39OgIPt4R4ZO0jpIv6SWvtsrOMSS8c3XDUMQDEz8c7nHPfdOJgrZ1BoC8AeE5cRiBlbufrIhM9drwNue3cbUimM3P89aN/e/QJY83947YoSr9quk3vKWNQPpN/H0CY8BLvOHFzScgSEZx17zHxXGPNGl3U6/wOf1pduu64KOtW1o0ZHmZWiWJitiqqbWeDISJrus2Np4xByc4kws0h/A7/+lpce0krPWCd9UVkaiFTeLpULYXE1OSsu4aI+jTrt5y4fTOmzNjet7bvlEKZX5WfVg7LS6zY7wJoO2eTz2qx6TE7uZ3hetyJJj/VmQKBWCt9bVANHnJw0xIq8ZKnvVdiG/vW2icdXIti9R8RqXPifAiOAqhopYcBkLGmXiBpAA3jHO3TvfO4bJUHTh6txxwUcytyC6tx9QZn3cWseFkqkWoLwmC7Zj1gnZ0KQuCca7iQNMDg7e73buk5Z3sAaMw2vu3ENTm4orHGIyLvq37FBxARSF0IDBE9fgzmZA2OvWy4vR771+4vr1+1frXH3i1EhHJY/jkzf2ydrQehZp3Nf00Qq5RaLFvlAQDIr8xjtGv0/NYxiWICUXeE6fdM94eHhp8B0BjbeJFiVbLO5iZaVxWrXbbHLjp51TPhhZXf6aO2+ctkqNv1t524TqXU3NjEi88DwhLR6wDe0krviJ+PewFAFRVst/16G7QzxL48p0ITZqyzs4moBYImZm7VWjdBILGNB0WkD8CgiHzOxAN5Px8MdQ2d90rvf5LI8eoQ7hrVAAAAAElFTkSuQmCC
[btb-link]: https://ki.se/forskning/barntumorbanken-0
[circleci-badge]: https://img.shields.io/circleci/project/github/nf-core/sarek.svg?logo=circleci
[docker-sarek-badge]: https://img.shields.io/docker/automated/nfcore/sarek.svg?logo=docker
[docker-snpeff-badge]: https://img.shields.io/docker/automated/nfcore/sareksnpeff.svg?logo=docker
[docker-vep-badge]: https://img.shields.io/docker/automated/nfcore/sarekvep.svg?logo=docker
[nbis-link]: https://nbis.se
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A519.04.0-brightgreen.svg?logo=data:image/png;base64,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
[nf-core-badge]: https://img.shields.io/badge/nf--core-pipeline-brightgreen.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEAAAABACAYAAACqaXHeAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAAOxAAADsQBlSsOGwAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAAg4SURBVHiczVt7UFTXGf99Zx88XRFFKkiAXUEEjLW2hCE1SKZptI/RgYmmbTJFoJi2No59mnYybtI6aaaZaSZJU19IahunxQBJnJppawdtdExMOzoKiBFWfBaKPNzw3N17vv6B6IrC7t57d29+/yx7zvd953d+3HPOPd85S/gUYcWKFWbZ3f1lElgH8BKAsgBEAegC+ByBmsYgXj/e0tKnV5tmvQJpRUlu7gPyevduEsgfLyH/6lSAUhl42MLcA2CHXu1+KgQozs1dL4l3gSEZVC/ApxiwAtgCwHLbkg7EeTxv6Nm24QKU5OU8KsG7AHQKptVNra0tE3Ur8nNjmflHALUB8vkjLWf/DID1bJ8Cm4QPeXl51jngswAnksJLD7e1dU62eXhJdurs7I+79u+HYgDF8KI4P/fp4rxFvCI/98dGcRBGNQwAYN4IYHiMaY9RFAwToGRxzjIAWQAf0nNZCxWGCSAligEAREeN4gAYOgTE4pt/tBvHwUABGJwGAMSixygOgIECEHDf+KdP13U9VBi5CiQBAEgkGMghci9Cy3ZUW3pSTfHmYWXGnNNtppjGnhYAMb7M6Neul2WPeG3xORCYR4wUBrwMKAT4AFwF6BwI/7Ky5+9tpXt79eSlvwBOp7Dff6mIiIpYogCEHADJAOZMmAivMhD3YVdz3MkeyGyyIoruH05L+mjQMW95gOg+AjWQUJ5rX1Pbqgdd3QTIaCxPMMG0kSUqQMgM1k/4lP7YjmunY7r689y56ZfHkmYuDcLNC8JW15qaX4O07Q10EcDRWLGKGW8ANFdtDOFV+uJdXWfcC1OLg/ei7a7S3d9V2yaggwCOhsrHGXgTekyoDAaFxolBGy+U7v6d2iY1CWBvqPgiQP/E+N7dKLgVNuVcLNv5XzXOqv9ry3ZUWwDaDmM7DwA2AZ9TrbNqAQaSlG8CyFPrryeIaF1GbXm0Gl/VAjBQodY3DJhpsolVahxVCZDRWJ4BINCaHVkIPKjOTQVImpbD4HTaZDDTZ9X4qUqKCqAglLeP/NjojyuSE64lWyzmfwwM0r7egWUehVWN2SlByFLjpkoAJiwKxs4qaGy7PfXD/Nio5QCyASDrM4momJtw7lvnr8RdHvPODxRjUXTU+SJbTHdKlFXp8/rEmeFR2zH3cJ4ymTvfftUOBaoeY3tDZTMCrwDdhbGxT2zPTD/lsYzaIcX3AHx7otIHdH61tdPW71MSJzuaibxPz5t9vCzRlmYR5P9a3UXAnzq8nuPrWq9kEvEW+O0xHoybFf/HR18aCqUvas8FkgEaAXgvMTWxkD4CrWbGExgX9byZeeW+la+69o3bXwdwordj20lmevlmwxk1C1I/KG27VOgfuCA+tvm3GfNiLQIP+RVLEL/4ifA+n5npHJ0onF9XWWs10zaANwAg6wwlCkBIAoT8BOTVPWYdMdn2QTH90LV25yX/OnvDd74GKJuFjx9vX1t7z0xPb8e2I8x0q3PbrvSceKfPXQAAP5ufdKQs0VaEO06DwEz0ZJLjmTen4uRorFzDjNqFKbHp7xW+6g6lPxGfyXs7Xihj5rcmvnuYO4rPuNJfcaQc+0JczF0bIQJ2zF7w86cCxbW/VZ41q9fa+Z8NO72h8Im4ANzitPZGWbsB3MoEjUpuixaUc5cxYUgK84K5mT/tChefiKfEKM/pAeGQf9k9Ow8AEnvD2XnAqJwg05FgzCTzznBTMUQAAT4ZyIYJF+Zm/+JU+LkYAJ/J3BHIhiQOR4CKMQIkZQz/D5j+uJsJuiQ9A8EQAYicEsC0y5UANUeCi5EHI3K6ShLyWiRIGCIAc50J47e/poQk8kSCiyEC9HR2JgEwGdH2ZBgzB/i8KYFshBSqtrehwqBJcOIu4NSQkKoPWUKBMXMA0ZLARpwaASoGPQGMIM7/KDv8TAwQgJucZgCfC8I0qLSbVkRcgL40cwGAmUGYLgw3F8CIIUCmR4K0TO0977SFlQsMEECy/FKQpiSFOezzQEQF6Gl7cQaBHgjWnuQUiZKbcLy9Pi2zofI9LZw0CZD8tyfjQrEns68MdyY8AzhgSgEy6ysfYik+ImClva76vlB4+EOTAHFDlnfHM8GBwewUAH1/UvHoPY1vQZTcVVT3mMneULWFxtNqyQAgBateMTQOAcoF5LuZjVWvZB+onvbVta/dsgnA5xnoPzs8+n616+rZkuYL3raRsWmuynJRz/lffX3i24K31+c6zLajAL8AvyfJRDJDdQ/UOgKAvaHSi9uHK5+AUc+gvzKLZhlt6rYKtyJGouYsjotdum7uzKcODQwmHr4xtNjHfMcwKJwRc+bZ+cljSRbTUty9SRo80Htjwy+vXl8F4Bv3qAeYN7vK9ryspg8aBagaBjhGSwx/JFhM/Y/MiDvniIkaswjC5VGv+ejgUEr7iGfaW2cM/ORCac1LatrU+JMZHgSgmwADXmXW/j53YWDLO0EM1ckTbXMA4aImf51AZpxQ66tNAKZ/a/LXB6c6VteovnKvSQAWfFCLvx4QIKc2fw244HEfBBAwxx8uEOO19tLd72iJoW0IrN2vALxVUwyVYMJf0hLTNmuNo8vpsL2+shaEcj1iBQE3GM+4Smt+r/WiNKDTZsiluKsI2BWkuSTCs2aLKUlhUwqBNgG4EtiNrgLYaoXX7iqreV2PzgM63w+4eUPkN8CUuzgfAVUdpTV/8C9ccPAHUTwyVM1EmwA4bhYPgugDSByTpDR1nk5/H07ntIcpaqD/BQkGORqrSkDyK5KpiDC+oyPiE8R4rr10z/Hp3NPrq+fBIqMvjt64ND7HhBf/B8qHpAjy2+0qAAAAAElFTkSuQmCC
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://scilifelab.se
[singularity-badge]: https://img.shields.io/badge/use%20with-singularity-purple.svg?logo=data:image/png;base64,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
[slack-badge]: https://img.shields.io/badge/slack-nfcore/sarek-blue.svg?logo=slack
[travis-badge]: https://img.shields.io/travis/nf-core/sarek.svg?logo=travis
