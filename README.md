[![](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/images/Sarek_logo.png "Sarek")][sarek-site-link]

# Sarek

[![sarek version][version-badge]][version-link]
[![License][license-badge]][license-link]
[![nextflow version][nextflow-badge]][nextflow-link]
[![Join the chat at https://gitter.im/SciLifeLab/Sarek][gitter-badge]][gitter-link]
[![Travis status][travis-badge]][travis-link]
[![DOI][zenodo-badge]][zenodo-link]

[![](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/images/CAW_logo.png "CAW")][sarek-site-link]

Previously known as the Cancer Analysis Workflow (CAW), Sarek is a complete open source pipeline to detect germline, or somatic variants from WGS data developed at the [National Genomics Infastructure][ngi-link] at [SciLifeLab Stockholm][scilifelab-stockholm-link] and [National Bioinformatics Infastructure Sweden][nbis-link] at [SciLifeLab][scilifelab-link].

The pipeline uses [Nextflow][nextflow-link], a bioinformatics domain specific language for workflow building and [Singularity](http://singularity.lbl.gov/), a container technology specific for high-performance computing.

This pipeline is primarily used with cluster on the Swedish [UPPMAX systems](https://www.uppmax.uu.se/).
However, the pipeline should be able to run on any system that supports Nextflow.
The pipeline comes with some configuration for different systems.
See the [documentation](#documentation) for more information.

Sarek is based on [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/) to align, realign and recalibrate short-read data (which is done in parallel for tumor/normal pair samples).
After these preprocessing steps, several variant callers scan the resulting BAM files:
[GATK HaplotyeCaller][gatk-link] and [Strelka][strelka-link] are used to find germline SNVs and small indels (also used on tumor samples).
[MuTect1][mutect1-link], [MuTect2][gatk-link], [Freebayes][freebayes-link] and [Strelka][strelka-link] are used to find somatic SNVs and small indels.
For structural variants (germline and somatic) we use [Manta][manta-link].
Furthermore, we are applying [ASCAT][ascat-link] to estimate sample heterogeneity, ploidy and CNVs.

The pipeline is prepared to process normal or tumor/normal pairs (and several relapse samples).
It can begin the analysis either from raw FASTQ files, only from the realignment step, or directly with any subset of variant callers using recalibrated BAM files.
At the end of the analysis the resulting VCF files and results from each caller are also retained.
And snpEff and/or VEP can be used to annotate them.

The flow is capable of accommodating additional variant calling software or CNV callers.

Besides variant calls, the workflow provides quality controls presented by [MultiQC][multiqc-link].

The [containers](containers) directory contains building rules for containers for all Sarek processes.

This pipeline is listed on [Elixir - Tools and Data Services Registry](https://bio.tools/Sarek).

## Documentation

The Sarek pipeline comes with documentation about the pipeline, found in the `doc/` directory:

01. [Installation documentation](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/INSTALL.md)
02. [Installation documentation specific for `rackham`](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/INSTALL_RACKHAM.md)
03. [Installation documentation specific for `bianca`](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/INSTALL_BIANCA.md)
04. [Tests documentation](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/TESTS.md)
05. [Reference files documentation](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/REFERENCES.md)
06. [Configuration and profiles documentation](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/CONFIG.md)
07. [Intervals documentation](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/INTERVALS.md)
08. [Running the pipeline](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/USAGE.md)
09. [Examples](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/USE_CASES.md)
10. [TSV file documentation](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/TSV.md)
11. [Processes documentation](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/PROCESS.md)
12. [Documentation about containers](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/CONTAINERS.md)
13. [Documentation about building](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/BUILD.md)
14. [More information about ASCAT](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/ASCAT.md)
15. [Folder structure](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/FOLDER.md)

## Contributions & Support

- [Contributions guidelines](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/.github/CONTRIBUTING.md)
For further information/help, don't hesitate to get in touch on [Gitter][gitter-link] or contact us: maxime.garcia@scilifelab.se, szilveszter.juhos@scilifelab.se

## Authors

- [Sebastian DiLorenzo](https://github.com/Sebastian-D)
- [Jesper Eisfeldt](https://github.com/J35P312)
- [Phil Ewels](https://github.com/ewels)
- [Maxime Garcia](https://github.com/MaxUlysse)
- [Szilveszter Juhos](https://github.com/szilvajuhos)
- [Max Käller](https://github.com/gulfshores)
- [Malin Larsson](https://github.com/malinlarsson)
- [Marcel Martin](https://github.com/marcelm)
- [Björn Nystedt](https://github.com/bjornnystedt)
- [Pall Olason](https://github.com/pallolason)
- [Pelin Sahlén](https://github.com/pelinakan)

--------------------------------------------------------------------------------

[![](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/images/NGI_logo.png "NGI")][ngi-link]
[![](https://raw.githubusercontent.com/SciLifeLab/Sarek/master/doc/images/NBIS_logo.png "NBIS")][nbis-link]

[ascat-link]: https://github.com/Crick-CancerGenomics/ascat
[freebayes-link]: https://github.com/ekg/freebayes
[gatk-link]: https://github.com/broadgsa/gatk-protected
[gitter-badge]: https://badges.gitter.im/SciLifeLab/Sarek.svg
[gitter-link]: https://gitter.im/SciLifeLab/Sarek
[license-badge]: https://img.shields.io/github/license/SciLifeLab/Sarek.svg
[license-link]: https://github.com/SciLifeLab/Sarek/blob/master/LICENSE
[manta-link]: https://github.com/Illumina/manta
[multiqc-link]: https://github.com/ewels/MultiQC/
[mutect1-link]: https://github.com/broadinstitute/mutect
[nbis-link]: https://www.nbis.se/
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A50.25.0-brightgreen.svg
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[sarek-site-link]: http://opensource.scilifelab.se/projects/sarek/
[scilifelab-link]: https://www.scilifelab.se/
[scilifelab-stockholm-link]: https://www.scilifelab.se/facilities/ngi-stockholm/
[strelka-link]: https://github.com/Illumina/strelka
[travis-badge]: https://api.travis-ci.org/SciLifeLab/Sarek.svg
[travis-link]: https://travis-ci.org/SciLifeLab/Sarek
[version-badge]: https://img.shields.io/github/release/SciLifeLab/Sarek.svg
[version-link]: https://github.com/SciLifeLab/Sarek/releases/latest
[zenodo-badge]: https://zenodo.org/badge/54024046.svg
[zenodo-link]: https://zenodo.org/badge/latestdoi/54024046
