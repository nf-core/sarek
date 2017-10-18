[![](https://raw.githubusercontent.com/SciLifeLab/CAW/master/doc/images/CAW_logo.png "CAW")][caw-site-link]

# Cancer Analysis Workflow

[![caw version][version-badge]][version-link] [![License][license-badge]][license-link] [![nextflow version][nextflow-badge]][nextflow-link] [![Join the chat at https://gitter.im/SciLifeLab/CAW][gitter-badge]][gitter-link] [![Travis status][travis-badge]][travis-link]

CAW is a complete open source pipeline to detect somatic variants from WGS data developed at the [National Genomics Infastructure][ngi-link] at [SciLifeLab Stockholm][scilifelab-stockholm-link], Sweden and [National Bioinformatics Infastructure Sweden][nbis-link] at [SciLifeLab][scilifelab-link].

The pipeline uses [Nextflow][nextflow-link], a bioinformatics domain specific language for workflow building and [Singularity](http://singularity.lbl.gov/), a container technology specific for high-performance computing.

This pipeline is primarily used with cluster on the Swedish [UPPMAX systems](https://www.uppmax.uu.se/). However, the pipeline should be able to run on any system that supports Nextflow. The pipeline comes with some configuration for different systems. See the [documentation](#documentation) for more information.

We utilize [GATK best practices](https://software.broadinstitute.org/gatk/best-practices/) to align, realign and recalibrate short-read data in parallel for both normal and tumor sample. After these preprocessing steps, several somatic variant callers scan the resulting BAM files: [MuTect1][mutect1-link], [MuTect2][gatk-link] and [Strelka][strelka-link] are used to find somatic SNVs and small indels, also [GATK HaplotyeCaller][gatk-link] for both the normal and the tumor sample. For structural variants we use [Manta][manta-link]. Furthermore, we are applying [ASCAT][ascat-link] to estimate sample heterogeneity, ploidy and CNVs.

The pipeline can begin the analysis either from raw FASTQ files, only from the realignment step, or directly with any subset of variant callers using recalibrated BAM files. At the end of the analysis the resulting VCF files are merged to facilitate further downstream processing, though results from each caller are also retained. The flow is capable of accommodating additional variant calling software or CNV callers. It is also prepared to process normal, tumor and several relapse samples.

Besides variant calls, the workflow provides quality controls presented by [MultiQC][multiqc-link].

The [containers](containers) directory contains building rules for containers for all CAW processes.

This pipeline is listed on [Elixir - Tools and Data Services Registry](https://bio.tools/CAW).

## Documentation

The CAW pipeline comes with documentation about the pipeline, found in the `doc/` directory:

01. [Installation documentation](doc/INSTALL.md)
02. [Installation documentation specific for `milou`](doc/INSTALL_MILOU.md)
03. [Installation documentation specific for `bianca`](doc/INSTALL_BIANCA.md)
04. [Tests documentation](doc/TESTS.md)
05. [Reference files documentation](doc/REFERENCES.md)
06. [Configuration and profiles documentation](doc/CONFIG.md)
07. [Intervals documentation](doc/INTERVALS.md)
08. [Running the pipeline](doc/USAGE.md)
09. [Examples](doc/USE_CASES.md)
10. [TSV file documentation](doc/TSV.md)
11. [Processes documentation](doc/PROCESS.md)
12. [Documentation about containers](doc/CONTAINERS.md)
13. [Documentation about building](doc/BUILD.md)
14. [More information about ASCAT](doc/ASCAT.md)
15. [Folder structure](doc/FOLDER.md)

For further information/help contact: maxime.garcia@scilifelab.se, szilveszter.juhos@scilifelab.se or join the gitter chat: [gitter.im/SciLifeLab/CAW][gitter-link].

## Authors

- [Sebastian DiLorenzo](https://github.com/Sebastian-D)
- [Jesper Eisfeldt](https://github.com/J35P312)
- [Maxime Garcia](https://github.com/MaxUlysse)
- [Szilveszter Juhos](https://github.com/szilvajuhos)
- [Max Käller](https://github.com/gulfshores)
- [Malin Larsson](https://github.com/malinlarsson)
- [Marcel Martin](https://github.com/marcelm)
- [Björn Nystedt](https://github.com/bjornnystedt)
- [Pall Olason](https://github.com/pallolason)
- [Pelin Sahlén](https://github.com/pelinakan)

--------------------------------------------------------------------------------

[![](https://raw.githubusercontent.com/SciLifeLab/CAW/master/doc/images/SciLifeLab_logo.png "SciLifeLab")][scilifelab-link]
[![](https://raw.githubusercontent.com/SciLifeLab/CAW/master/doc/images/NGI_logo.png "NGI")][ngi-link]
[![](https://raw.githubusercontent.com/SciLifeLab/CAW/master/doc/images/NBIS_logo.png "NBIS")][nbis-link]

[ascat-link]: https://github.com/Crick-CancerGenomics/ascat
[caw-site-link]: http://opensource.scilifelab.se/projects/caw/
[gatk-link]: https://github.com/broadgsa/gatk-protected
[gitter-badge]: https://badges.gitter.im/SciLifeLab/CAW.svg
[gitter-link]: https://gitter.im/SciLifeLab/CAW
[license-badge]: https://img.shields.io/github/license/SciLifeLab/CAW.svg
[license-link]: https://github.com/SciLifeLab/CAW/blob/master/LICENSE
[manta-link]: https://github.com/Illumina/manta
[multiqc-link]: https://github.com/ewels/MultiQC/
[mutect1-link]: https://github.com/broadinstitute/mutect
[nbis-link]: https://www.nbis.se/
[nextflow-badge]: https://img.shields.io/badge/nextflow-%E2%89%A50.25.0-brightgreen.svg
[nextflow-link]: https://www.nextflow.io/
[ngi-link]: https://ngisweden.scilifelab.se/
[scilifelab-link]: https://www.scilifelab.se/
[scilifelab-stockholm-link]: https://www.scilifelab.se/platforms/ngi/
[strelka-link]: https://github.com/Illumina/strelka
[travis-badge]: https://api.travis-ci.org/SciLifeLab/CAW.svg
[travis-link]: https://travis-ci.org/SciLifeLab/CAW
[version-badge]: https://img.shields.io/github/release/SciLifeLab/CAW.svg
[version-link]: https://github.com/SciLifeLab/CAW/releases/latest
