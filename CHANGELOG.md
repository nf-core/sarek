# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Python wrapper script

## [2.0.0] - 2018-03-23
### `Added`
- basic wrapper script
- Abstract, posters and figures
- ROI selector and FreeBayes sanitizer scripts
- New logo and icon for the project
- check for existing tumor/normal channel
- `lib/SarekUtils.groovy` with `checkParams`, `checkParameterList`, `checkParameterExistence` and `isAllowedParams` functions
- some `runOptions` for `docker` (prevent some user right problem)
- This `CHANGELOG`

### `Changed`
- `CAW` is now `Sarek`
- Dissect Workflow in 5 new scripts: `annotate.nf`, `main.nf`, `germlineVC.nf`, `runMultiQC.nf` and `somaticVC.nf`
- `report.html`, `timeline.html` and `trace.html` are generated in `Reports/`
- `--version` is now used to define the workflow version
- most params are now defined in the base.config file instead of in the scripts
- update RELEASE_CHECKLIST.md
- `checkParams`, `checkParameterList`, `checkParameterExistence` and `isAllowedParams` in script functions are now called within `SarekUtils`
- `nf_required_version` is now `params.nfRequiredVersion`
- in `buildReferences.nf` script, channels now begin by `ch_`, and files by `f_`
- use `PublishDir mode: 'link'`` instead of `copy`
- `directoryMap` now contains `params.outDir`
- use Nextflow support of scratch (close #539)
- reordered Travis CI tests
- update documentation
- `MultiQC` version in container from v`1.4` to v`1.5`
- `vepgrch37` container base image from `release_90.6` to `release_92`
- `vepgrch38` container base image from `release_90.6` to `release_92`
- `VEP` version in containers from v`90` to v`91`
- `nucleotidesPerSecond` is now `params.nucleotidesPerSecond`
- default `params.tag` is now `latest` instead of current version, so --tag needs to be specified with the right version to be sure of using the `containers` corresponding

### `Deprecated`
- `standard` profile
- `uppmax-localhost.config` file

### `Removed`
- `scripts/skeleton_batch.sh`
- old data and tsv files
- UPPMAX directories from containers
- `--step` in `annotate.nf`, `germlineVC.nf` and `somatic.nf`
- some `runOptions` for Singularity (binding not needed anymore on UPPMAX)
- `download` profile

### `Fixed`
- Replace `VEP` `--pick` option by `--per_gene` (fix #533)
- use `$PWD` for default `outDir` (fix #530)

## [1.2.5] - 2018-01-18

### `Added`
- Zenodo for DOI
- Delivery README
- Document use of the `--sampleDir` option
- Contributing Guidelines
- Issue Templates
- Release Checklist
- `--outDir`
- `awsbatch` profile
- `aws-batch.config` config file
- `--noBAMQC` params (failing sometimes on Bianca)

### `Changed`
- Update `Nextflow` to `0.26.0` (new fancy report + AWS Batch)
- Extra time on Travis CI testing
- Replace `bundleDir` by `params.genome_base`
- Update `MultiQC` to `1.3` (MEGAQC FTW)
- Move and rename some test files

### `Fixed`
- Version of COSMIC GRCh37 v83
- Write an error message when `--sampleDir` does not find any FASTQ files
- `base.config` for ConcatVCF process
- File specification for recalibrationReport in RecalibrateBam process (got error on AWS Batch)

## [1.2.4] - 2017-10-27

### `Fixed`
- Better CPU requirements for `ConcatVCF` (fix #488)
- Exception handling for `ASCAT` (close #489)
- CPU requirements for `runSingleStrelka` and `runSingleManta` (fix #490)

## [1.2.3] - 2017-10-18

### `Fixed`
- 16 cpus for local executor (fix #475)
- `ASCAT` works for GRCh38 (fix #357)
- Running `Singularity` on /scratch (fix #471)
- No tsv for step `annotate` (fix #480)

## [1.2.2] - 2017-10-06

### `Fixed`
 - Typo in `uppmax-localhost.config` (fix #479)

## [1.2.1] - 2017-10-06

### `Changed`
- `runascat` and `runconvertallelecounts` containers are now replaced by `r-base`
- `willmclaren/ensembl-vep:release_90.5` is now base for `vepgrch37` and `vepgrch38`

### `Removed`
- `vep` container
- `strelka_config.ini` file

### `Fixed`
- Running `Singularity` on /scratch (fix #471)
- Update function to check Nextflow version (fix #472)
- Remove `returnMin()` function (fix #473)

## [1.2.0] - 2017-10-02

### `Changed`
- Fix version for Manuscript

## [1.1] - 2017-09-15

### `Added`
- Singularity possibilities

### `Changed`
- Reports made by default
- Intervals file can be a bed file
- Normal sample preprocessing + HaplotypeCaller is possible
- Better Travis CI tests

### `Fixed`
- Memory requirements

## [1.0] - 2017-02-16

### `Added`
- Docker possibilities

## [0.9] - 2016-11-16

## [0.8] - 2016-11-16

## [0.1] - 2016-04-05
