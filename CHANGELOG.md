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
