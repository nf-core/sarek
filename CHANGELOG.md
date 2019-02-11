# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### `Added`
-   [#712](https://github.com/SciLifeLab/Sarek/pull/712), [#718](https://github.com/SciLifeLab/Sarek/pull/718) - Added possibilities to run Sarek with `conda`

### `Changed`

-   [#710](https://github.com/SciLifeLab/Sarek/pull/710) - Improve release checklist and script
-   [#711](https://github.com/SciLifeLab/Sarek/pull/711) - Improve configuration priorities
-   [#719](https://github.com/SciLifeLab/Sarek/pull/719) - `vepCacheVersion` is now defined in `conf/genomes.config` or `conf/igenomes.config`
-   [#719](https://github.com/SciLifeLab/Sarek/pull/719) - `snpeff` and `vep` containers are now built with conda
-   [#716](https://github.com/SciLifeLab/Sarek/pull/716) - Update paths to containers and iGenomes

### `Added`
-   [#719](https://github.com/SciLifeLab/Sarek/pull/719) - Possibility to use cache wen annotating with `snpEff` and `VEP`
-   [#719](https://github.com/SciLifeLab/Sarek/pull/719) - New `--annotation_cache`, `--snpEff_cache`, `--vep_cache` parameters
-   [#719](https://github.com/SciLifeLab/Sarek/pull/719) - Helper script to download `snpeff` and `VEP` cache files
-   [#719](https://github.com/SciLifeLab/Sarek/pull/719) - Annotation documentation

### `Removed`
-   [#715](https://github.com/SciLifeLab/Sarek/pull/715) - Remove `defReferencesFiles` function from `buildReferences.nf`
-   [#719](https://github.com/SciLifeLab/Sarek/pull/719) - `snpEff` base container is no longer used
-   [#721](https://github.com/SciLifeLab/Sarek/pull/721) - Remove COSMIC docs

### `Fixed`
-   [#720](https://github.com/SciLifeLab/Sarek/pull/720) - bamQC is now run on the recalibrated bams, and not after MarkDuplicates

## [2.2.2] - 2018-12-19

### `Added`

-   [#671](https://github.com/SciLifeLab/Sarek/pull/671) - New `publishDirMode` param and docs
-   [#673](https://github.com/SciLifeLab/Sarek/pull/673), [#675](https://github.com/SciLifeLab/Sarek/pull/675),  [#676](https://github.com/SciLifeLab/Sarek/pull/676) - Profiles for BinAC and CFC clusters in Tübingen
-   [#679](https://github.com/SciLifeLab/Sarek/pull/679) - Add container for `CreateIntervalBeds`
-   [#692](https://github.com/SciLifeLab/Sarek/pull/692), [#697](https://github.com/SciLifeLab/Sarek/pull/697) - Add AWS iGenomes possibilities (within `conf/igenomes.conf`)
-   [#694](https://github.com/SciLifeLab/Sarek/pull/694) - Add monochrome and grey logos for light or dark background
-   [#698](https://github.com/SciLifeLab/Sarek/pull/698) - Add btb profile for munin server
-   [#702](https://github.com/SciLifeLab/Sarek/pull/702) - Add font-ttf-dejavu-sans-mono `2.37` and fontconfig `2.12.6` to container
-   [#722](https://github.com/SciLifeLab/Sarek/pull/722) - Update `Sarek-data` submodule
-   [#722](https://github.com/SciLifeLab/Sarek/pull/722) - Add path to ASCAT `.gc` file in `igenomes.config`

### `Changed`

-   [#678](https://github.com/SciLifeLab/Sarek/pull/678) - Changing VEP to v92 and adjusting CPUs for VEP
-   [#663](https://github.com/SciLifeLab/Sarek/pull/663) - Update `do_release.sh` script
-   [#671](https://github.com/SciLifeLab/Sarek/pull/671) - publishDir modes are now params
-   [#677](https://github.com/SciLifeLab/Sarek/pull/677), [#698](https://github.com/SciLifeLab/Sarek/pull/698), [#703](https://github.com/SciLifeLab/Sarek/pull/703) - Update docs
-   [#679](https://github.com/SciLifeLab/Sarek/pull/679) - Update old awsbatch configuration
-   [#682](https://github.com/SciLifeLab/Sarek/pull/682) - Specifications for memory and cpus for awsbatch
-   [#693](https://github.com/SciLifeLab/Sarek/pull/693) - Qualimap bamQC is now ran after mapping and after recalibration for better QC
-   [#700](https://github.com/SciLifeLab/Sarek/pull/700) - Update GATK to `4.0.9.0`
-   [#702](https://github.com/SciLifeLab/Sarek/pull/702) - Update FastQC to `0.11.8`
-   [#705](https://github.com/SciLifeLab/Sarek/pull/705) - Change `--TMP_DIR` by `--tmp-dir` for GATK `4.0.9.0` BaseRecalibrator
-   [#706](https://github.com/SciLifeLab/Sarek/pull/706) - Update TravisCI testing

### `Fixed`

-   [#665](https://github.com/SciLifeLab/Sarek/pull/665) - Input bam file now has always the same name (whether it is from a single fastq pair or multiple) in the MarkDuplicates process, so metrics too
-   [#672](https://github.com/SciLifeLab/Sarek/pull/672) - process `PullSingularityContainers` from `buildContainers.nf` now expect a file with the correct `.simg` extension for singularity images, and no longer the `.img` one.
-   [#679](https://github.com/SciLifeLab/Sarek/pull/679) - Add publishDirMode for `germlineVC.nf`
-   [#700](https://github.com/SciLifeLab/Sarek/pull/700) - Fix [#699](https://github.com/SciLifeLab/Sarek/issues/699) missing DP in the FORMAT column VCFs for MuTect2
-   [#702](https://github.com/SciLifeLab/Sarek/pull/702) - Fix [#701](https://github.com/SciLifeLab/Sarek/issues/701)
-   [#705](https://github.com/SciLifeLab/Sarek/pull/705) - Fix [#704](https://github.com/SciLifeLab/Sarek/issues/704)

## [2.2.1] - 2018-10-04

### `Changed`

-   [#646](https://github.com/SciLifeLab/Sarek/pull/646) - Update [`pathfindr`](https://github.com/NBISweden/pathfindr) submodule
-   [#659](https://github.com/SciLifeLab/Sarek/pull/659) - Update Nextflow to `0.32.0`
-   [#660](https://github.com/SciLifeLab/Sarek/pull/660) - Update docs

### `Fixed`

-   [#657](https://github.com/SciLifeLab/Sarek/pull/657) - Fix `RunMultiQC.nf` bug
-   [#659](https://github.com/SciLifeLab/Sarek/pull/659) - Fix bugs due to updating Nextflow

## [2.2.0] - Skårki - 2018-09-21

### `Added`

-   [#613](https://github.com/SciLifeLab/Sarek/pull/613) - Add Issue Templates (bug report and feature request)
-   [#614](https://github.com/SciLifeLab/Sarek/pull/614) - Add PR Template
-   [#615](https://github.com/SciLifeLab/Sarek/pull/615) - Add presentation
-   [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Update documentation
-   [#620](https://github.com/SciLifeLab/Sarek/pull/620) - Add `tmp/` to `.gitignore`
-   [#625](https://github.com/SciLifeLab/Sarek/pull/625) - Add [`pathfindr`](https://github.com/NBISweden/pathfindr) as a submodule
-   [#639](https://github.com/SciLifeLab/Sarek/pull/639) - Add a complete example analysis to docs
-   [#635](https://github.com/SciLifeLab/Sarek/pull/635) - To process targeted sequencing with a target BED
-   [#640](https://github.com/SciLifeLab/Sarek/pull/640), [#642](https://github.com/SciLifeLab/Sarek/pull/642) - Add helper script for changing version number

### `Changed`

-   [#608](https://github.com/SciLifeLab/Sarek/pull/608) - Update Nextflow required version
-   [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Update CHANGELOG
-   [#615](https://github.com/SciLifeLab/Sarek/pull/615) - Use `splitCsv` instead of `readlines`
-   [#621](https://github.com/SciLifeLab/Sarek/pull/621), [#638](https://github.com/SciLifeLab/Sarek/pull/638) - Improve install script
-   [#621](https://github.com/SciLifeLab/Sarek/pull/621), [#638](https://github.com/SciLifeLab/Sarek/pull/638) - Simplify tests
-   [#627](https://github.com/SciLifeLab/Sarek/pull/627), [#629](https://github.com/SciLifeLab/Sarek/pull/629), [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Refactor docs
-   [#629](https://github.com/SciLifeLab/Sarek/pull/629) - Refactor config
-   [#632](https://github.com/SciLifeLab/Sarek/pull/632) - Use 2 threads and 2 cpus FastQC processes
-   [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Update tool version gathering
-   [#638](https://github.com/SciLifeLab/Sarek/pull/638) - Use correct `.simg` extension for Singularity images
-   [#639](https://github.com/SciLifeLab/Sarek/pull/639) - Smaller refactoring of the docs
-   [#640](https://github.com/SciLifeLab/Sarek/pull/640) - Update RELEASE_CHECKLIST
-   [#642](https://github.com/SciLifeLab/Sarek/pull/642) - Update conda channel order priorities
-   [#642](https://github.com/SciLifeLab/Sarek/pull/642) - MultiQC 1.5 -> 1.6
-   [#642](https://github.com/SciLifeLab/Sarek/pull/642) - Qualimap 2.2.2a -> 2.2.2b
-   [#642](https://github.com/SciLifeLab/Sarek/pull/642) - VCFanno 0.2.8 -> 0.3.0
-   [#642](https://github.com/SciLifeLab/Sarek/pull/642) - VCFtools 0.1.15 -> 0.1.16

### `Removed`

-   [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Remove old Issue Template
-   [#629](https://github.com/SciLifeLab/Sarek/pull/629) - Remove old Dockerfiles
-   [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Remove old comments

### `Fixed`

-   [#621](https://github.com/SciLifeLab/Sarek/pull/621) - Fix VEP tests
-   [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Fix links in MD files

## [2.1.0] - Ruotes - 2018-08-14

### `Added`

-   [#555](https://github.com/SciLifeLab/Sarek/pull/555) - `snpEff` output into `VEP`
-   [#556](https://github.com/SciLifeLab/Sarek/pull/556) - `Strelka` Best Practices
-   [#563](https://github.com/SciLifeLab/Sarek/pull/563) - Use `SnpEFF` reports in `MultiQC`
-   [#568](https://github.com/SciLifeLab/Sarek/pull/568) - `VCFTools` process `RunVcftools` for QC
-   [#574](https://github.com/SciLifeLab/Sarek/pull/574), [#580](https://github.com/SciLifeLab/Sarek/pull/580) - Abstracts for NPMI, JOBIM and  EACR25
-   [#577](https://github.com/SciLifeLab/Sarek/pull/577) - New repository for testing: [Sarek-data](https://github.com/SciLifeLab/Sarek-data)
-   [#595](https://github.com/SciLifeLab/Sarek/pull/595) - New library `QC` for functions `bamQC`, `bcftools`, `samtoolsStats`, `vcftools`, `getVersionBCFtools`, `getVersionGATK`, `getVersionManta`, `getVersionSnpEFF`, `getVersionStrelka`, `getVersionVCFtools`, `getVersionVEP`
-   [#595](https://github.com/SciLifeLab/Sarek/pull/595) - New Processes `GetVersionBCFtools`, `GetVersionGATK`, `GetVersionManta`, `GetVersionSnpEFF`, `GetVersionStrelka`, `GetVersionVCFtools`, `GetVersionVEP`
-   [#595](https://github.com/SciLifeLab/Sarek/pull/595) - new Python script `bin/scrape_tool_versions.py` inspired by @ewels and @apeltzer
-   [#595](https://github.com/SciLifeLab/Sarek/pull/595) - New QC Process `RunVcftools`
-   [#596](https://github.com/SciLifeLab/Sarek/pull/596) - New profile for BinAC cluster
-   [#597](https://github.com/SciLifeLab/Sarek/pull/597) - New function `sarek_ascii()` in `SarekUtils`
-   [#599](https://github.com/SciLifeLab/Sarek/pull/599), [#602](https://github.com/SciLifeLab/Sarek/pull/602) - New Process `CompressVCF`
-   [#601](https://github.com/SciLifeLab/Sarek/pull/601), [#603](https://github.com/SciLifeLab/Sarek/pull/603) - Container for GATK4
-   [#606](https://github.com/SciLifeLab/Sarek/pull/606) - Add test data as a submodule from [`Sarek-data`](https://github.com/SciLifeLab/Sarek-data)
-   [#608](https://github.com/SciLifeLab/Sarek/pull/608) - Add documentation on how to install Nextflow on `bianca`

### `Changed`

-   [#557](https://github.com/SciLifeLab/Sarek/pull/557), [#583](https://github.com/SciLifeLab/Sarek/pull/583), [#585](https://github.com/SciLifeLab/Sarek/pull/585), [#588](https://github.com/SciLifeLab/Sarek/pull/588) - Update help
-   [#560](https://github.com/SciLifeLab/Sarek/pull/560) - GitHub langage for the repository is now `Nextflow`
-   [#561](https://github.com/SciLifeLab/Sarek/pull/561) - `do_all.sh` build only containers for one genome reference (default `GRCh38`) only
-   [#571](https://github.com/SciLifeLab/Sarek/pull/571) - Only one container for all QC tools
-   [#582](https://github.com/SciLifeLab/Sarek/pull/582), [#587](https://github.com/SciLifeLab/Sarek/pull/587) - Update figures
-   [#595](https://github.com/SciLifeLab/Sarek/pull/595) - Function `defineDirectoryMap()` is now part of `SarekUtils`
-   [#595](https://github.com/SciLifeLab/Sarek/pull/595) - Process `GenerateMultiQCconfig` replace by function `createMultiQCconfig()`
-   [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Move `checkFileExtension()`, `checkParameterExistence()`, `checkParameterList()`, `checkReferenceMap()`, `checkRefExistence()`, `extractBams()`, `extractGenders()`, `returnFile()`, `returnStatus()` and `returnTSV()` functions to `SarekUtils`
-   [#597](https://github.com/SciLifeLab/Sarek/pull/597) - `extractBams()` now takes an extra parameter.
-   [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Replace depreciated operator `phase` by `join`.
-   [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Reduce data footprint for Process `CreateRecalibrationTable`
-   [#599](https://github.com/SciLifeLab/Sarek/pull/599) - Merge is tested with `ANNOTATEALL`
-   [#604](https://github.com/SciLifeLab/Sarek/pull/604) - Synching `GRCh38` `wgs_calling_regions` bedfiles
-   [#607](https://github.com/SciLifeLab/Sarek/pull/607) - Update to GATK4
-   [#607](https://github.com/SciLifeLab/Sarek/pull/607) - One container approach
-   [#608](https://github.com/SciLifeLab/Sarek/pull/608) - Update Nextflow required version
-   [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Update CHANGELOG
-   [#617](https://github.com/SciLifeLab/Sarek/pull/617) - Replace depreciated $name syntax with withName

### `Fixed`

-   [#560](https://github.com/SciLifeLab/Sarek/pull/560) - Display message for `repository` and `containerPath`
-   [#566](https://github.com/SciLifeLab/Sarek/pull/566) - `slurmDownload` profile
-   [#579](https://github.com/SciLifeLab/Sarek/pull/579), [#584](https://github.com/SciLifeLab/Sarek/pull/584) - `Manta` output reorganized after modification for `Strelka Best Practices` process
-   [#585](https://github.com/SciLifeLab/Sarek/pull/583) - Trace file is plain txt
-   [#590](https://github.com/SciLifeLab/Sarek/pull/590), [#593](https://github.com/SciLifeLab/Sarek/pull/593) - Fix Singularity installation in Travis CI testing
-   [#598](https://github.com/SciLifeLab/Sarek/pull/598), [#601](https://github.com/SciLifeLab/Sarek/pull/601) - Fixes for  Python script `selectROI.py` to work with CLC viewer

### `Removed`

-   [#607](https://github.com/SciLifeLab/Sarek/pull/607) - Remove Mutect1

## [2.0.0] - 2018-03-23

### `Added`

-   basic wrapper script
-   Abstract, posters and figures
-   ROI selector and FreeBayes sanitizer scripts
-   New logo and icon for the project
-   check for existing tumor/normal channel
-   `SarekUtils` with `checkParams()`, `checkParameterList()`, `checkParameterExistence()` and `isAllowedParams()` functions
-   some `runOptions` for `docker` (prevent some user right problem)
-   This `CHANGELOG`

### `Changed`

-   `CAW` is now `Sarek`
-   Dissect Workflow in 5 new scripts: `annotate.nf`, `main.nf`, `germlineVC.nf`, `runMultiQC.nf` and `somaticVC.nf`
-   `report.html`, `timeline.html` and `trace.html` are generated in `Reports/`
-   `--version` is now used to define the workflow version
-   most params are now defined in the base.config file instead of in the scripts
-   update RELEASE_CHECKLIST.md
-   `checkParams()`, `checkParameterList()`, `checkParameterExistence()` and `isAllowedParams()` in script functions are now called within `SarekUtils`
-   `nf_required_version` is now `params.nfRequiredVersion`
-   in `buildReferences.nf` script, channels now begin by `ch_`, and files by `f_`
-   use `PublishDir mode: 'link'` instead of `copy`
-   `directoryMap` now contains `params.outDir`
-   [#539](https://github.com/SciLifeLab/Sarek/issues/539) - use Nextflow support of scratch
-   reordered Travis CI tests
-   update documentation
-   `MultiQC` version in container from v`1.4` to v`1.5`
-   `vepgrch37` container base image from `release_90.6` to `release_92`
-   `vepgrch38` container base image from `release_90.6` to `release_92`
-   `VEP` version in containers from v`90` to v`91`
-   `nucleotidesPerSecond` is now `params.nucleotidesPerSecond`
-   default `params.tag` is now `latest` instead of current version, so --tag needs to be specified with the right version to be sure of using the `containers` corresponding

### `Deprecated`

-   `standard` profile
-   `uppmax-localhost.config` file

### `Removed`

-   `scripts/skeleton_batch.sh`
-   old data and tsv files
-   UPPMAX directories from containers
-   `--step` in `annotate.nf`, `germlineVC.nf` and `somatic.nf`
-   some `runOptions` for Singularity (binding not needed anymore on UPPMAX)
-   `download` profile

### `Fixed`

-   [#533](https://github.com/SciLifeLab/Sarek/issues/533) - Replace `VEP` `--pick` option by `--per_gene`
-   [#530](https://github.com/SciLifeLab/Sarek/issues/530) - use `$PWD` for default `outDir`

## [1.2.5] - 2018-01-18

### `Added`

-   Zenodo for DOI
-   Delivery README
-   Document use of the `--sampleDir` option
-   Contributing Guidelines
-   Issue Templates
-   Release Checklist
-   `--outDir`
-   `awsbatch` profile
-   `aws-batch.config` config file
-   `--noBAMQC` params (failing sometimes on Bianca)

### `Changed`

-   Update `Nextflow` to `0.26.0` (new fancy report + AWS Batch)
-   Extra time on Travis CI testing
-   Replace `bundleDir` by `params.genome_base`
-   Update `MultiQC` to `1.3` (MEGAQC FTW)
-   Move and rename some test files

### `Fixed`

-   Version of COSMIC GRCh37 v83
-   Write an error message when `--sampleDir` does not find any FASTQ files
-   `base.config` for ConcatVCF process
-   File specification for recalibrationReport in RecalibrateBam process (got error on AWS Batch)

## [1.2.4] - 2017-10-27

### `Fixed`

-   [#488](https://github.com/SciLifeLab/Sarek/issues/488) - Better CPU requirements for `ConcatVCF`
-   [#489](https://github.com/SciLifeLab/Sarek/issues/489) - Exception handling for `ASCAT`
-   [#490](https://github.com/SciLifeLab/Sarek/issues/490) - CPU requirements for `runSingleStrelka` and `runSingleManta`

## [1.2.3] - 2017-10-18

### `Fixed`

-   [#475](https://github.com/SciLifeLab/Sarek/issues/475) - 16 cpus for local executor
-   [#357](https://github.com/SciLifeLab/Sarek/issues/357) - `ASCAT` works for GRCh38
-   [#471](https://github.com/SciLifeLab/Sarek/issues/471) - Running `Singularity` on `/scratch`
-   [#480](https://github.com/SciLifeLab/Sarek/issues/480) - No `tsv` file needed for step `annotate`

## [1.2.2] - 2017-10-06

### `Fixed`

-   [#479](https://github.com/SciLifeLab/Sarek/issues/479) - Typo in `uppmax-localhost.config`

## [1.2.1] - 2017-10-06

### `Changed`

-   `runascat` and `runconvertallelecounts` containers are now replaced by `r-base`
-   `willmclaren/ensembl-vep:release_90.5` is now base for `vepgrch37` and `vepgrch38`

### `Removed`

-   `vep` container
-   `strelka_config.ini` file

### `Fixed`

-   [#471](https://github.com/SciLifeLab/Sarek/issues/471) - Running `Singularity` on /scratch
-   [#472](https://github.com/SciLifeLab/Sarek/issues/472) - Update function to check Nextflow version
-   [#473](https://github.com/SciLifeLab/Sarek/issues/473) - Remove `returnMin()` function

## [1.2.0] - 2017-10-02

### `Changed`

-   Fix version for Manuscript

## [1.1] - 2017-09-15

### `Added`

-   Singularity possibilities

### `Changed`

-   Reports made by default
-   Intervals file can be a bed file
-   Normal sample preprocessing + HaplotypeCaller is possible
-   Better Travis CI tests

### `Fixed`

-   Memory requirements

## [1.0] - 2017-02-16

### `Added`

-   Docker possibilities

## [0.9] - 2016-11-16

## [0.8] - 2016-11-16

## [0.1] - 2016-04-05
