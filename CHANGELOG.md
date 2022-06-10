# nf-core/sarek: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.7.2](https://github.com/nf-core/sarek/releases/tag/2.7.2) - Áhkká

Áhkká is one of the massifs just outside of the Sarek National Park.

### Fixed

- [#566](https://github.com/nf-core/sarek/pull/566) - Fix caching bug affecting a variable number of `MapReads` jobs due to non-deterministic state of `statusMap` during caching evaluation


## [2.7.1](https://github.com/nf-core/sarek/releases/tag/2.7.1) - Pårtejekna

Pårtejekna is one of glaciers of the Pårte Massif.

### Added

- [#353](https://github.com/nf-core/sarek/pull/353) - Add support for task retries with exit code 247 (exhibited by `Picard MarkDuplicates`)
- [#354](https://github.com/nf-core/sarek/pull/354) - Add tumor only mode for `Mutect2` and `MSIsensor`
- [#356](https://github.com/nf-core/sarek/pull/356) - Add `--cf_contamination_adjustment` params to adjust contamination with `Control-FREEC`
- [#372](https://github.com/nf-core/sarek/pull/372) - Add `--cf_contamination` params to specify contamination value with `Control-FREEC`

### Changed

- [#373](https://github.com/nf-core/sarek/pull/373) - Sync `TEMPLATE` with `tools` 1.14
- [#376](https://github.com/nf-core/sarek/pull/376) - Better logo on Github dark Mode
- [#387](https://github.com/nf-core/sarek/pull/387) - Fix tables for TSV file content

### Fixed

- [#375](https://github.com/nf-core/sarek/pull/375), [#381](https://github.com/nf-core/sarek/pull/381), [#382](https://github.com/nf-core/sarek/pull/382), [#385](https://github.com/nf-core/sarek/pull/385) - Fix bugs due to `TEMPLATE` sync from [#373](https://github.com/nf-core/sarek/pull/373)
- [#378](https://github.com/nf-core/sarek/pull/378) - Fix `Spark` related issue due to `Docker` settings in `nextflow.config`
- [#378](https://github.com/nf-core/sarek/pull/378) - Fix `Spark` related issue due to `Docker` settings in `nextflow.config`

### Deprecated

### Removed

- [#368](https://github.com/nf-core/sarek/pull/368) - Remove social preview image to use GitHub OpenGraph

## [2.7](https://github.com/nf-core/sarek/releases/tag/2.7) - Pårte

Pårte is one of the main massif in the Sarek National Park.

### Added

- [#145](https://github.com/nf-core/sarek/pull/145) - Add `UMI annotation and consensus` functionality to `Sarek`
- [#230](https://github.com/nf-core/sarek/pull/230) - Add `ignore_soft_clipped_bases` option for `GATK Mutect2` [#218](https://github.com/nf-core/sarek/issues/218)
- [#253](https://github.com/nf-core/sarek/pull/253) - Add `UMI` `CI` testing
- [#262](https://github.com/nf-core/sarek/pull/262) - Add `nextflow_schema.json`
- [#237](https://github.com/nf-core/sarek/pull/237), [#282](https://github.com/nf-core/sarek/pull/282) - Add `--aligner` to choose between `bwa` and `bwa-mem2`
- [#294](https://github.com/nf-core/sarek/pull/294) - Add `Troubleshooting` section to `docs/usage.md`
- [#302](https://github.com/nf-core/sarek/pull/302), [#304](https://github.com/nf-core/sarek/pull/304) - Add WES and tumor-only mode for `Control-FREEC`

### Changed

- [#253](https://github.com/nf-core/sarek/pull/253), [#255](https://github.com/nf-core/sarek/pull/255), [#326](https://github.com/nf-core/sarek/pull/326), [#329](https://github.com/nf-core/sarek/pull/329) - Update docs
- [#260](https://github.com/nf-core/sarek/pull/260), [#262](https://github.com/nf-core/sarek/pull/262), [#278](https://github.com/nf-core/sarek/pull/278), [#322](https://github.com/nf-core/sarek/pull/322) - Sync with `TEMPLATE` updated from [nf-core/tools](https://github.com/nf-core/tools) [`1.10.2`](https://github.com/nf-core/tools/releases/tag/1.10.2)
- [#262](https://github.com/nf-core/sarek/pull/262) - Update issue templates to fit the recommended community standards
- [#278](https://github.com/nf-core/sarek/pull/278), [#322](https://github.com/nf-core/sarek/pull/322) - Refactor docs
- [#284](https://github.com/nf-core/sarek/pull/284) - Update F1000Research publication to version 2
- [#284](https://github.com/nf-core/sarek/pull/284) - Update Scilifelab logo
- [#317](https://github.com/nf-core/sarek/pull/317) - Update `README.md` (Add: QBiC + Friederike/Gisela)
- [#320](https://github.com/nf-core/sarek/pull/278) - Set `MarkDuplicates MAX_RECORDS_IN_RAM` to default value

### Fixed

- [#229](https://github.com/nf-core/sarek/pull/229) - Fix `Control-FREEC` restart issue [#225](https://github.com/nf-core/sarek/issues/225)
- [#236](https://github.com/nf-core/sarek/pull/236) - Fix `GATK Mutect2` typo issue [#227](https://github.com/nf-core/sarek/issues/227)
- [#271](https://github.com/nf-core/sarek/pull/271) - Fix `ConcatVCF_Mutect2` `SIGPIPE` issue [#268](https://github.com/nf-core/sarek/issues/268)
- [#272](https://github.com/nf-core/sarek/pull/272) - Fix annotation `--tools merge` issue
- [#279](https://github.com/nf-core/sarek/pull/279) - Fix issue with `--step prepare_recalibration` [#267](https://github.com/nf-core/sarek/issues/267)
- [#280](https://github.com/nf-core/sarek/pull/280) - Use HTML codes instead of `<` and `>` in docs
- [#288](https://github.com/nf-core/sarek/pull/288) - Fix `test_annotation` profile
- [#289](https://github.com/nf-core/sarek/pull/289) - Random string added to `extractFastqFromDir` to avoid name collition
- [#290](https://github.com/nf-core/sarek/pull/290), [#323](https://github.com/nf-core/sarek/pull/323) - Faster solving of `Conda` environment
- [#293](https://github.com/nf-core/sarek/pull/293) - Fix typo issue when printing infos [#292](https://github.com/nf-core/sarek/issues/292)
- [#309](https://github.com/nf-core/sarek/pull/309) - Fixed concatenation of many VCF files
- [#310](https://github.com/nf-core/sarek/pull/310) - Fix Github Actions not running after November 16, 2020 (deprecated Github Actions API [#739](https://github.com/nf-core/tools/issues/739)
- [#329](https://github.com/nf-core/sarek/pull/329) - Simplify `Control-FREEC` usage
- [#331](https://github.com/nf-core/sarek/pull/331) - Replace `spread` operator by `combine` to remove `Nextflow` deprecation warning

### Removed

- [#234](https://github.com/nf-core/sarek/pull/243) - Removing obsolete script [#92](https://github.com/nf-core/sarek/issues/92)
- [#262](https://github.com/nf-core/sarek/pull/262) - Removing deprecated params: `annotateTools`, `annotateVCF`, `cadd_InDels`, `cadd_InDels_tbi`, `cadd_WG_SNVs`, `cadd_WG_SNVs_tbi`, `maxMultiqcEmailFileSize`, `noGVCF`, `noReports`, `noStrelkaBP`, `nucleotidesPerSecond`, `publishDirMode`, `sample`, `sampleDir`, `saveGenomeIndex`, `skipQC`, `snpEff_cache`, `targetBed`
- [#262](https://github.com/nf-core/sarek/pull/262) - Removing warning message about deprecated and obsolete params
- [#324](https://github.com/nf-core/sarek/pull/324) - `--no_gatk_spark` is now removed, use `--use_gatk_spark` instead
- [#324](https://github.com/nf-core/sarek/pull/324) - `--no_gvcf` is now removed, use `--generate_gvcf` instead

## [2.6.1](https://github.com/nf-core/sarek/releases/tag/2.6.1) - Gådokgaskatjåhkkå

Gådokgaskatjåhkkå is the highest peak in the Piellorieppe massif.

### Changed

- [#208](https://github.com/nf-core/sarek/pull/208) - Merge changes from the release PR
- [#208](https://github.com/nf-core/sarek/pull/208) - Bump version to `3.0dev`
- [#214](https://github.com/nf-core/sarek/pull/214) - Update `GATK` from `4.1.6.0` to `4.1.7.0`
- [#219](https://github.com/nf-core/sarek/pull/219) - Added `awsfulltest.yml` GitHub Actions workflow
- [#222](https://github.com/nf-core/sarek/pull/222) - Bump version to `2.6.1` and minor release
- [#223](https://github.com/nf-core/sarek/pull/223) - Apply comments from the release PR

### Fixed

- [#211](https://github.com/nf-core/sarek/pull/211) - Extend timeout for pushing to DockerHub for VEP containers
- [#212](https://github.com/nf-core/sarek/pull/212) - No AWS test on forks
- [#214](https://github.com/nf-core/sarek/pull/214) - Fix channels collision between `Freebayes` and `GATK Mutect2` [#200](https://github.com/nf-core/sarek/issues/200)
- [#214](https://github.com/nf-core/sarek/pull/214) - Fix warning Invalid tag value for `CreateIntervalBeds` [#209](https://github.com/nf-core/sarek/issues/209)
- [#214](https://github.com/nf-core/sarek/pull/214) - Fix `GATK Mutect2` issue [#210](https://github.com/nf-core/sarek/issues/210)
- [#219](https://github.com/nf-core/sarek/pull/219) - Updated `awstest.yml` GitHub actions workflow
- [#221](https://github.com/nf-core/sarek/pull/221) - Fix issue with `tmp_dir` in `BaseRecalibrator` process

## [2.6](https://github.com/nf-core/sarek/releases/tag/2.6) - Piellorieppe

Piellorieppe is one of the main massif in the Sarek National Park.

### Added

- [#76](https://github.com/nf-core/sarek/pull/76) - Add `GATK Spark` possibilities to Sarek
- [#87](https://github.com/nf-core/sarek/pull/87) - Add `GATK BaseRecalibrator` plot to `MultiQC` report
- [#115](https://github.com/nf-core/sarek/pull/115) - Add [@szilvajuhos](https://github.com/szilvajuhos) abstract for ESHG2020
- [#117](https://github.com/nf-core/sarek/pull/117) - Add `Trim Galore` possibilities to Sarek
- [#141](https://github.com/nf-core/sarek/pull/141) - Add containers for `WBcel235`
- [#150](https://github.com/nf-core/sarek/pull/150), [#151](https://github.com/nf-core/sarek/pull/151), [#154](https://github.com/nf-core/sarek/pull/154) - Add AWS mega test GitHub Actions
- [#153](https://github.com/nf-core/sarek/pull/153) - Add `CNVkit` possibilities to Sarek
- [#158](https://github.com/nf-core/sarek/pull/158) - Added `ggplot2` version `3.3.0`
- [#163](https://github.com/nf-core/sarek/pull/163) - Add [MSIsensor](https://github.com/ding-lab/msisensor) in tools and container
- [#164](https://github.com/nf-core/sarek/pull/164) - Add `--no_gatk_spark` params and tests
- [#167](https://github.com/nf-core/sarek/pull/167) - Add `--markdup_java_options` documentation
- [#169](https://github.com/nf-core/sarek/pull/169) - Add `RELEASE_CHECKLIST.md` document
- [#174](https://github.com/nf-core/sarek/pull/174) - Add `variant_calling.md` documentation
- [#175](https://github.com/nf-core/sarek/pull/175) - Add `Sentieon` documentation
- [#176](https://github.com/nf-core/sarek/pull/176) - Add empty `custom` genome in `genomes.config` to allow genomes that are not in `AWS iGenomes`
- [#179](https://github.com/nf-core/sarek/pull/179), [#201](https://github.com/nf-core/sarek/pull/201) - Add `FreeBayes` germline variant calling
- [#180](https://github.com/nf-core/sarek/pull/180) - Now saving Mapped BAMs (and creating TSV) in minimal setting
- [#182](https://github.com/nf-core/sarek/pull/182) - Add possibility to run `HaplotypeCaller` without `dbsnp` so it can be used to actually generate vcfs to build a set of known sites (cf [gatkforums](https://gatkforums.broadinstitute.org/gatk/discussion/1247/what-should-i-use-as-known-variants-sites-for-running-tool-x))
- [#195](https://github.com/nf-core/sarek/pull/195) - Now creating TSV for duplicates marked BAMs in minimal setting
- [#195](https://github.com/nf-core/sarek/pull/195), [#202](https://github.com/nf-core/sarek/pull/202) - Add `--save_bam_mapped` params to save mapped BAMs.
- [#197](https://github.com/nf-core/sarek/pull/197) - Add step `prepare_recalibration` to allow restart from DuplicatesMarked BAMs
- [#204](https://github.com/nf-core/sarek/pull/204) - Add step `Control-FREEC` to allow restart from pileup files
- [#205](https://github.com/nf-core/sarek/pull/205) - Add `--skip_markduplicates` to allow skipping the `MarkDuplicates` process

### Changed

- [#76](https://github.com/nf-core/sarek/pull/76) - Use `MarkDuplicatesSpark` instead of `MarkDuplicates`
- [#76](https://github.com/nf-core/sarek/pull/76) - Use `gatk4-spark` instead of `gatk4` in `environment.yml`
- [#80](https://github.com/nf-core/sarek/pull/80) - Re-bump `dev` branch
- [#85](https://github.com/nf-core/sarek/pull/85) - Use new merged vcf files for known indels to simplify setting up channel
- [#104](https://github.com/nf-core/sarek/pull/104) - Update Figure 1
- [#107](https://github.com/nf-core/sarek/pull/107) - Switch params to snake_case
- [#109](https://github.com/nf-core/sarek/pull/109) - Update publication with F1000Research preprint
- [#113](https://github.com/nf-core/sarek/pull/113) - Move social preview image
- [#120](https://github.com/nf-core/sarek/pull/120) - Sync TEMPLATE
- [#121](https://github.com/nf-core/sarek/pull/121) - Update `MultiQC` to `1.8`
- [#126](https://github.com/nf-core/sarek/pull/126), [#131](https://github.com/nf-core/sarek/pull/131) - Update docs
- [#131](https://github.com/nf-core/sarek/pull/131) - Use `nfcore/base:1.9` as base for containers
- [#131](https://github.com/nf-core/sarek/pull/131) - Update `Control-FREEC` to `11.5`
- [#131](https://github.com/nf-core/sarek/pull/131) - Update `FastQC` to `0.11.9`
- [#131](https://github.com/nf-core/sarek/pull/131) - Update `FreeBayes` to `1.3.2`
- [#131](https://github.com/nf-core/sarek/pull/131) - Update `Manta` to `1.6.0`
- [#131](https://github.com/nf-core/sarek/pull/131) - Update `Qualimap` to `2.2.2d`
- [#131](https://github.com/nf-core/sarek/pull/131) - Update `VEP` to `99.2`
- [#141](https://github.com/nf-core/sarek/pull/141) - Update `snpEff` cache version from `75` to `87` for `GRCh37`
- [#141](https://github.com/nf-core/sarek/pull/141) - Update `snpEff` cache version from `86` to `92` for `GRCh38`
- [#141](https://github.com/nf-core/sarek/pull/141) - Update `VEP` databases to `99`
- [#143](https://github.com/nf-core/sarek/pull/143) - Revert `snpEff` cache version to `75` for `GRCh37`
- [#143](https://github.com/nf-core/sarek/pull/143) - Revert `snpEff` cache version to `86` for `GRCh38`
- [#152](https://github.com/nf-core/sarek/pull/152), [#158](https://github.com/nf-core/sarek/pull/158), [#164](https://github.com/nf-core/sarek/pull/164), [#174](https://github.com/nf-core/sarek/pull/174), [#194](https://github.com/nf-core/sarek/pull/194), [#198](https://github.com/nf-core/sarek/pull/198), [#204](https://github.com/nf-core/sarek/pull/204) - Update docs
- [#164](https://github.com/nf-core/sarek/pull/164) - Update `gatk4-spark` from `4.1.4.1` to `4.1.6.0`
- [#180](https://github.com/nf-core/sarek/pull/180), [#195](https://github.com/nf-core/sarek/pull/195) - Improve minimal setting
- [#183](https://github.com/nf-core/sarek/pull/183), [#204](https://github.com/nf-core/sarek/pull/204) - Update `input.md` documentation
- [#197](https://github.com/nf-core/sarek/pull/197) - Output directory `DuplicateMarked` is now replaced by `DuplicatesMarked`
- [#204](https://github.com/nf-core/sarek/pull/204) - Output directory `controlFREEC` is now replaced by `Control-FREEC`

### Fixed

- [#83](https://github.com/nf-core/sarek/pull/83) - Fix some typos in `docs/input.md`
- [#107](https://github.com/nf-core/sarek/pull/107) - Fix linting
- [#110](https://github.com/nf-core/sarek/pull/110) - Fix `snpEff` report issue cf [#106](https://github.com/nf-core/sarek/issues/106)
- [#126](https://github.com/nf-core/sarek/pull/126) - Fix `iGenomes` paths
- [#127](https://github.com/nf-core/sarek/pull/127), [#128](https://github.com/nf-core/sarek/pull/128) - Fix `ASCAT`
- [#129](https://github.com/nf-core/sarek/pull/129) - Fix issue with Channel `channel ch_software_versions_yaml`
- [#129](https://github.com/nf-core/sarek/pull/129) - Apply @drpatelh fix for `mardown_to_html.py` compatibility with Python 2
- [#129](https://github.com/nf-core/sarek/pull/129) - Removed `Python` `3.7.3` from conda environment due to incompatibility
- [#129](https://github.com/nf-core/sarek/pull/129) - Change ascii characters that were not supported from the `output.md` docs
- [#140](https://github.com/nf-core/sarek/pull/140) - Fix extra T/N combinations for `ASCAT` cf [#136](https://github.com/nf-core/sarek/issues/136)
- [#141](https://github.com/nf-core/sarek/pull/141) - Fix `download_cache.nf` script to download cache for `snpEff` and `VEP`
- [#143](https://github.com/nf-core/sarek/pull/143) - Fix annotation CI testing with `snpEff` and `VEP`
- [#144](https://github.com/nf-core/sarek/pull/144) - Fix CircleCI for building `VEP` containers
- [#146](https://github.com/nf-core/sarek/pull/146) - Fix `--no_intervals` for `GATK Mutect2` cf [#135](https://github.com/nf-core/sarek/issues/135)
- [#156](https://github.com/nf-core/sarek/pull/156) - Fix typos
- [#156](https://github.com/nf-core/sarek/pull/156) - Fix issues with `dbsnp` files while using only `Sention` tools
- [#158](https://github.com/nf-core/sarek/pull/158) - Fix typo with `params.snpeff_cache` to decide containers for `snpEff`
- [#164](https://github.com/nf-core/sarek/pull/164) - Fix issues when running with `Sentieon`
- [#164](https://github.com/nf-core/sarek/pull/164) - Add more VCFs to annotation
- [#167](https://github.com/nf-core/sarek/pull/167) - Add `--markdup_java_options` documentation to fix [#166](https://github.com/nf-core/sarek/issues/166)
- [#178](https://github.com/nf-core/sarek/pull/178) - Fix `Sentieon` variant calling, now using deduped bam files
- [#188](https://github.com/nf-core/sarek/pull/188) - Fix input/output channels for process `IndexBamFile` to match actual files in the `mapped.tsv` files
- [#189](https://github.com/nf-core/sarek/pull/189) - Fix `no_intervals` for process `HaplotypeCaller` (the file just need to actually exists...)
- [#197](https://github.com/nf-core/sarek/pull/197) - Fix issue with `--step recalibrate`
- [#197](https://github.com/nf-core/sarek/pull/197) - Fix typo in output directory `DuplicateMarked` -> `DuplicatesMarked`

### Deprecated

- [#107](https://github.com/nf-core/sarek/pull/107) - `--annotateTools` is now deprecated, use `--annotate_tools` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--cadd_InDels` is now deprecated, use `--cadd_indels` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--cadd_InDels_tbi` is now deprecated, use `--cadd_indels_tbi` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--cadd_WG_SNVs` is now deprecated, use `--cadd_wg_snvs` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--cadd_WG_SNVs_tbi` is now deprecated, use `--cadd_wg_snvs_tbi` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--maxMultiqcEmailFileSize` is now deprecated, use `--max_multiqc_email_size` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--noGVCF` is now deprecated, use `--no_gvcf` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--noStrelkaBP` is now deprecated, use `--no_strelka_bp` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--nucleotidesPerSecond` is now deprecated, use `--nucleotides_per_second` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--publishDirMode` is now deprecated, use `--publish_dir_mode` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--saveGenomeIndex` is now deprecated, use `--save_reference` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--skipQC` is now deprecated, use `--skip_qc` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--snpEff_cache` is now deprecated, use `--snpeff_cache` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--targetBed` is now deprecated, use `--target_bed` instead

### Removed

- [#107](https://github.com/nf-core/sarek/pull/107) - `--acLociGC` is now removed, use `--ac_loci_gc` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--acLoci` is now removed, use `--ac_loci` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--bwaIndex` is now removed, use `--bwa` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--chrDir` is now removed, use `--chr_dir` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--chrLength` is now removed, use `--chr_length` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--dbsnpIndex` is now removed, use `--dbsnp_index` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--fastaFai` is now removed, use `--fasta_fai` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--genomeDict` is now removed, use `--dict` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--genomeFile` is now removed, use `--fasta` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--genomeIndex` is now removed, use `--fasta_fai` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--germlineResourceIndex` is now removed, use `--germline_resource_index` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--germlineResource` is now removed, use `--germline_resource` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--igenomesIgnore` is now removed, use `--igenomes_ignore` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--knownIndelsIndex` is now removed, use `--known_indels_index` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--knownIndels` is now removed, use `--known_indels` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--singleCPUMem` is now removed, use `--single_cpu_mem` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--snpeffDb` is now removed, use `--snpeff_db` instead
- [#107](https://github.com/nf-core/sarek/pull/107) - `--vepCacheVersion` is now removed, use `--vep_cache_version` instead
- [#152](https://github.com/nf-core/sarek/pull/152) - Removed `Jenkinsfile`
- [#169](https://github.com/nf-core/sarek/pull/169) - Removed omicX from README
- [#181](https://github.com/nf-core/sarek/pull/181) - Remove duplicate code in `nextflow.config`

## [2.5.2](https://github.com/nf-core/sarek/releases/tag/2.5.2) - Jåkkåtjkaskajekna

Jåkkåtjkaskajekna is one of the two glaciers of the Ålkatj Massif.

### Added

- [#45](https://github.com/nf-core/sarek/pull/45) - Include Workflow figure in `README.md`
- [#46](https://github.com/nf-core/sarek/pull/46) - Add location to abstracts
- [#52](https://github.com/nf-core/sarek/pull/52) - Add support for mouse data `GRCm38`
- [#60](https://github.com/nf-core/sarek/pull/60) - Add `no_intervals` params
- [#60](https://github.com/nf-core/sarek/pull/60) - Add automatic generation of `intervals` file with `BuildIntervals` process
- [#60](https://github.com/nf-core/sarek/pull/60) - Add minimal support for minimal genome (only `fasta`, or `fasta` + `knownIndels`)
- [#60](https://github.com/nf-core/sarek/pull/60) - Add new processes (`IndexBamFile`, `IndexBamRecal`) to deal with optional usage of interval files and minimal genome
- [#60](https://github.com/nf-core/sarek/pull/60) - Add tests for minimal genome usage
- [#60](https://github.com/nf-core/sarek/pull/60) - Add new minimal genomes (`TAIR10`, `EB2`, `UMD3.1`, `bosTau8`, `WBcel235`, `ce10`, `CanFam3.1`, `canFam3`, `GRCz10`, `danRer10`, `BDGP6`, `dm6`, `EquCab2`, `equCab2`, `EB1`, `Galgal4`, `galGal4`, `Gm01`, `hg38`, `hg19`, `Mmul_1`, `mm10`, `IRGSP-1.0`, `CHIMP2.1.4`, `panTro4`, `Rnor_6.0`, `rn6`, `R64-1-1`, `sacCer3`, `EF2`, `Sbi1`, `Sscrofa10.2`, `susScr3`, `AGPv3`) to `igenomes.config`
- [#61](https://github.com/nf-core/sarek/pull/61) - Add params `split_fastq`
- [#61](https://github.com/nf-core/sarek/pull/61) - Add test `SPLITFASTQ`
- [#66](https://github.com/nf-core/sarek/pull/66) - Add `Sentieon` possibilities to Sarek

### Changed

- [#54](https://github.com/nf-core/sarek/pull/54) - Bump version to `2.5.2dev`
- [#60](https://github.com/nf-core/sarek/pull/60) - Some process (`BaseRecalibrator`, `ApplyBQSR`, `Mpileup`) have now optional usage of interval files
- [#60](https://github.com/nf-core/sarek/pull/60) - Update documentation
- [#71](https://github.com/nf-core/sarek/pull/71) - Update `README`
- [#71](https://github.com/nf-core/sarek/pull/71) - Update `CHANGELOG`
- [#74](https://github.com/nf-core/sarek/pull/74) - Update docs
- [#74](https://github.com/nf-core/sarek/pull/74) - Improve CI tests (both Jenkins and GitHub actions tests)
- [#74](https://github.com/nf-core/sarek/pull/74) - Move all CI from `ci-extra.yml` to `ci.yml`

### Removed

- [#46](https://github.com/nf-core/sarek/pull/46) - Remove mention of old `build.nf` script which was included in `main.nf`
- [#74](https://github.com/nf-core/sarek/pull/74) - Remove `download_image.sh` and `run_tests.sh` scripts
- [#76](https://github.com/nf-core/sarek/pull/76) - Remove `runOptions = "-u \$(id -u):\$(id -g)"` in `nextflow.config` to enable `Spark` possibilities

### Fixed

- [#40](https://github.com/nf-core/sarek/pull/40) - Fix issue with `publishDirMode` within `test` profile
- [#42](https://github.com/nf-core/sarek/pull/42) - Fix typos, and minor updates in `README.md`
- [#43](https://github.com/nf-core/sarek/pull/43) - Fix automated `VEP` builds with circleCI
- [#54](https://github.com/nf-core/sarek/pull/54) - Apply fixes from release `2.5.1`
- [#58](https://github.com/nf-core/sarek/pull/58) - Fix issue with `.interval_list` file from the `GATK` bundle [#56](https://github.com/nf-core/sarek/issues/56) that was not recognized in the `CreateIntervalsBed` process
- [#71](https://github.com/nf-core/sarek/pull/71) - Fix typos in `CHANGELOG`
- [#73](https://github.com/nf-core/sarek/pull/73) - Fix issue with label `memory_max` for `BaseRecalibrator` process [#72](https://github.com/nf-core/sarek/issues/72)

## [2.5.1](https://github.com/nf-core/sarek/releases/tag/2.5.1) - Årjep-Ålkatjjekna

Årjep-Ålkatjjekna is one of the two glaciers of the Ålkatj Massif.

### Added

- [#53](https://github.com/nf-core/sarek/pull/53) - Release `2.5.1`

### Fixed

- [#48](https://github.com/nf-core/sarek/issues/48) - Fix `singularity.autoMounts` issue.
- [#49](https://github.com/nf-core/sarek/issues/49) - Use correct tag for annotation containers.
- [#50](https://github.com/nf-core/sarek/issues/50) - Fix paths for scripts.

## [2.5](https://github.com/nf-core/sarek/releases/tag/2.5) - Ålkatj

Ålkatj is one of the main massif in the Sarek National Park.

Initial release of `nf-core/sarek`, created with the [nf-core](http://nf-co.re/) template.

### Added

- [#2](https://github.com/nf-core/sarek/pull/2) - Create `nf-core/sarek` `environment.yml` file
- [#2](https://github.com/nf-core/sarek/pull/2), [#3](https://github.com/nf-core/sarek/pull/3), [#4](https://github.com/nf-core/sarek/pull/4), [#5](https://github.com/nf-core/sarek/pull/5), [#7](https://github.com/nf-core/sarek/pull/7), [#9](https://github.com/nf-core/sarek/pull/9), [#10](https://github.com/nf-core/sarek/pull/10), [#11](https://github.com/nf-core/sarek/pull/11), [#12](https://github.com/nf-core/sarek/pull/12) - Add CI for `nf-core/sarek`
- [#3](https://github.com/nf-core/sarek/pull/3) - Add preprocessing to `nf-core/sarek`
- [#4](https://github.com/nf-core/sarek/pull/4) - Add variant calling to `nf-core/sarek` with `HaplotypeCaller`, and single mode `Manta` and `Strelka`
- [#5](https://github.com/nf-core/sarek/pull/5), [#34](https://github.com/nf-core/sarek/pull/34) - Add variant calling to `nf-core/sarek` with `Manta`, `Strelka`, `Strelka Best Practices`, `GATK Mutect2`, `FreeBayes`, `ASCAT`, `ControlFREEC`
- [#6](https://github.com/nf-core/sarek/pull/6) - Add default containers for annotation to `nf-core/sarek`
- [#7](https://github.com/nf-core/sarek/pull/7) - Add `MultiQC`
- [#7](https://github.com/nf-core/sarek/pull/7) - Add annotation
- [#7](https://github.com/nf-core/sarek/pull/7) - Add social preview image in `png` and `svg` format
- [#7](https://github.com/nf-core/sarek/pull/7), [#8](https://github.com/nf-core/sarek/pull/8), [#11](https://github.com/nf-core/sarek/pull/11), [#21](https://github.com/nf-core/sarek/pull/21) - Add helper script `run_tests.sh` to run different tests
- [#7](https://github.com/nf-core/sarek/pull/7), [#8](https://github.com/nf-core/sarek/pull/8), [#9](https://github.com/nf-core/sarek/pull/9) - Add automatic build of specific containers for annotation for `GRCh37`, `GRCh38` and `GRCm38` using `CircleCI`
- [#7](https://github.com/nf-core/sarek/pull/7), [#8](https://github.com/nf-core/sarek/pull/8), [#9](https://github.com/nf-core/sarek/pull/9), [#11](https://github.com/nf-core/sarek/pull/11) - Add helper script `build_reference.sh` to build small reference from [nf-core/test-datasets:sarek](https://github.com/nf-core/test-datasets/tree/sarek)
- [#7](https://github.com/nf-core/sarek/pull/7), [#9](https://github.com/nf-core/sarek/pull/9), [#11](https://github.com/nf-core/sarek/pull/11), [#12](https://github.com/nf-core/sarek/pull/12) - Add helper script `download_image.sh` to download containers for testing
- [#8](https://github.com/nf-core/sarek/pull/8) - Add test configuration for easier testing
- [#9](https://github.com/nf-core/sarek/pull/9), [#11](https://github.com/nf-core/sarek/pull/11) - Add scripts for `ASCAT`
- [#10](https://github.com/nf-core/sarek/pull/10) - Add `TIDDIT` to detect structural variants
- [#11](https://github.com/nf-core/sarek/pull/11) - Add automatic build of specific containers for annotation for `CanFam3.1` using `CircleCI`
- [#11](https://github.com/nf-core/sarek/pull/11), [#12](https://github.com/nf-core/sarek/pull/12) - Add posters and abstracts
- [#12](https://github.com/nf-core/sarek/pull/12) - Add helper script `make_snapshot.sh` to make an archive for usage on a secure cluster
- [#12](https://github.com/nf-core/sarek/pull/12) - Add helper scripts `filter_locifile.py` and `selectROI.py`
- [#12](https://github.com/nf-core/sarek/pull/12) - Use `label` for processes configuration
- [#13](https://github.com/nf-core/sarek/pull/13) - Add Citation documentation
- [#13](https://github.com/nf-core/sarek/pull/13) - Add `BamQC` process
- [#13](https://github.com/nf-core/sarek/pull/13) - Add `CompressVCFsnpEff` and `CompressVCFvep` processes
- [#18](https://github.com/nf-core/sarek/pull/18) - Add `--no-reports` option for tests + add snpEff,VEP,merge to MULTIPLE test
- [#18](https://github.com/nf-core/sarek/pull/18) - Add logo to `MultiQC` report
- [#18](https://github.com/nf-core/sarek/pull/18), [#29](https://github.com/nf-core/sarek/pull/29) - Add params `--skipQC` to skip specified QC tools
- [#18](https://github.com/nf-core/sarek/pull/18) - Add possibility to download other genome for `sareksnpeff` and `sarekvep` containers
- [#20](https://github.com/nf-core/sarek/pull/20) - Add `markdownlint` config file
- [#21](https://github.com/nf-core/sarek/pull/21) - Add tests for latest `Nextflow` version as well
- [#21](https://github.com/nf-core/sarek/pull/21) - Add `genomes.config` for genomes without `AWS iGenomes`
- [#24](https://github.com/nf-core/sarek/pull/24) - Added `GATK4 Mutect2` calling and filtering
- [#27](https://github.com/nf-core/sarek/pull/27), [#30](https://github.com/nf-core/sarek/pull/30) - Use Github actions for CI, linting and branch protection
- [#31](https://github.com/nf-core/sarek/pull/31) - Add `nf-core lint`
- [#31](https://github.com/nf-core/sarek/pull/31) - Add extra CI to `GitHub Actions` nf-core extra CI
- [#35](https://github.com/nf-core/sarek/pull/35) - Building indexes from [nf-core/test-datasets:sarek](https://github.com/nf-core/test-datasets/tree/sarek) for CI and small tests

### Changed

- [#1](https://github.com/nf-core/sarek/pull/1), [#2](https://github.com/nf-core/sarek/pull/2), [#3](https://github.com/nf-core/sarek/pull/3), [#4](https://github.com/nf-core/sarek/pull/4), [#5](https://github.com/nf-core/sarek/pull/5), [#6](https://github.com/nf-core/sarek/pull/6), [#7](https://github.com/nf-core/sarek/pull/7), [#8](https://github.com/nf-core/sarek/pull/8), [#9](https://github.com/nf-core/sarek/pull/9), [#10](https://github.com/nf-core/sarek/pull/10), [#11](https://github.com/nf-core/sarek/pull/11), [#12](https://github.com/nf-core/sarek/pull/12), [#18](https://github.com/nf-core/sarek/pull/18), [#20](https://github.com/nf-core/sarek/pull/20), [#21](https://github.com/nf-core/sarek/pull/21), [#23](https://github.com/nf-core/sarek/pull/23), [#29](https://github.com/nf-core/sarek/pull/29) - Update docs
- [#4](https://github.com/nf-core/sarek/pull/4) - Update `cancerit-allelecount` from `2.1.2` to `4.0.2`
- [#4](https://github.com/nf-core/sarek/pull/4) - Update `gatk4` from `4.1.1.0` to `4.1.2.0`
- [#7](https://github.com/nf-core/sarek/pull/7), [#23](https://github.com/nf-core/sarek/pull/23) - `--sampleDir` is now deprecated, use `--input` instead
- [#7](https://github.com/nf-core/sarek/pull/8), [#23](https://github.com/nf-core/sarek/pull/23) - `--annotateVCF` is now deprecated, use `--input` instead
- [#8](https://github.com/nf-core/sarek/pull/8), [#12](https://github.com/nf-core/sarek/pull/12) - Improve helper script `build.nf` for downloading and building reference files
- [#9](https://github.com/nf-core/sarek/pull/9) - `ApplyBQSR` is now parallelized
- [#9](https://github.com/nf-core/sarek/pull/9) - Fastq files are named following "${idRun}_R1.fastq.gz" in the `FastQC` output for easier reporting
- [#9](https://github.com/nf-core/sarek/pull/9) - Status is now a map with `idpatient`, `idsample` as keys (ie: `status = statusMap[idPatient, idSample]`)
- [#9](https://github.com/nf-core/sarek/pull/9) - Use `ensembl-vep` `95.2` instead of `96.0`
- [#11](https://github.com/nf-core/sarek/pull/11) - Summary HTML from `VEP` is now in the `Reports` directory
- [#12](https://github.com/nf-core/sarek/pull/12) - Update configuration files
- [#12](https://github.com/nf-core/sarek/pull/12) - Disable `Docker` in `singularity` profile
- [#12](https://github.com/nf-core/sarek/pull/12) - Disable `Singularity` in `docker` profile
- [#12](https://github.com/nf-core/sarek/pull/12) - Disable `Docker` and `Singularity` in `conda` profile
- [#12](https://github.com/nf-core/sarek/pull/12) - Simplify `check_max()` function
- [#13](https://github.com/nf-core/sarek/pull/13) - Merge `BamQCmapped` and `BamQCrecalibrated` processes into `BamQC` process
- [#13](https://github.com/nf-core/sarek/pull/13) - Split `CompressVCF` process into `CompressVCFsnpEff` and `CompressVCFvep` processes
- [#16](https://github.com/nf-core/sarek/pull/16) - Make scripts in `bin/` and `scripts/` executable
- [#18](https://github.com/nf-core/sarek/pull/18) - Use `--no-reports` for TravisCI testing
- [#18](https://github.com/nf-core/sarek/pull/18) - Add `--no-reports` for all tests but MULTIPLE in Jenkins
- [#18](https://github.com/nf-core/sarek/pull/18), [#29](https://github.com/nf-core/sarek/pull/29) - `--noReports` is now `--skipQC all`
- [#18](https://github.com/nf-core/sarek/pull/18), [#21](https://github.com/nf-core/sarek/pull/21) - Update logo
- [#21](https://github.com/nf-core/sarek/pull/21) - Moved `smallGRCh37` path to `genomes.config`
- [#23](https://github.com/nf-core/sarek/pull/23) - Rename `genomeFile`, `genomeIndex` and  `genomeDict` by `fasta`, `fastaFai` and `dict`
- [#23](https://github.com/nf-core/sarek/pull/23) - `--sample` is now deprecated, use `--input` instead
- [#23](https://github.com/nf-core/sarek/pull/23) - `--genomeFile` is now deprecated, use `--fasta` instead
- [#23](https://github.com/nf-core/sarek/pull/23) - `--genomeIndex` is now deprecated, use `--fastaFai` instead
- [#23](https://github.com/nf-core/sarek/pull/23) - `--genomeDict` is now deprecated, use `--dict` instead
- [#24](https://github.com/nf-core/sarek/pull/24) - `AWS iGenomes` config now contains germline resource for `GATK4 Mutect2`
- [#30](https://github.com/nf-core/sarek/pull/30) - Simplify code for `MapReads` process
- [#24](https://github.com/nf-core/sarek/pull/24) - `AWS iGenomes` config now contains germline resource for `GATK4 Mutect2`
- [#31](https://github.com/nf-core/sarek/pull/31) - Move extra CI to `GitHub Actions` nf-core extra CI
- [#32](https://github.com/nf-core/sarek/pull/32), [#33](https://github.com/nf-core/sarek/pull/33) - Install `ASCAT` with `conda` in the `environment.yml` file
- [#33](https://github.com/nf-core/sarek/pull/33) - Use `workflow.manifest.version` to specify workflow version in path to scripts for `ControlFREEC` and `VEP` processes
- [#35](https://github.com/nf-core/sarek/pull/35) - Building indexes is now done in `main.nf`
- [#35](https://github.com/nf-core/sarek/pull/35) - `build.nf` script now only download cache, so renamed to `downloadcache.nf`
- [#35](https://github.com/nf-core/sarek/pull/35) - Use `tabix` instead of `IGVtools` to build vcf indexes
- [#35](https://github.com/nf-core/sarek/pull/35) - Refactor references handling
- [#35](https://github.com/nf-core/sarek/pull/35) - Use Channel values instead of `referenceMap`
- [#37](https://github.com/nf-core/sarek/pull/37) - Bump version for Release
- [#38](https://github.com/nf-core/sarek/pull/38) - File names before merge is based on `${idSample}_${idRun}` instead of `${idRun}`

### Removed

- [#9](https://github.com/nf-core/sarek/pull/9) - Removed `relatedness2` graph from `vcftools stats`
- [#13](https://github.com/nf-core/sarek/pull/13) - Removed `BamQCmapped` and `BamQCrecalibrated` processes
- [#13](https://github.com/nf-core/sarek/pull/13) - Removed `CompressVCF`
- [#18](https://github.com/nf-core/sarek/pull/18) - Removed params `--noReports`
- [#24](https://github.com/nf-core/sarek/pull/18) - Removed `GATK3.X Mutect2`
- [#31](https://github.com/nf-core/sarek/pull/31) - Remove extra CI from `Travis CI` and `GitHub Actions` nf-core CI
- [#32](https://github.com/nf-core/sarek/pull/32), [#35](https://github.com/nf-core/sarek/pull/35) - Clean up `environment.yml` file
- [#35](https://github.com/nf-core/sarek/pull/35) - Remove building indexes from `build.nf` script
- [#35](https://github.com/nf-core/sarek/pull/35) - Remove helper script `build_reference.sh`
- [#35](https://github.com/nf-core/sarek/pull/35) - Remove `IGVtools`
- [#35](https://github.com/nf-core/sarek/pull/35) - Remove `GATK Mutect2` from `MULTIPLE` test
- [#35](https://github.com/nf-core/sarek/pull/35) - Remove `referenceMap` and `defineReferenceMap()` and use Channel values instead

### Fixed

- [#3](https://github.com/nf-core/sarek/pull/3) - Fix `Docker` ownership
- [#11](https://github.com/nf-core/sarek/pull/11) - Fix `MergeMpileup` PublishDir
- [#13](https://github.com/nf-core/sarek/pull/13) - Fix merge in annotation
- [#14](https://github.com/nf-core/sarek/pull/14) - Fix output name for vcf files
- [#16](https://github.com/nf-core/sarek/pull/16) - Fix path to `Rscript`
- [#18](https://github.com/nf-core/sarek/pull/18) - Improve cpu usage
- [#18](https://github.com/nf-core/sarek/pull/18) - Use same font for `nf-core` and `sarek` in ascii art
- [#20](https://github.com/nf-core/sarek/pull/20) - Use new logo in README
- [#20](https://github.com/nf-core/sarek/pull/20) - Fix path to references genomes
- [#22](https://github.com/nf-core/sarek/pull/22) - Fix `--singleCPUMem` issue
- [#30](https://github.com/nf-core/sarek/pull/30) - Fix choice between `inputPairReadsFastQC` and `inputBAMFastQC` channels
- [#31](https://github.com/nf-core/sarek/pull/31) - Fix badges according to nf-core lint
- [#31](https://github.com/nf-core/sarek/pull/31) - Fix `rcolorbrewer` version according to nf-core lint
- [#33](https://github.com/nf-core/sarek/pull/33) - Fix MD Linting
- [#38](https://github.com/nf-core/sarek/pull/38) - Avoid collision in `MultiQC`
- [#39](https://github.com/nf-core/sarek/pull/39) - Fix `ch_dbsnp` channel

### Deprecated

- [#23](https://github.com/nf-core/sarek/pull/23) - `--sample` is now deprecated, use `--input` instead
- [#23](https://github.com/nf-core/sarek/pull/23) - `--genomeFile` is now deprecated, use `--fasta` instead
- [#23](https://github.com/nf-core/sarek/pull/23) - `--genomeIndex` is now deprecated, use `--fastaFai` instead
- [#23](https://github.com/nf-core/sarek/pull/23) - `--genomeDict` is now deprecated, use `--dict` instead
- [#29](https://github.com/nf-core/sarek/pull/29) - `--noReports` is now deprecated, use `--skipQC all`

## [2.3.FIX1](https://github.com/SciLifeLab/Sarek/releases/tag/2.3.FIX1) - 2019-03-04

### Fixed

- [#742](https://github.com/SciLifeLab/Sarek/pull/742) - Fix output dirs (`HaplotypeCaller` that was not recognized by `annotate.nf` introduced by [#728](https://github.com/SciLifeLab/Sarek/pull/728))

## [2.3](https://github.com/SciLifeLab/Sarek/releases/tag/2.3) - Äpar - 2019-02-27

Äpar is one of the main massif in the Sarek National Park.

### Added

- [#628](https://github.com/SciLifeLab/Sarek/pull/628), [#722](https://github.com/SciLifeLab/Sarek/pull/722) - `ASCAT` now use `.gc` file
- [#712](https://github.com/SciLifeLab/Sarek/pull/712), [#718](https://github.com/SciLifeLab/Sarek/pull/718) - Added possibilities to run Sarek with `conda`
- [#719](https://github.com/SciLifeLab/Sarek/pull/719) - Annotation documentation
- [#719](https://github.com/SciLifeLab/Sarek/pull/719) - Helper script to download `snpeff` and `VEP` cache files
- [#719](https://github.com/SciLifeLab/Sarek/pull/719) - New `--annotation_cache`, `--snpEff_cache`, `--vep_cache` parameters
- [#719](https://github.com/SciLifeLab/Sarek/pull/719) - Possibility to use cache wen annotating with `snpEff` and `VEP`
- [#722](https://github.com/SciLifeLab/Sarek/pull/722) - Add path to ASCAT `.gc` file in `igenomes.config`
- [#728](https://github.com/SciLifeLab/Sarek/pull/728) - Update `Sarek-data` submodule with multiple patients TSV file
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Add `cadd_WG_SNVs`, `cadd_WG_SNVs_tbi`, `cadd_InDels`, `cadd_InDels_tbi` and `cadd_cache` params
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Add `tabix` indexed cache for `VEP`
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - New `DownloadCADD` process to download CADD files
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Specify values for `cadd_WG_SNVs`, `cadd_WG_SNVs_tbi`, `cadd_InDels`, `cadd_InDels_tbi` and `cadd_cache` params in `munin.conf` file
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Use `cadd_cache` param for optional use of CADD VEP plugin in `annotate.nf`
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - `VEP` cache has now fasta files for `--HGVS`
- [#735](https://github.com/SciLifeLab/Sarek/pull/735) - Added `--exome` for `Manta`, and for `StrelkaBP`
- [#735](https://github.com/SciLifeLab/Sarek/pull/735) - Added `Travis CI` test for targeted

### Changed

- [#710](https://github.com/SciLifeLab/Sarek/pull/710) - Improve release checklist and script
- [#711](https://github.com/SciLifeLab/Sarek/pull/711) - Improve configuration priorities
- [#716](https://github.com/SciLifeLab/Sarek/pull/716) - Update paths to containers and `AWS iGenomes`
- [#717](https://github.com/SciLifeLab/Sarek/pull/717) - `checkFileExtension` has changed to `hasExtension`, and now only verify if file has extension
- [#717](https://github.com/SciLifeLab/Sarek/pull/717) - `fastqFiles` renamed to `inputFiles`
- [#717](https://github.com/SciLifeLab/Sarek/pull/717) - `mapping` step can now map BAM files too
- [#717](https://github.com/SciLifeLab/Sarek/pull/717) - `MapReads` can now convert BAM to FASTQ and feed it to BWA on the fly
- [#717](https://github.com/SciLifeLab/Sarek/pull/717), [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Update documentation
- [#719](https://github.com/SciLifeLab/Sarek/pull/719) - `snpeff` and `vep` containers are now built with `conda`
- [#719](https://github.com/SciLifeLab/Sarek/pull/719) - `vepCacheVersion` is now defined in `conf/genomes.config` or `conf/igenomes.config`
- [#722](https://github.com/SciLifeLab/Sarek/pull/722) - Add path to ASCAT `.gc` file in `igenomes.config`
- [#722](https://github.com/SciLifeLab/Sarek/pull/722) - Update `Sarek-data` submodule
- [#723](https://github.com/SciLifeLab/Sarek/pull/723), [#725](https://github.com/SciLifeLab/Sarek/pull/725) - Update docs
- [#724](https://github.com/SciLifeLab/Sarek/pull/724) - Improved `AWS batch` configuration
- [#728](https://github.com/SciLifeLab/Sarek/pull/728) - Improved usage of `targetBED` params
- [#728](https://github.com/SciLifeLab/Sarek/pull/728) - `Strelka` Best Practices output is now prefixed with `StrelkaBP_`
- [#728](https://github.com/SciLifeLab/Sarek/pull/728) - VCFs and Annotated VCFs are now ordered by Patient, then tools
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Merge `buildContainers.nf` and `buildReferences.nf` in `build.nf`
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Reduce number of CPUs for `RunVEP` to `4` cf: [VEP docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#faster)
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Update `VEP` from `95.1` to `95.2`

### Removed

- [#715](https://github.com/SciLifeLab/Sarek/pull/715) - Remove `defReferencesFiles` function from `buildReferences.nf`
- [#719](https://github.com/SciLifeLab/Sarek/pull/719) - `snpEff` base container is no longer used
- [#721](https://github.com/SciLifeLab/Sarek/pull/721) - Remove `COSMIC` docs
- [#728](https://github.com/SciLifeLab/Sarek/pull/728) - Remove `defineDirectoryMap()`
- [#732](https://github.com/SciLifeLab/Sarek/pull/732) - Remove `--database` option for VEP cf: [VEP docs](https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#faster)

### Fixed

- [#720](https://github.com/SciLifeLab/Sarek/pull/720) - `bamQC` is now run on the recalibrated BAMs, and not after `MarkDuplicates`
- [#726](https://github.com/SciLifeLab/Sarek/pull/726) - Fix `Ascat` ref file input (one file can't be a set)
- [#727](https://github.com/SciLifeLab/Sarek/pull/727) - `bamQC` outputs are no longer overwritten (name of dir is now the file instead of sample)
- [#728](https://github.com/SciLifeLab/Sarek/pull/728) - Fix issue with annotation that was consuming `cache` channels
- [#728](https://github.com/SciLifeLab/Sarek/pull/728) - Fix multi sample TSV file [#691](https://github.com/SciLifeLab/Sarek/issues/691)
- [#733](https://github.com/SciLifeLab/Sarek/pull/733) - Fix the possibility to specify reference files on the command line

## [2.2.2](https://github.com/SciLifeLab/Sarek/releases/tag/2.2.2) - 2018-12-19

### Added

- [#671](https://github.com/SciLifeLab/Sarek/pull/671) - New `publishDirMode` param and docs
- [#673](https://github.com/SciLifeLab/Sarek/pull/673), [#675](https://github.com/SciLifeLab/Sarek/pull/675),  [#676](https://github.com/SciLifeLab/Sarek/pull/676) - Profiles for BinAC and CFC clusters in Tübingen
- [#679](https://github.com/SciLifeLab/Sarek/pull/679) - Add container for `CreateIntervalBeds`
- [#692](https://github.com/SciLifeLab/Sarek/pull/692), [#697](https://github.com/SciLifeLab/Sarek/pull/697) - Add `AWS iGenomes` possibilities (within `conf/igenomes.conf`)
- [#694](https://github.com/SciLifeLab/Sarek/pull/694) - Add monochrome and grey logos for light or dark background
- [#698](https://github.com/SciLifeLab/Sarek/pull/698) - Add btb profile for munin server
- [#702](https://github.com/SciLifeLab/Sarek/pull/702) - Add `font-ttf-dejavu-sans-mono` `2.37` and `fontconfig` `2.1dev` to container

### Changed

- [#663](https://github.com/SciLifeLab/Sarek/pull/663) - Update `do_release.sh` script
- [#671](https://github.com/SciLifeLab/Sarek/pull/671) - `publishDir` modes are now params
- [#677](https://github.com/SciLifeLab/Sarek/pull/677), [#698](https://github.com/SciLifeLab/Sarek/pull/698), [#703](https://github.com/SciLifeLab/Sarek/pull/703) - Update docs
- [#678](https://github.com/SciLifeLab/Sarek/pull/678) - Changing `VEP` to `v92` and adjusting CPUs for `VEP`
- [#679](https://github.com/SciLifeLab/Sarek/pull/679) - Update old `awsbatch` configuration
- [#682](https://github.com/SciLifeLab/Sarek/pull/682) - Specifications for memory and cpus for `awsbatch`
- [#693](https://github.com/SciLifeLab/Sarek/pull/693) - `Qualimap bamQC` is now ran after mapping and after recalibration for better QC
- [#700](https://github.com/SciLifeLab/Sarek/pull/700) - Update `GATK` to `4.0.9.0`
- [#702](https://github.com/SciLifeLab/Sarek/pull/702) - Update `FastQC` to `0.11.8`
- [#705](https://github.com/SciLifeLab/Sarek/pull/705) - Change `--TMP_DIR` by `--tmp-dir` for `GATK` `4.0.9.0` `BaseRecalibrator`
- [#706](https://github.com/SciLifeLab/Sarek/pull/706) - Update `Travis CI` testing

### Fixed

- [#665](https://github.com/SciLifeLab/Sarek/pull/665) - Input bam file now has always the same name (whether it is from a single fastq pair or multiple) in the `MarkDuplicates` process, so metrics too
- [#672](https://github.com/SciLifeLab/Sarek/pull/672) - Process `PullSingularityContainers` from `buildContainers.nf` now expect a file with the correct `.simg` extension for singularity images, and no longer the `.img` one.
- [#679](https://github.com/SciLifeLab/Sarek/pull/679) - Add `publishDirMode` for `germlineVC.nf`
- [#700](https://github.com/SciLifeLab/Sarek/pull/700) - Fix [#699](https://github.com/SciLifeLab/Sarek/issues/699) missing DP in the FORMAT column VCFs for Mutect2
- [#702](https://github.com/SciLifeLab/Sarek/pull/702) - Fix [#701](https://github.com/SciLifeLab/Sarek/issues/701)
- [#705](https://github.com/SciLifeLab/Sarek/pull/705) - Fix [#704](https://github.com/SciLifeLab/Sarek/issues/704)

## [2.2.1](https://github.com/SciLifeLab/Sarek/releases/tag/2.2.1) - 2018-10-04

### Changed

- [#646](https://github.com/SciLifeLab/Sarek/pull/646) - Update [`pathfindr`](https://github.com/NBISweden/pathfindr) submodule
- [#659](https://github.com/SciLifeLab/Sarek/pull/659) - Update `Nextflow` to `0.32.0`
- [#660](https://github.com/SciLifeLab/Sarek/pull/660) - Update docs

### Fixed

- [#657](https://github.com/SciLifeLab/Sarek/pull/657) - Fix `RunMultiQC.nf` bug
- [#659](https://github.com/SciLifeLab/Sarek/pull/659) - Fix bugs due to updating `Nextflow`

## [2.2.0](https://github.com/SciLifeLab/Sarek/releases/tag/2.2.0) - Skårki - 2018-09-21

Skårki is one of the main massif in the Sarek National Park.

### Added

- [#613](https://github.com/SciLifeLab/Sarek/pull/613) - Add Issue Templates (bug report and feature request)
- [#614](https://github.com/SciLifeLab/Sarek/pull/614) - Add PR Template
- [#615](https://github.com/SciLifeLab/Sarek/pull/615) - Add presentation
- [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Update documentation
- [#620](https://github.com/SciLifeLab/Sarek/pull/620) - Add `tmp/` to `.gitignore`
- [#625](https://github.com/SciLifeLab/Sarek/pull/625) - Add [`pathfindr`](https://github.com/NBISweden/pathfindr) as a submodule
- [#635](https://github.com/SciLifeLab/Sarek/pull/635) - To process targeted sequencing with a target BED
- [#639](https://github.com/SciLifeLab/Sarek/pull/639) - Add a complete example analysis to docs
- [#640](https://github.com/SciLifeLab/Sarek/pull/640), [#642](https://github.com/SciLifeLab/Sarek/pull/642) - Add helper script for changing version number

### Changed

- [#608](https://github.com/SciLifeLab/Sarek/pull/608) - Update `Nextflow` required version
- [#615](https://github.com/SciLifeLab/Sarek/pull/615) - Use `splitCsv` instead of `readlines`
- [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Update `CHANGELOG`
- [#621](https://github.com/SciLifeLab/Sarek/pull/621), [#638](https://github.com/SciLifeLab/Sarek/pull/638) - Improve install script
- [#621](https://github.com/SciLifeLab/Sarek/pull/621), [#638](https://github.com/SciLifeLab/Sarek/pull/638) - Simplify tests
- [#627](https://github.com/SciLifeLab/Sarek/pull/627), [#629](https://github.com/SciLifeLab/Sarek/pull/629), [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Refactor docs
- [#629](https://github.com/SciLifeLab/Sarek/pull/629) - Refactor config
- [#632](https://github.com/SciLifeLab/Sarek/pull/632) - Use 2 threads and 2 cpus `FastQC` processes
- [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Update tool version gathering
- [#638](https://github.com/SciLifeLab/Sarek/pull/638) - Use correct `.simg` extension for Singularity images
- [#639](https://github.com/SciLifeLab/Sarek/pull/639) - Smaller refactoring of the docs
- [#640](https://github.com/SciLifeLab/Sarek/pull/640) - Update RELEASE_CHECKLIST
- [#642](https://github.com/SciLifeLab/Sarek/pull/642) - `MultiQC` 1.5 -> 1.6
- [#642](https://github.com/SciLifeLab/Sarek/pull/642) - `Qualimap` 2.2.2a -> 2.2.2b
- [#642](https://github.com/SciLifeLab/Sarek/pull/642) - Update `conda` channel order priorities
- [#642](https://github.com/SciLifeLab/Sarek/pull/642) - `VCFanno` 0.2.8 -> 0.3.0
- [#642](https://github.com/SciLifeLab/Sarek/pull/642) - `VCFtools` 0.1.15 -> 0.1.16

### Removed

- [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Remove old Issue Template
- [#629](https://github.com/SciLifeLab/Sarek/pull/629) - Remove old Dockerfiles
- [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Remove old comments

### Fixed

- [#621](https://github.com/SciLifeLab/Sarek/pull/621) - Fix `VEP` tests
- [#637](https://github.com/SciLifeLab/Sarek/pull/637) - Fix links in MD files

## [2.1.0](https://github.com/SciLifeLab/Sarek/releases/tag/2.1.0) - Ruotes - 2018-08-14

Ruotes is one of the main massif in the Sarek National Park.

### Added

- [#555](https://github.com/SciLifeLab/Sarek/pull/555) - `snpEff` output into `VEP`
- [#556](https://github.com/SciLifeLab/Sarek/pull/556) - `Strelka` Best Practices
- [#563](https://github.com/SciLifeLab/Sarek/pull/563) - Use `SnpEFF` reports in `MultiQC`
- [#568](https://github.com/SciLifeLab/Sarek/pull/568) - `VCFTools` process `RunVcftools` for QC
- [#574](https://github.com/SciLifeLab/Sarek/pull/574), [#580](https://github.com/SciLifeLab/Sarek/pull/580) - Abstracts for `NPMI`, `JOBIM` and  `EACR25`
- [#577](https://github.com/SciLifeLab/Sarek/pull/577) - New repository for testing: [Sarek-data](https://github.com/SciLifeLab/Sarek-data)
- [#595](https://github.com/SciLifeLab/Sarek/pull/595) - New library `QC` for functions `bamQC`, `bcftools`, `samtoolsStats`, `vcftools`, `getVersionBCFtools`, `getVersionGATK`, `getVersionManta`, `getVersionSnpEFF`, `getVersionStrelka`, `getVersionVCFtools`, `getVersionVEP`
- [#595](https://github.com/SciLifeLab/Sarek/pull/595) - New Processes `GetVersionBCFtools`, `GetVersionGATK`, `GetVersionManta`, `GetVersionSnpEFF`, `GetVersionStrelka`, `GetVersionVCFtools`, `GetVersionVEP`
- [#595](https://github.com/SciLifeLab/Sarek/pull/595) - New `Python` script `bin/scrape_tool_versions.py` inspired by @ewels and @apeltzer
- [#595](https://github.com/SciLifeLab/Sarek/pull/595) - New QC Process `RunVcftools`
- [#596](https://github.com/SciLifeLab/Sarek/pull/596) - New profile for `BinAC` cluster
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - New function `sarek_ascii()` in `SarekUtils`
- [#599](https://github.com/SciLifeLab/Sarek/pull/599), [#602](https://github.com/SciLifeLab/Sarek/pull/602) - New Process `CompressVCF`
- [#601](https://github.com/SciLifeLab/Sarek/pull/601), [#603](https://github.com/SciLifeLab/Sarek/pull/603) - Container for `GATK4`
- [#606](https://github.com/SciLifeLab/Sarek/pull/606) - Add test data as a submodule from [`Sarek-data`](https://github.com/SciLifeLab/Sarek-data)
- [#608](https://github.com/SciLifeLab/Sarek/pull/608) - Add documentation on how to install Nextflow on `bianca`

### Changed

- [#557](https://github.com/SciLifeLab/Sarek/pull/557), [#583](https://github.com/SciLifeLab/Sarek/pull/583), [#585](https://github.com/SciLifeLab/Sarek/pull/585), [#588](https://github.com/SciLifeLab/Sarek/pull/588) - Update help
- [#560](https://github.com/SciLifeLab/Sarek/pull/560) - `GitHub` langage for the repository is now `Nextflow`
- [#561](https://github.com/SciLifeLab/Sarek/pull/561) - `do_all.sh` build only containers for one genome reference (default `GRCh38`) only
- [#571](https://github.com/SciLifeLab/Sarek/pull/571) - Only one container for all QC tools
- [#582](https://github.com/SciLifeLab/Sarek/pull/582), [#587](https://github.com/SciLifeLab/Sarek/pull/587) - Update figures
- [#595](https://github.com/SciLifeLab/Sarek/pull/595) - Function `defineDirectoryMap()` is now part of `SarekUtils`
- [#595](https://github.com/SciLifeLab/Sarek/pull/595) - Process `GenerateMultiQCconfig` replace by function `createMultiQCconfig()`
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - `extractBams()` now takes an extra parameter.
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Move `checkFileExtension()`, `checkParameterExistence()`, `checkParameterList()`, `checkReferenceMap()`, `checkRefExistence()`, `extractBams()`, `extractGenders()`, `returnFile()`, `returnStatus()` and `returnTSV()` functions to `SarekUtils`
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Reduce data footprint for Process `CreateRecalibrationTable`
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Replace deprecated operator `phase` by `join`.
- [#599](https://github.com/SciLifeLab/Sarek/pull/599) - Merge is tested with `ANNOTATEALL`
- [#604](https://github.com/SciLifeLab/Sarek/pull/604) - Synching `GRCh38` `wgs_calling_regions` bedfiles
- [#607](https://github.com/SciLifeLab/Sarek/pull/607) - One container approach
- [#607](https://github.com/SciLifeLab/Sarek/pull/607) - Update to `GATK4`
- [#608](https://github.com/SciLifeLab/Sarek/pull/608) - Update `Nextflow` required version
- [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Update `CHANGELOG`
- [#617](https://github.com/SciLifeLab/Sarek/pull/617) - Replace deprecated `Nextflow ``$name` syntax with `withName`

### Fixed

- [#560](https://github.com/SciLifeLab/Sarek/pull/560) - Display message for `repository` and `containerPath`
- [#566](https://github.com/SciLifeLab/Sarek/pull/566) - `slurmDownload` profile
- [#579](https://github.com/SciLifeLab/Sarek/pull/579), [#584](https://github.com/SciLifeLab/Sarek/pull/584) - `Manta` output reorganized after modification for `Strelka Best Practices` process
- [#585](https://github.com/SciLifeLab/Sarek/pull/583) - Trace file is plain txt
- [#590](https://github.com/SciLifeLab/Sarek/pull/590), [#593](https://github.com/SciLifeLab/Sarek/pull/593) - Fix `Singularity` installation in `Travis CI` testing
- [#598](https://github.com/SciLifeLab/Sarek/pull/598), [#601](https://github.com/SciLifeLab/Sarek/pull/601) - Fixes for  `Python` script `selectROI.py` to work with `CLC` viewer

### Removed

- [#607](https://github.com/SciLifeLab/Sarek/pull/607) - Remove `Mutect1`

## [2.0.0](https://github.com/SciLifeLab/Sarek/releases/tag/2.0.0) - 2018-03-23

First release under the `Sarek` name, from the National Park in Northern Sweden

### Added

- Basic wrapper script
- Abstract, posters and figures
- ROI selector and `FreeBayes` sanitizer scripts
- New logo and icon for the project
- Check for existing tumor/normal channel
- `SarekUtils` with `checkParams()`, `checkParameterList()`, `checkParameterExistence()` and `isAllowedParams()` functions
- Some `runOptions` for `docker` (prevent some user right problem)
- This `CHANGELOG`

### Changed

- `CAW` is now `Sarek`
- Dissect Workflow in 5 new scripts: `annotate.nf`, `main.nf`, `germlineVC.nf`, `runMultiQC.nf` and `somaticVC.nf`
- `report.html`, `timeline.html` and `trace.html` are generated in `Reports/`
- `--version` is now used to define the workflow version
- Most params are now defined in the `base.config` file instead of in the scripts
- Update `RELEASE_CHECKLIST.md`
- `checkParams()`, `checkParameterList()`, `checkParameterExistence()` and `isAllowedParams()` in script functions are now called within `SarekUtils`
- `nf_required_version` is now `params.nfRequiredVersion`
- In `buildReferences.nf` script, channels now begin by `ch_`, and files by `f_`
- Use `PublishDir mode: 'link'` instead of `copy`
- `directoryMap` now contains `params.outDir`
- [#539](https://github.com/SciLifeLab/Sarek/issues/539) - Use Nextflow support of scratch
- Reordered `Travis CI` tests
- Update documentation
- `MultiQC` version in container from v`1.4` to v`1.5`
- `vepgrch37` container base image from `release_90.6` to `release_92`
- `vepgrch38` container base image from `release_90.6` to `release_92`
- `VEP` version in containers from v`90` to v`91`
- `nucleotidesPerSecond` is now `params.nucleotidesPerSecond`
- Default `params.tag` is now `latest` instead of current version, so `--tag` needs to be specified with the right version to be sure of using the `containers` corresponding

### Deprecated

- `standard` profile
- `uppmax-localhost.config` file

### Removed

- `scripts/skeleton_batch.sh`
- Old data and tsv files
- `UPPMAX` directories from containers
- `--step` in `annotate.nf`, `germlineVC.nf` and `somatic.nf`
- Some `runOptions` for `Singularity` (binding not needed anymore on `UPPMAX`)
- `download` profile

### Fixed

- [#530](https://github.com/SciLifeLab/Sarek/issues/530) - Use `$PWD` for default `outDir`
- [#533](https://github.com/SciLifeLab/Sarek/issues/533) - Replace `VEP` `--pick` option by `--per_gene`

## [1.2.5](https://github.com/SciLifeLab/Sarek/releases/tag/1.2.5) - 2018-01-18

### Added

- `Zenodo` for DOI
- Delivery README
- Document use of the `--sampleDir` option
- Contributing Guidelines
- Issue Templates
- Release Checklist
- `--outDir`
- `awsbatch` profile
- `aws-batch.config` config file
- `--noBAMQC` params (failing sometimes on `Bianca`)

### Changed

- Update `Nextflow` to `0.26.0` (new fancy report + `AWS Batch`)
- Extra time on `Travis CI` testing
- Replace `bundleDir` by `params.genome_base`
- Update `MultiQC` to `1.3` (`MEGAQC` FTW)
- Move and rename some test files

### Fixed

- Version of `COSMIC` `GRCh37` `v83`
- Write an error message when `--sampleDir` does not find any FASTQ files
- `base.config` for `ConcatVCF` process
- File specification for `recalibrationReport` in `RecalibrateBam` process (got error on `AWS Batch`)

## [1.2.4](https://github.com/SciLifeLab/Sarek/releases/tag/1.2.4) - 2017-10-27

### Fixed

- [#488](https://github.com/SciLifeLab/Sarek/issues/488) - Better CPU requirements for `ConcatVCF`
- [#489](https://github.com/SciLifeLab/Sarek/issues/489) - Exception handling for `ASCAT`
- [#490](https://github.com/SciLifeLab/Sarek/issues/490) - CPU requirements for `runSingleStrelka` and `runSingleManta`

## [1.2.3](https://github.com/SciLifeLab/Sarek/releases/tag/1.2.3) - 2017-10-18

### Fixed

- [#357](https://github.com/SciLifeLab/Sarek/issues/357) - `ASCAT` works for `GRCh38`
- [#471](https://github.com/SciLifeLab/Sarek/issues/471) - Running `Singularity` on `/scratch`
- [#475](https://github.com/SciLifeLab/Sarek/issues/475) - 16 cpus for local executor
- [#480](https://github.com/SciLifeLab/Sarek/issues/480) - No `tsv` file needed for step `annotate`

## [1.2.2](https://github.com/SciLifeLab/Sarek/releases/tag/1.2.2) - 2017-10-06

### Fixed

- [#479](https://github.com/SciLifeLab/Sarek/issues/479) - Typo in `uppmax-localhost.config`

## [1.2.1](https://github.com/SciLifeLab/Sarek/releases/tag/1.2.1) - 2017-10-06

### Changed

- `runascat` and `runconvertallelecounts` containers are now replaced by `r-base`
- `willmclaren/ensembl-vep:release_90.5` is now base for `vepgrch37` and `vepgrch38`

### Removed

- `vep` container
- `strelka_config.ini` file

### Fixed

- [#471](https://github.com/SciLifeLab/Sarek/issues/471) - Running `Singularity` on /scratch
- [#472](https://github.com/SciLifeLab/Sarek/issues/472) - Update function to check `Nextflow` version
- [#473](https://github.com/SciLifeLab/Sarek/issues/473) - Remove `returnMin()` function

## [1.2.0](https://github.com/SciLifeLab/Sarek/releases/tag/1.2.0) - 2017-10-02

### Changed

- Fix version for Manuscript

## [1.1](https://github.com/SciLifeLab/Sarek/releases/tag/1.1) - 2017-09-15

### Added

- `Singularity` possibilities

### Changed

- Reports made by default
- Intervals file can be a bed file
- Normal sample preprocessing + `HaplotypeCaller` is possible
- Better `Travis CI` tests

### Fixed

- Memory requirements

## [1.0](https://github.com/SciLifeLab/Sarek/releases/tag/1.0) - 2017-02-16

### Added

- `Docker` possibilities

## [0.9](https://github.com/SciLifeLab/Sarek/releases/tag/0.9) - 2016-11-16

## [0.8](https://github.com/SciLifeLab/Sarek/releases/tag/0.8) - 2016-11-16

## [0.1] - 2016-04-05
