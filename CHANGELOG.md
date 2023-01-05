# nf-core/sarek: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.2](https://github.com/nf-core/sarek/releases/tag/3.1.2) - Lesser Lule River

Lesser Lule River is English for Lilla Luleälven

### Added

### Changed

### Fixed

- [#905](https://github.com/nf-core/sarek/pull/905) - Remove usages of deprecated `Channel.from` method

### Deprecated

### Removed

### Dependencies

## [3.1.1](https://github.com/nf-core/sarek/releases/tag/3.1.1) - Lilla Luleälven

Lilla Luleälven river's main affluent is Rapaätno.

### Added

- [#856](https://github.com/nf-core/sarek/pull/856) - Add annotation for `R64-1-1` and `UMD3.1`

### Changed

- [#855](https://github.com/nf-core/sarek/pull/855) - Speed up duplicate marking by using `samtools` for CRAM conversion
- [#858](https://github.com/nf-core/sarek/pull/858) - Prepare release `3.1.1`

### Fixed

- [#851](https://github.com/nf-core/sarek/pull/851) - Fix `schema` definition `None` for `cf_chrom_len`

### Deprecated

### Removed

### Dependencies

## [3.1](https://github.com/nf-core/sarek/releases/tag/3.1) - Rapaätno

Rapaätno is the river you can see from the Skierfe mountain.

### Added

- [#735](https://github.com/nf-core/sarek/pull/735) - GATK Markduplicates now natively supports CRAM output
- [#774](https://github.com/nf-core/sarek/pull/774) - Add logo for Danish National Genome Center
- [#783](https://github.com/nf-core/sarek/pull/783) - Add paths for chr length used by controlfreec to GRCh38 config
- [#820](https://github.com/nf-core/sarek/pull/820) - Improve documentation on scatter/gather effects
- [#833](https://github.com/nf-core/sarek/pull/833) - Add name to CI tests to avoid confusion between runs

### Changed

- [#735](https://github.com/nf-core/sarek/pull/735) - `--save_mapped` now saves mapping output in CRAM format
- [#762](https://github.com/nf-core/sarek/pull/762) - Back to dev
- [#762](https://github.com/nf-core/sarek/pull/762) - Update deepvariant module
- [#773](https://github.com/nf-core/sarek/pull/773) - Sync `TEMPLATE` with `tools` `2.6`
- [#782](https://github.com/nf-core/sarek/pull/782) - Reduce scatter/gather for full size tests on AWS
- [#785](https://github.com/nf-core/sarek/pull/785) - Update description of `bcftools stats`
- [#784](https://github.com/nf-core/sarek/pull/784) - Update all subworkflows names thanks to @scorreard
- [#806](https://github.com/nf-core/sarek/pull/806) - Refactor all tests
- [#806](https://github.com/nf-core/sarek/pull/806) - Split up `modules.config` file
- [#810](https://github.com/nf-core/sarek/pull/810) - Update CHANGELOG
- [#821](https://github.com/nf-core/sarek/pull/821) - Change `replace` to `putIfAbsent` for automatic search of `input` if none is provided to avoid overwriting values
- [#822](https://github.com/nf-core/sarek/pull/822) - Update modules with `nf-core modules update -a`: Update GATK version to 4.3.0
- [#827](https://github.com/nf-core/sarek/pull/827) - Add `--genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader` to `GenomicsDB` parameters to speed up the analysis
- [#842](https://github.com/nf-core/sarek/pull/842) - Increase default memory for samtools stats
- [#844](https://github.com/nf-core/sarek/pull/844) - All small scale tests are run on PR to `master`

### Fixed

- [#762](https://github.com/nf-core/sarek/pull/762) - Polish CHANGELOG + figures
- [#766](https://github.com/nf-core/sarek/pull/766) - Align box description in subway map
- [#768](https://github.com/nf-core/sarek/pull/768) - Use double quotes to fix import of singularity images for deepvariant module
- [#770](https://github.com/nf-core/sarek/pull/770) - Use double quotes to fix import of singularity images for gatk4/cnnscorevariants module
- [#771](https://github.com/nf-core/sarek/pull/771) - update to new modules syntax
- [#777](https://github.com/nf-core/sarek/pull/777) - Fix mixed up aws full size tests output paths
- [#790](https://github.com/nf-core/sarek/pull/790) - Fix issue [#789](https://github.com/nf-core/sarek/issues/789) somatic mutect2 test
- [#793](https://github.com/nf-core/sarek/pull/793) - Remove DeepVariant GVCF from annotation
- [#794](https://github.com/nf-core/sarek/pull/794) - Fix publishing for unzipped reference files
- [#807](https://github.com/nf-core/sarek/pull/807) - Fix read group when uBAMs are provided (see issue [#732](https://github.com/nf-core/sarek/issues/732))
- [#813](https://github.com/nf-core/sarek/pull/813) - Fix input validation when launching from website (see issue [#694](https://github.com/nf-core/sarek/issues/694))
- [#814](https://github.com/nf-core/sarek/pull/814) - Fix readgroups when using DragMap together with FreeBayes or Mutect2 (see issue [#780](https://github.com/nf-core/sarek/issues/780))
- [#817](https://github.com/nf-core/sarek/pull/817) - Fix CNVKit run on tumor-only sample to be run on all samples
- [#828](https://github.com/nf-core/sarek/pull/817) - Fix issue [#763](https://github.com/nf-core/sarek/issues/763) to run variantcalling when starting form step recalibration
- [#837](https://github.com/nf-core/sarek/pull/837) - Fix Freebayes config selector after subworkflow renaming
- [#839](https://github.com/nf-core/sarek/pull/839) - Remove `copyTo` method that fails on S3 when the source and destination buckets are in different regions
- [#841](https://github.com/nf-core/sarek/pull/841) - Fix path priority for `cf_chrom_len`

### Deprecated

### Removed

### Dependencies

| Dependency    | Old version | New version |
| ------------- | ----------- | ----------- |
| `bcftools`    | 1.15.1      | 1.16        |
| `deepvariant` | 1.3.0       | 1.4.0       |
| `freebayes`   | 1.3.5       | 1.3.6       |
| `gatk4`       | 4.2.6.1     | 4.3.0.0     |
| `samtools`    | 1.15.1      | 1.16.1      |
| `tiddit`      | 3.1.0       | 3.3.2       |

## [3.0.2](https://github.com/nf-core/sarek/releases/tag/3.0.2) - Lájtávrre

Lájtávrre is a lake you can see from the Skierfe mountain, formed by the Rapaätno river.

### Added

- [#691](https://github.com/nf-core/sarek/pull/691) - Enable `PROFILE=conda`, `PROFILE=docker` and `PROFILE=singularity` for pytest
- [#716](https://github.com/nf-core/sarek/pull/716) - Add documentation for Azure recommended config vm_size
- [#752](https://github.com/nf-core/sarek/pull/752) - Add tracking of all dependencies starting 3.0

### Changed

- [#679](https://github.com/nf-core/sarek/pull/679) - Back to `dev`
- [#685](https://github.com/nf-core/sarek/pull/685) - Updating the nf-core modules used by Sarek
- [#691](https://github.com/nf-core/sarek/pull/691) - To run the same pytest as before locally, use `PROFILE=docker`
- [#692](https://github.com/nf-core/sarek/pull/692) - Use `params.tools=strelka` in profile `test`
- [#696](https://github.com/nf-core/sarek/pull/696) - Adding check of md5-sums in CI-tests
- [#719](https://github.com/nf-core/sarek/pull/719) - Added boxes to subway map
- [#720](https://github.com/nf-core/sarek/pull/720) - Sync `TEMPLATE` with `tools` `2.5`
- [#723](https://github.com/nf-core/sarek/pull/723) - Sync `TEMPLATE` with `tools` `2.5.1`
- [#726](https://github.com/nf-core/sarek/pull/726) - Adapt resource requests
- [#730](https://github.com/nf-core/sarek/pull/730) - Reduce number of tests
- [#731](https://github.com/nf-core/sarek/pull/731) - Run the somatic test as default on `-profile test_full`, the germline can be tested with `-profile test_full_germline`
- [#733](https://github.com/nf-core/sarek/pull/733) - Add description for params.cf_chrom_len
- [#734](https://github.com/nf-core/sarek/pull/734) - nf-core modules update -a
- [#736](https://github.com/nf-core/sarek/pull/736) - More extensive CI for default test
- [#742](https://github.com/nf-core/sarek/pull/742) - Requiring the Haplotypecaller to be specified as one of the tools for joint germline genotyping
- [#752](https://github.com/nf-core/sarek/pull/752) - Code polishing

### Fixed

- [#679](https://github.com/nf-core/sarek/pull/679) - Fixed typos in subway maps
- [#681](https://github.com/nf-core/sarek/pull/681) - Fixed intermediate files published cf [#680](https://github.com/nf-core/sarek/issues/680)
- [#688](https://github.com/nf-core/sarek/pull/688) - Fixed VEP plugins issue cf [#687](https://github.com/nf-core/sarek/issues/687)
- [#689](https://github.com/nf-core/sarek/pull/689) - Fixed when clause for non `BWA mem` building mapping indexes
- [#704](https://github.com/nf-core/sarek/pull/704) - Fixed `cf_ploidy` to string instead of number
- [#705](https://github.com/nf-core/sarek/pull/705) - Fix publishing for processes in `alignment_to_fastq` subworkflow; prevent tabix computation for `known_snps` when present; publish `umi` processed files into `preprocessing/umi` subdirectory
- [#706](https://github.com/nf-core/sarek/pull/706) - Fixed `vep_version` not found error when running `--vep_loftee`
- [#724](https://github.com/nf-core/sarek/pull/724) - Fixed prettier issue
- [#727](https://github.com/nf-core/sarek/pull/727) - Allow `.list` interval files; remove `seconds` from GRCh38 file to allow `--nucleotides_per_second` to be used
- [#728](https://github.com/nf-core/sarek/pull/728) - Circumvent issue with controlfreec and length file containing regions not in intervals file
- [#729](https://github.com/nf-core/sarek/pull/729) - Trailing commas in `--tools`, `--skip_tools` and `--use_gatk_spark` now raise failure cf [#722](https://github.com/nf-core/sarek/issues/722)
- [#741](https://github.com/nf-core/sarek/pull/741) - Fix prefix for `bcftools sort` for joint germline variant calling
- [#743](https://github.com/nf-core/sarek/pull/743) - Remove profile definitions in profile to avoid issues with Tower
- [#758](https://github.com/nf-core/sarek/pull/758) - Fix Zenodo batch
- [#760](https://github.com/nf-core/sarek/pull/760) - Fix CHANGELOG dependencies
- [#761](https://github.com/nf-core/sarek/pull/761) - Fix font in subway map and workflow image

### Deprecated

### Removed

- [#742](https://github.com/nf-core/sarek/pull/742) - Removed some lines from the usage-doc as Sarek no longer support input supplied as a list of multiple csv-files
- [#757](https://github.com/nf-core/sarek/pull/757) - Remove `errorStrategy` in `conf/modules.config`

## [3.0.1](https://github.com/nf-core/sarek/releases/tag/3.0.1) - Saiva

Saiva is a lake in the Sarek national park, just below the Skierfe mountain.

### Fixed

- [#708](https://github.com/nf-core/sarek/pull/708) - Fixes mpileup bug. Update nf-core module `samtools/mpileup` to subset CRAM file by intervals

## [3.0](https://github.com/nf-core/sarek/releases/tag/3.0) - Skierfe

Skierfe is a mountain in the Sarek national park, and the inspiration for the logo.

### Added

- [#388](https://github.com/nf-core/sarek/pull/388) - Add cram support + read splitting with `SeqKit` for speedup
- [#394](https://github.com/nf-core/sarek/pull/394) - Add `DeepVariant`
- [#411](https://github.com/nf-core/sarek/pull/411) - cram in csv samplesheet
- [#448](https://github.com/nf-core/sarek/pull/448) - Allow to skip base quality recalibration with `--skip_bqsr`
- [#449](https://github.com/nf-core/sarek/pull/449) - [@FriederikeHanssen](https://github.com/FriederikeHanssen) is now a `CODEOWNERS`
- [#460](https://github.com/nf-core/sarek/pull/460) - Add posters
- [#463](https://github.com/nf-core/sarek/pull/463) - Add dark/light logo versions
- [#464](https://github.com/nf-core/sarek/pull/464), [#514](https://github.com/nf-core/sarek/pull/514) - Add `DRAGMAP` as a possible aligner
- [#479](https://github.com/nf-core/sarek/pull/479) - Add more subworkflows
- [#485](https://github.com/nf-core/sarek/pull/485) - `--skip_qc`, `--skip_markduplicates` and `--skip_bqsr` is now `--skip_tools`
- [#507](https://github.com/nf-core/sarek/pull/507), [#537](https://github.com/nf-core/sarek/pull/537) - Subway map for building indexes
- [#512](https://github.com/nf-core/sarek/pull/512), [#531](https://github.com/nf-core/sarek/pull/531), [#537](https://github.com/nf-core/sarek/pull/537) - Subway map for pipeline
- [#522](https://github.com/nf-core/sarek/pull/522) - Add QC for vcf files & MultiQC
- [#533](https://github.com/nf-core/sarek/pull/533) - Add param `--only_paired_variant_calling` to allow skipping of germline variantcalling for paired samples
- [#536](https://github.com/nf-core/sarek/pull/536) - Add `--step markduplicates` to start from duplicate marking, `--step prepare_recalibration` now ONLY starts at process `BaseRecalibrator` & adding `bam` and `cram` input support for `--step` `markduplicates`, `prepare_recalibration`, `recalibrate`, and `variant_calling`
- [#538](https://github.com/nf-core/sarek/pull/538) - Add param `--seq_platform`, default: `ILLUMINA`
- [#545](https://github.com/nf-core/sarek/pull/545) - Add modules and subworkflows for `cnvkit` tumor_only mode
- [#540](https://github.com/nf-core/sarek/pull/540) - Add modules and subworkflows for `cnvkit` somatic mode
- [#557](https://github.com/nf-core/sarek/pull/557) - Add `Haplotypecaller` single sample mode together with `CNNScoreVariants` and `FilterVariantTranches`
- [#576](https://github.com/nf-core/sarek/pull/576) - Add modules and subworkflows for `cnvkit` germline mode
- [#582](https://github.com/nf-core/sarek/pull/582) - Added option `--vep_out_format` for setting the format of the output-file from VEP to `json`, `tab` or `vcf` (default)
- [#594](https://github.com/nf-core/sarek/pull/594) - Add parameter `--save_output_as_bam` to allow output of result files in BAM format
- [#595](https://github.com/nf-core/sarek/pull/595) - Added Haplotypecaller joint germline calling
- [#597](https://github.com/nf-core/sarek/pull/597) - Added tiddit for tumor variant calling
- [#600](https://github.com/nf-core/sarek/pull/600) - Added description for UMI related params in schema
- [#604](https://github.com/nf-core/sarek/pull/604), [#617](https://github.com/nf-core/sarek/pull/617) - Added full size tests WGS 30x NA12878
- [#613](https://github.com/nf-core/sarek/pull/613) - Added params `--dbnsfp_fields` to allow configuration of fields for the `dbnsfp` `VEP` plugin
- [#613](https://github.com/nf-core/sarek/pull/613) - Added params `--dbnsfp_consequence` to allow configuration of consequence for the `dbnsfp` `VEP` plugin
- [#613](https://github.com/nf-core/sarek/pull/613) - Added params `--vep_version` to allow more configuration on the vep container definition
- [#620](https://github.com/nf-core/sarek/pull/620) - Added checks for sex information when running a CNV tools
- [#623](https://github.com/nf-core/sarek/pull/623) - Additional checks of data in the input sample sheet
- [#629](https://github.com/nf-core/sarek/pull/629) - Added checks to catch inconsistency between supplied samples and requested tools
- [#632](https://github.com/nf-core/sarek/pull/632) - Added params `--snpeff_version` to allow more configuration on the snpeff container definition
- [#632](https://github.com/nf-core/sarek/pull/632) - Added params `--vep_include_fasta` to use the fasta file for annotation
- [#639](https://github.com/nf-core/sarek/pull/639) - Adding genes-txt-file and summary-html-file to the published output from snpEff
- [#647](https://github.com/nf-core/sarek/pull/647) - Update resource requests for preprocessing based on what worked for 5 ICGC matched WGS samples
- [#652](https://github.com/nf-core/sarek/pull/652) - Added full size somatic test profile

### Changed

- [#580](https://github.com/nf-core/sarek/pull/580) - changed the test_full config to real public WXS data. 1 sample WXS germline, 1 Tumor/Normal pair. https://doi.org/10.1038/sdata.2016.25 and https://doi.org/10.1038/s41587-021-00994-5
- [#383](https://github.com/nf-core/sarek/pull/383), [#528](https://github.com/nf-core/sarek/pull/528) - Update `CHANGELOG`
- [#390](https://github.com/nf-core/sarek/pull/390) - Update `nextflow_schema.json`
- [#408](https://github.com/nf-core/sarek/pull/408) - Sync `TEMPLATE` with `tools` `2.0.1`
- [#416](https://github.com/nf-core/sarek/pull/416) - Sync `TEMPLATE` with `tools` `2.1`
- [#417](https://github.com/nf-core/sarek/pull/417) - Merge `dsl2` and `dev` branches
- [#419](https://github.com/nf-core/sarek/pull/419) - Improve preprocessing
- [#420](https://github.com/nf-core/sarek/pull/420), [#455](https://github.com/nf-core/sarek/pull/455), [#459](https://github.com/nf-core/sarek/pull/459), [#633](https://github.com/nf-core/sarek/pull/633) - `nf-core modules update --all`
- [#427](https://github.com/nf-core/sarek/pull/427) - Update `DeepVariant`
- [#462](https://github.com/nf-core/sarek/pull/462) - Update modules and `modules.config`
- [#465](https://github.com/nf-core/sarek/pull/465) - Improve `test_data.config`
- [#466](https://github.com/nf-core/sarek/pull/466), [#478](https://github.com/nf-core/sarek/pull/478), [#492](https://github.com/nf-core/sarek/pull/492), [#521](https://github.com/nf-core/sarek/pull/521) - Move some local modules to `nf-core/modules`
- [#466](https://github.com/nf-core/sarek/pull/466), [#485](https://github.com/nf-core/sarek/pull/485), [#492](https://github.com/nf-core/sarek/pull/492), [#494](https://github.com/nf-core/sarek/pull/494), [#515](https://github.com/nf-core/sarek/pull/515) - Improve preprocessing subworkflows
- [#474](https://github.com/nf-core/sarek/pull/474), [#475](https://github.com/nf-core/sarek/pull/475) - Sync `TEMPLATE` with `tools` `2.2`
- [#487](https://github.com/nf-core/sarek/pull/487), [#489](https://github.com/nf-core/sarek/pull/489), [#492](https://github.com/nf-core/sarek/pull/492), [#497](https://github.com/nf-core/sarek/pull/497), [#522](https://github.com/nf-core/sarek/pull/522), [#583](https://github.com/nf-core/sarek/pull/583) - Improve variant calling subworkflows
- [#498](https://github.com/nf-core/sarek/pull/498) - Update docs
- [#501](https://github.com/nf-core/sarek/pull/501) - Sync `TEMPLATE` with `tools` `2.3`
- [#511](https://github.com/nf-core/sarek/pull/511) - Sync `TEMPLATE` with `tools` `2.3.2`
- [#520](https://github.com/nf-core/sarek/pull/520) - Improve annotation subworkflows
- [#537](https://github.com/nf-core/sarek/pull/537) - Update workflow figure
- [#539](https://github.com/nf-core/sarek/pull/539) - Update `CITATIONS.md`
- [#544](https://github.com/nf-core/sarek/pull/544) - `Mutect2` is no longer compatible with `--no_intervals`
- [#551](https://github.com/nf-core/sarek/pull/551) - Sync `TEMPLATE` with `tools` `2.4`
- [#562](https://github.com/nf-core/sarek/pull/562) - Restart from `--step annotate` is now also requiring a CSV file
- [#563](https://github.com/nf-core/sarek/pull/563) - Updated subway map
- [#570](https://github.com/nf-core/sarek/pull/570) - Extract mpileup into its own subworkflow; zip mpileup files
- [#571](https://github.com/nf-core/sarek/pull/571) - Including and using GATK4's mergeVcfs
- [#572](https://github.com/nf-core/sarek/pull/572) - Adjusted subway map svg for firefox compatibility
- [#577](https://github.com/nf-core/sarek/pull/577) - Update `RELEASE_CHECKLIST`
- [#578](https://github.com/nf-core/sarek/pull/578) - Updated module deeptools/bamcoverage
- [#585](https://github.com/nf-core/sarek/pull/585) - Remove explicit BAM to CRAM conversion after MarkduplicatesSpark; tool does it internally
- [#581](https://github.com/nf-core/sarek/pull/581) - `TIDDIT` is updated to `3.1.0`
- [#593](https://github.com/nf-core/sarek/pull/593) - update `ensembl-vep` cache version and module
- [#600](https://github.com/nf-core/sarek/pull/600) - Remove `TODO` in awsfulltest
- [#606](https://github.com/nf-core/sarek/pull/606) - Updated `ASCAT` to version `3.0` as module
- [#608](https://github.com/nf-core/sarek/pull/608) - Prevent candidate VCFs from getting published in manta
- [#618](https://github.com/nf-core/sarek/pull/618) - Update `multiqc` module
- [#618](https://github.com/nf-core/sarek/pull/618) - Update test yml files
- [#620](https://github.com/nf-core/sarek/pull/620) - `gender` is now `sex` in the samplesheet
- [#630](https://github.com/nf-core/sarek/pull/630) - Update citations file
- [#632](https://github.com/nf-core/sarek/pull/632) - Update `snpEff` version to `5.1` and cache up to `105`
- [#632](https://github.com/nf-core/sarek/pull/632) - Update `VEP` version to `106.1` and cache up to `106`
- [#618](https://github.com/nf-core/sarek/pull/618) - Update `multiqc` module update test yml files
- [#618](https://github.com/nf-core/sarek/pull/618) - Update test yml files
- [#633](https://github.com/nf-core/sarek/pull/633) - Update `BCFTOOLS` version to `1.15.1`
- [#644](https://github.com/nf-core/sarek/pull/644) - Use `-Y` for `bwa-mem(2)` and remove `-M`
- [#645](https://github.com/nf-core/sarek/pull/645) - Merge `tests/nextflow.config` in `conf/test.config`
- [#646](https://github.com/nf-core/sarek/pull/646) - Update `nextflow_schema.json` to reflect new parameters and functions, removes `--annotation_cache`, removes `--ascat_chromosomes`
- [#649](https://github.com/nf-core/sarek/pull/649) - Update, simplify and add more files to all `test_*.yml` files
- [#651](https://github.com/nf-core/sarek/pull/651) - Added TIDDIT_SOMATIC subworkflow
- [#653](https://github.com/nf-core/sarek/pull/653) - Coherent results subfolder structure between preprocessing, variantcalling and reporting
- [#659](https://github.com/nf-core/sarek/pull/659) - Update usage.md docu section on `How to run ASCAT with WES`
- [#661](https://github.com/nf-core/sarek/pull/661) - Add cnvkit reference creation to index subway map
- [#662](https://github.com/nf-core/sarek/pull/662) - Add bgzipped and indexed GATKBundle reference files for `GATK.GRCh37` and replace germline-resources with GATKBundle one
- [#663](https://github.com/nf-core/sarek/pull/663) - Add separate parameters for `ASCAT` and `ControlFREEC` back in
- [#668](https://github.com/nf-core/sarek/pull/668) - Update annotation documentation
- [#674](https://github.com/nf-core/sarek/pull/664) - Default value for splitting is `50000000`

### Fixed

- [#234](https://github.com/nf-core/sarek/pull/234) - Switching to DSL2
- [#234](https://github.com/nf-core/sarek/pull/234), [#238](https://github.com/nf-core/sarek/pull/238) - Add modules and sub workflow for building indices
- [#234](https://github.com/nf-core/sarek/pull/234), [#252](https://github.com/nf-core/sarek/pull/252), [#256](https://github.com/nf-core/sarek/pull/256), [#283](https://github.com/nf-core/sarek/pull/283), [#334](https://github.com/nf-core/sarek/pull/334) - Update Nextflow `19.10.0` -> `20.11.0-edge`
- [#239](https://github.com/nf-core/sarek/pull/239) - Restore Sarek ascii art to header
- [#241](https://github.com/nf-core/sarek/pull/241), [#248](https://github.com/nf-core/sarek/pull/248), [#250](https://github.com/nf-core/sarek/pull/250), [#257](https://github.com/nf-core/sarek/pull/257), [#259](https://github.com/nf-core/sarek/pull/259) - Add modules and sub workflow for preprocessing
- [#242](https://github.com/nf-core/sarek/pull/242), [#244](https://github.com/nf-core/sarek/pull/244), [#245](https://github.com/nf-core/sarek/pull/245), [#246](https://github.com/nf-core/sarek/pull/246), [#247](https://github.com/nf-core/sarek/pull/247), [#249](https://github.com/nf-core/sarek/pull/249), [#252](https://github.com/nf-core/sarek/pull/252), [#256](https://github.com/nf-core/sarek/pull/256), [#263](https://github.com/nf-core/sarek/pull/263), [#264](https://github.com/nf-core/sarek/pull/264), [#283](https://github.com/nf-core/sarek/pull/283), [#285](https://github.com/nf-core/sarek/pull/285), [#338](https://github.com/nf-core/sarek/pull/338) - Refactor `dsl2` branch
- [#257](https://github.com/nf-core/sarek/pull/257) - Use a params modules config file
- [#266](https://github.com/nf-core/sarek/pull/266), [#285](https://github.com/nf-core/sarek/pull/285), [#297](https://github.com/nf-core/sarek/pull/297) - Add modules and sub workflow for variant calling
- [#333](https://github.com/nf-core/sarek/pull/333) - Bump `Sarek` version to `3.0dev`
- [#334](https://github.com/nf-core/sarek/pull/334) - Sync `dsl2` and `dev` branches
- [#342](https://github.com/nf-core/sarek/pull/342) - Update `README.md`
- [#386](https://github.com/nf-core/sarek/pull/386) - Annotation is back
- [#410](https://github.com/nf-core/sarek/pull/410), [#412](https://github.com/nf-core/sarek/pull/412), [#584](https://github.com/nf-core/sarek/pull/584) - Update `CI` tests
- [#418](https://github.com/nf-core/sarek/pull/418) - Fix `known_sites` channels
- [#432](https://github.com/nf-core/sarek/pull/432), [#457](https://github.com/nf-core/sarek/pull/457) - Sort before `tabix index`
- [#454](https://github.com/nf-core/sarek/pull/454) - Input is optional (can actually be found automatically by `Sarek` if previously run)
- [#463](https://github.com/nf-core/sarek/pull/463), [#468](https://github.com/nf-core/sarek/pull/468) - Fix `nf-core lint`
- [#513](https://github.com/nf-core/sarek/pull/513), [#527](https://github.com/nf-core/sarek/pull/527) - CNV is back
- [#529](https://github.com/nf-core/sarek/pull/529) - Do not save `versions.yml` files
- [#524](https://github.com/nf-core/sarek/pull/524) - Fix intervals usage by counting the actual list of scatter/gather files produced and not overall number of intervals
- [#549](https://github.com/nf-core/sarek/pull/549) - Fix unique lanes required for Freebayes: issue [#311](https://github.com/nf-core/sarek/issues/311), replaces `meta.clone()` with actual copy of map to avoid issues with <https://nfcore.slack.com/archives/C027CM7P08M/p1644241819942339>
- [#567](https://github.com/nf-core/sarek/pull/567) - Fix interval name resolving during scatter/gather by moving logic to modules.config causing name to be correctly resolved on process execution; also fixed duplicate naming when variant callers produce multiple vcf files by adding field `type` to `meta` map
- [#585](https://github.com/nf-core/sarek/pull/585) - Fix Spark usage for GATK4 modules
- [#587](https://github.com/nf-core/sarek/pull/587) - Fix issue with VEP extra files
- [#581](https://github.com/nf-core/sarek/pull/581) - `TIDDIT` is back
- [#590](https://github.com/nf-core/sarek/pull/590) - Fix empty folders during scatter/gather
- [#592](https://github.com/nf-core/sarek/pull/592) - Fix optional resources for Mutect2, GetPileupSummaries, and HaplotypeCaller: issue [#299](https://github.com/nf-core/sarek/issues/299), [#359](https://github.com/nf-core/sarek/issues/359), [#367](https://github.com/nf-core/sarek/issues/367)
- [#598](https://github.com/nf-core/sarek/pull/598), [#614](https://github.com/nf-core/sarek/pull/614), [#626](https://github.com/nf-core/sarek/pull/626) - Remove WARNING message for config selector not matching
- [#599](https://github.com/nf-core/sarek/pull/599) - Add checks for correct data type for `params.step`
- [#599](https://github.com/nf-core/sarek/pull/599) - Add checks for no empty `--tools` with `--step variant_calling` or `--step annotate`
- [#600](https://github.com/nf-core/sarek/pull/600) - Remove `nf-core lint` warnings
- [#602](https://github.com/nf-core/sarek/pull/602) - Fixed bug in `alignment_to_fastq` and added tests
- [#609](https://github.com/nf-core/sarek/pull/609) - Remove unused intervals code, reorganize combined intervals file
- [#613](https://github.com/nf-core/sarek/pull/613) - Fixed filenames for `dbnsfp` and `SpliceAI` `VEP` plugin
- [#615](https://github.com/nf-core/sarek/pull/615) - Fix ASCAT igenomes file paths
- [#619](https://github.com/nf-core/sarek/pull/619) - Fix issue with checking samplesheet content with AWS
- [#628](https://github.com/nf-core/sarek/pull/628) - Fix issue with value converting to string before schema validation
- [#628](https://github.com/nf-core/sarek/pull/628) - Fix dbsnp check issue with `--step annotate`
- [#618](https://github.com/nf-core/sarek/pull/618) - Fix `bcftools/vcftools` sample labelling in multiqc report
- [#618](https://github.com/nf-core/sarek/pull/618) - Fix issue with tiddit [#621](https://github.com/nf-core/sarek/issues/621)
- [#618](https://github.com/nf-core/sarek/pull/618) - Fix channel issue with `targets.bed` in prepare_intervals
- [#634](https://github.com/nf-core/sarek/pull/634) - Fix issue with samtools/mosdepth plots in multiqc_report
- [#641](https://github.com/nf-core/sarek/pull/641) - Fix issue with duplicate substring in tools and skip_tools
- [#642](https://github.com/nf-core/sarek/pull/642) - Only unzip ref files if tool is run, only publish ref files if `--save_reference` and simplify CNKit logic
- [#650](https://github.com/nf-core/sarek/pull/650) - Fix intervals checks
- [#654](https://github.com/nf-core/sarek/pull/654) - Allow any step but annotation to start from BAM files
- [#655](https://github.com/nf-core/sarek/pull/655) - Fix `--intervals false` logic & add versioning for local modules
- [#658](https://github.com/nf-core/sarek/pull/658) - Fix split fastq names in multiqc-report
- [#666](https://github.com/nf-core/sarek/pull/666) - Simplify multiqc config channel input
- [#668](https://github.com/nf-core/sarek/pull/668) - Add `snpeff_version` and `vep_version` to `schema_ignore_params` to avoid issue when specifying on command line
- [#669](https://github.com/nf-core/sarek/pull/669) - Fix path to files when creating csv files

### Dependencies

| Dependency             | Old version | New version |
| ---------------------- | ----------- | ----------- |
| `ascat`                | 2.5.2       | 3.0.0       |
| `bcftools`             | 1.9         | 1.15.1      |
| `bwa-mem2`             | 2.0         | 2.2.1       |
| `bwa`                  | 0.7.17      | unchanged   |
| `cancerit-allelecount` | 4.0.2       | 4.3.0       |
| `cnvkit`               | 0.9.6       | 0.9.9       |
| `control-freec`        | 11.6        | unchanged   |
| `deepvariant`          | added       | 1.3.0       |
| `dragmap`              | added       | 1.2.1       |
| `ensembl-vep`          | 99.2        | 106.1       |
| `fastp`                | added       | 0.23.2      |
| `fastqc`               | 0.11.9      | unchanged   |
| `fgbio`                | 1.1.0       | 2.0.2       |
| `freebayes`            | 1.3.2       | 1.3.5       |
| `gatk4`                | 4.1.7.0     | 4.2.6.1     |
| `gawk`                 | added       | 5.1.0       |
| `genesplicer`          | 1.0         | removed     |
| `htslib`               | 1.9         | removed     |
| `llvm-openmp`          | 8.0.1       | removed     |
| `manta`                | 1.6.0       | unchanged   |
| `markdown`             | 3.1.1       | removed     |
| `mosdepth`             | 0.3.3       | unchanged   |
| `msisensor-pro`        | 1.1.a       | 1.2.0       |
| `msisensor`            | 0.5         | removed     |
| `multiqc`              | 1.8         | 1.13a       |
| `openjdk`              | added       | 8.0.312     |
| `openmp`               | 8.0.1       | removed     |
| `p7zip`                | added       | 15.09       |
| `pigz`                 | 2.3.4       | unchanged   |
| `pygments`             | 2.5.2       | removed     |
| `pymdown-extensions`   | 6.0         | removed     |
| `qualimap`             | 2.2.2d      | removed     |
| `r-ggplot2`            | 3.3.0       | removed     |
| `samblaster`           | 0.1.24      | 0.1.26      |
| `samtools`             | 1.9         | 1.15.1      |
| `sed`                  | added       | 4.7         |
| `snpeff`               | 4.3.1t      | 5.1         |
| `strelka`              | 2.9.10      | unchanged   |
| `svdb`                 | added       | 2.6.1       |
| `tabix`                | added       | 1.11        |
| `tiddit`               | 2.7.1       | 3.1.0       |
| `trim-galore`          | 0.6.5       | removed     |
| `vcfanno`              | 0.3.2       | removed     |
| `vcftools`             | 0.1.16      | unchanged   |

### Deprecated

### Removed

- [#485](https://github.com/nf-core/sarek/pull/485) - `--skip_qc`, `--skip_markduplicates` and `--skip_bqsr` is now `--skip_tools`
- [#538](https://github.com/nf-core/sarek/pull/538) - `--sequencing_center` is now `--seq_center`
- [#538](https://github.com/nf-core/sarek/pull/538) - `--markdup_java_options` has been removed
- [#539](https://github.com/nf-core/sarek/pull/539) - `--annotate_tools` has been removed
- [#539](https://github.com/nf-core/sarek/pull/539) - `--cadd_cache`, `--cadd_indels`, `--cadd_indels_tbi`, `--cadd_wg_snvs`, `--cadd_wg_snvs_tbi` have been removed
- [#539](https://github.com/nf-core/sarek/pull/539) - `--genesplicer` has been removed
- [#539](https://github.com/nf-core/sarek/pull/539) - `conf/genomes.config` and `params.genomes_base` have been removed
- [#562](https://github.com/nf-core/sarek/pull/562) - Restart from `--step annotate` from folder is removed. Use a `csv` file instead
- [#571](https://github.com/nf-core/sarek/pull/571) - Removed the local module `concat_vcf`
- [#605](https://github.com/nf-core/sarek/pull/605) - Removed Scatter/gather from GATK_SINGLE_SAMPLE_GERMLINE_VARIANT_CALLING, all intervals are processed together
- [#643](https://github.com/nf-core/sarek/pull/643) - Removed Sentieon parameters

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
- [#195](https://github.com/nf-core/sarek/pull/195), [#202](https://github.com/nf-core/sarek/pull/202) - Add `--save_bam_mapped` params to save mapped BAMs
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

- [#48](https://github.com/nf-core/sarek/issues/48) - Fix `singularity.autoMounts` issue
- [#49](https://github.com/nf-core/sarek/issues/49) - Use correct tag for annotation containers
- [#50](https://github.com/nf-core/sarek/issues/50) - Fix paths for scripts

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
- [#9](https://github.com/nf-core/sarek/pull/9) - Fastq files are named following "${idRun}\_R1.fastq.gz" in the `FastQC` output for easier reporting
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
- [#23](https://github.com/nf-core/sarek/pull/23) - Rename `genomeFile`, `genomeIndex` and `genomeDict` by `fasta`, `fastaFai` and `dict`
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
- [#673](https://github.com/SciLifeLab/Sarek/pull/673), [#675](https://github.com/SciLifeLab/Sarek/pull/675), [#676](https://github.com/SciLifeLab/Sarek/pull/676) - Profiles for BinAC and CFC clusters in Tübingen
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
- [#672](https://github.com/SciLifeLab/Sarek/pull/672) - Process `PullSingularityContainers` from `buildContainers.nf` now expect a file with the correct `.simg` extension for singularity images, and no longer the `.img` one
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
- [#574](https://github.com/SciLifeLab/Sarek/pull/574), [#580](https://github.com/SciLifeLab/Sarek/pull/580) - Abstracts for `NPMI`, `JOBIM` and `EACR25`
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
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - `extractBams()` now takes an extra parameter
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Move `checkFileExtension()`, `checkParameterExistence()`, `checkParameterList()`, `checkReferenceMap()`, `checkRefExistence()`, `extractBams()`, `extractGenders()`, `returnFile()`, `returnStatus()` and `returnTSV()` functions to `SarekUtils`
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Reduce data footprint for Process `CreateRecalibrationTable`
- [#597](https://github.com/SciLifeLab/Sarek/pull/597) - Replace deprecated operator `phase` by `join`
- [#599](https://github.com/SciLifeLab/Sarek/pull/599) - Merge is tested with `ANNOTATEALL`
- [#604](https://github.com/SciLifeLab/Sarek/pull/604) - Synching `GRCh38` `wgs_calling_regions` bedfiles
- [#607](https://github.com/SciLifeLab/Sarek/pull/607) - One container approach
- [#607](https://github.com/SciLifeLab/Sarek/pull/607) - Update to `GATK4`
- [#608](https://github.com/SciLifeLab/Sarek/pull/608) - Update `Nextflow` required version
- [#616](https://github.com/SciLifeLab/Sarek/pull/616) - Update `CHANGELOG`
- [#617](https://github.com/SciLifeLab/Sarek/pull/617) - Replace deprecated ` Nextflow ``$name ` syntax with `withName`

### Fixed

- [#560](https://github.com/SciLifeLab/Sarek/pull/560) - Display message for `repository` and `containerPath`
- [#566](https://github.com/SciLifeLab/Sarek/pull/566) - `slurmDownload` profile
- [#579](https://github.com/SciLifeLab/Sarek/pull/579), [#584](https://github.com/SciLifeLab/Sarek/pull/584) - `Manta` output reorganized after modification for `Strelka Best Practices` process
- [#585](https://github.com/SciLifeLab/Sarek/pull/583) - Trace file is plain txt
- [#590](https://github.com/SciLifeLab/Sarek/pull/590), [#593](https://github.com/SciLifeLab/Sarek/pull/593) - Fix `Singularity` installation in `Travis CI` testing
- [#598](https://github.com/SciLifeLab/Sarek/pull/598), [#601](https://github.com/SciLifeLab/Sarek/pull/601) - Fixes for `Python` script `selectROI.py` to work with `CLC` viewer

### Removed

- [#607](https://github.com/SciLifeLab/Sarek/pull/607) - Remove `Mutect1`

## [2.0.0](https://github.com/SciLifeLab/Sarek/releases/tag/2.0.0) - 2018-03-23

First release under the `Sarek` name, from the National Park in Northern Sweden.

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
