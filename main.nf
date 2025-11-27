#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/sarek
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Started March 2016.
    Ported to nf-core May 2019.
    Ported to DSL 2 July 2020.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/sarek:
        An open-source analysis pipeline to detect germline or somatic variants
        from whole genome or targeted sequencing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/sarek
    Website: https://nf-co.re/sarek
    Docs   : https://nf-co.re/sarek/usage
    Slack  : https://nfcore.slack.com/channels/sarek
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.ascat_alleles           = getGenomeAttribute('ascat_alleles')
params.ascat_genome            = getGenomeAttribute('ascat_genome')
params.ascat_loci              = getGenomeAttribute('ascat_loci')
params.ascat_loci_gc           = getGenomeAttribute('ascat_loci_gc')
params.ascat_loci_rt           = getGenomeAttribute('ascat_loci_rt')
params.bwa                     = getGenomeAttribute('bwa')
params.bwamem2                 = getGenomeAttribute('bwamem2')
params.cf_chrom_len            = getGenomeAttribute('cf_chrom_len')
params.chr_dir                 = getGenomeAttribute('chr_dir')
params.dbsnp                   = getGenomeAttribute('dbsnp')
params.dbsnp_tbi               = getGenomeAttribute('dbsnp_tbi')
params.dbsnp_vqsr              = getGenomeAttribute('dbsnp_vqsr')
params.dict                    = getGenomeAttribute('dict')
params.dragmap                 = getGenomeAttribute('dragmap')
params.fasta                   = getGenomeAttribute('fasta')
params.fasta_fai               = getGenomeAttribute('fasta_fai')
params.germline_resource       = getGenomeAttribute('germline_resource')
params.germline_resource_tbi   = getGenomeAttribute('germline_resource_tbi')
params.intervals               = getGenomeAttribute('intervals')
params.known_indels            = getGenomeAttribute('known_indels')
params.known_indels_tbi        = getGenomeAttribute('known_indels_tbi')
params.known_indels_vqsr       = getGenomeAttribute('known_indels_vqsr')
params.known_snps              = getGenomeAttribute('known_snps')
params.known_snps_tbi          = getGenomeAttribute('known_snps_tbi')
params.known_snps_vqsr         = getGenomeAttribute('known_snps_vqsr')
params.mappability             = getGenomeAttribute('mappability')
params.msisensor2_models       = getGenomeAttribute('msisensor2_models')
params.msisensorpro_scan       = getGenomeAttribute('msisensorpro_scan')
params.ngscheckmate_bed        = getGenomeAttribute('ngscheckmate_bed')
params.pon                     = getGenomeAttribute('pon')
params.pon_tbi                 = getGenomeAttribute('pon_tbi')
params.sentieon_dnascope_model = getGenomeAttribute('sentieon_dnascope_model')
params.snpeff_db               = getGenomeAttribute('snpeff_db')
params.vep_cache_version       = getGenomeAttribute('vep_cache_version')
params.vep_genome              = getGenomeAttribute('vep_genome')
params.vep_species             = getGenomeAttribute('vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK                           } from './workflows/sarek'
include { ANNOTATION_CACHE_INITIALISATION } from './subworkflows/local/annotation_cache_initialisation'
include { DOWNLOAD_CACHE_SNPEFF_VEP       } from './subworkflows/local/download_cache_snpeff_vep'
include { PIPELINE_COMPLETION             } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PIPELINE_INITIALISATION         } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PREPARE_GENOME                  } from './subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS               } from './subworkflows/local/prepare_intervals'
include { PREPARE_REFERENCE_CNVKIT        } from './subworkflows/local/prepare_reference_cnvkit'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Run main nf-core/sarek analysis pipeline
workflow NFCORE_SAREK {
    take:
    samplesheet

    main:
    versions = Channel.empty()

    // build indexes if needed
    PREPARE_GENOME(
        params.ascat_alleles,
        params.ascat_loci,
        params.ascat_loci_gc,
        params.ascat_loci_rt,
        params.bbsplit_fasta_list,
        params.bbsplit_index,
        params.bcftools_annotations,
        params.bcftools_annotations_tbi,
        params.bwa,
        params.bwamem2,
        params.chr_dir,
        params.dbsnp,
        params.dbsnp_tbi,
        params.dict,
        params.dragmap,
        params.fasta,
        params.fasta_fai,
        params.germline_resource,
        params.germline_resource_tbi,
        params.known_indels,
        params.known_indels_tbi,
        params.known_snps,
        params.known_snps_tbi,
        params.msisensor2_models,
        params.msisensorpro_scan,
        params.pon,
        params.pon_tbi,
        params.aligner,
        params.step,
        params.tools ?: 'no_tools',
        params.vep_include_fasta,
    )

    // Build intervals if needed
    PREPARE_INTERVALS(PREPARE_GENOME.out.fasta_fai, params.intervals, params.no_intervals, params.nucleotides_per_second, params.outdir, params.step)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined
    intervals_bed_combined_for_variant_calling = PREPARE_INTERVALS.out.intervals_bed_combined

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes
        ? intervals_bed_combined.map { it -> [[id: it.baseName], it] }.collect()
        : Channel.value([[id: 'null'], []])
    // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals = PREPARE_INTERVALS.out.intervals_bed
    // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi

    intervals_and_num_intervals = intervals.map { file, num_intervals ->
        [num_intervals < 1 ? [] : file, num_intervals]
    }

    intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map { file, num_intervals ->
        [num_intervals < 1 ? [] : file[0], num_intervals < 1 ? [] : file[1], num_intervals]
    }

    if (params.tools && params.tools.split(',').contains('cnvkit')) {
        if (params.cnvkit_reference) {
            cnvkit_reference = Channel.fromPath(params.cnvkit_reference).collect()
        }
        else {
            PREPARE_REFERENCE_CNVKIT(PREPARE_GENOME.out.fasta, intervals_bed_combined)
            cnvkit_reference = PREPARE_REFERENCE_CNVKIT.out.cnvkit_reference
            versions = versions.mix(PREPARE_REFERENCE_CNVKIT.out.versions)
        }
    }
    else {
        cnvkit_reference = Channel.value([])
    }
    // Gather used softwares versions
    versions = versions.mix(PREPARE_GENOME.out.versions)
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    // Fails when consensus calling is specified without normalization
    if (params.snv_consensus_calling && !params.normalize_vcfs){
        error("Consensus calling was specified without normalization. Set --normalize_vcfs in addition. See: https://www.biostars.org/p/307035/")
    }


    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, sarek will download the cache
        ensemblvep_info = Channel.of([[id: "${params.vep_cache_version}_${params.vep_genome}"], params.vep_genome, params.vep_species, params.vep_cache_version])
        snpeff_info = Channel.of([[id: "${params.snpeff_db}"], params.snpeff_db])
        DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)
        snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
        vep_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map { _meta, cache -> [cache] }

        versions = versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    }
    else {
        // Looks for cache information either locally or on the cloud
        ANNOTATION_CACHE_INITIALISATION(
            (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            params.snpeff_db,
            (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_custom_args,
            "Please refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.",
        )

        snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
        vep_cache = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }

    vep_extra_files = []

    if (params.dbnsfp && params.dbnsfp_tbi) {
        vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
        vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
    }
    else if (params.dbnsfp && !params.dbnsfp_tbi) {
        System.err.println("DBNSFP: ${params.dbnsfp} has been provided with `--dbnsfp, but no dbnsfp_tbi has")
        System.err.println("cf: https://nf-co.re/sarek/parameters/#dbnsfp")
        error("Execution halted due to dbnsfp inconsistency.")
    }

    if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
        vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
    }

    //
    // WORKFLOW: Run pipeline
    //
    SAREK(
        samplesheet,
        params.aligner,
        params.skip_tools ?: 'no_tools',
        params.step,
        params.tools ?: 'no_tools',
        PREPARE_GENOME.out.ascat_alleles,
        PREPARE_GENOME.out.ascat_loci,
        PREPARE_GENOME.out.ascat_loci_gc,
        PREPARE_GENOME.out.ascat_loci_rt,
        PREPARE_GENOME.out.bbsplit_index,
        PREPARE_GENOME.out.bcftools_annotations,
        PREPARE_GENOME.out.bcftools_annotations_tbi,
        params.bcftools_columns ? Channel.fromPath(params.bcftools_columns).collect() : Channel.value([]),
        params.bcftools_header_lines ? Channel.fromPath(params.bcftools_header_lines).collect() : Channel.empty(),
        params.cf_chrom_len ? Channel.fromPath(params.cf_chrom_len).collect() : [],
        PREPARE_GENOME.out.chr_dir,
        cnvkit_reference,
        PREPARE_GENOME.out.dbsnp,
        PREPARE_GENOME.out.dbsnp_tbi,
        params.dbsnp_vqsr ? Channel.value(params.dbsnp_vqsr) : Channel.empty(),
        PREPARE_GENOME.out.dict,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fasta_fai,
        PREPARE_GENOME.out.germline_resource,
        PREPARE_GENOME.out.germline_resource_tbi,
        PREPARE_GENOME.out.index_alignment,
        intervals_and_num_intervals,
        intervals_bed_combined,
        intervals_bed_combined_for_variant_calling,
        intervals_bed_gz_tbi_and_num_intervals,
        intervals_bed_gz_tbi_combined,
        intervals_for_preprocessing,
        params.known_indels_vqsr ? Channel.value(params.known_indels_vqsr) : Channel.empty(),
        PREPARE_GENOME.out.known_sites_indels,
        PREPARE_GENOME.out.known_sites_indels_tbi,
        PREPARE_GENOME.out.known_sites_snps,
        PREPARE_GENOME.out.known_sites_snps_tbi,
        params.known_snps_vqsr ? Channel.value(params.known_snps_vqsr) : Channel.empty(),
        params.mappability ? Channel.fromPath(params.mappability).collect() : Channel.value([]),
        PREPARE_GENOME.out.msisensor2_models,
        PREPARE_GENOME.out.msisensorpro_scan,
        params.ngscheckmate_bed ? Channel.value(params.ngscheckmate_bed) : Channel.empty(),
        PREPARE_GENOME.out.pon,
        PREPARE_GENOME.out.pon_tbi,
        params.sentieon_dnascope_model ? Channel.fromPath(params.sentieon_dnascope_model).collect() : Channel.value([]),
         // Set Varlociraptor reference files
        params.varlociraptor_scenario_germline   ? Channel.fromPath(params.varlociraptor_scenario_germline).map { it -> [[id: it.baseName - '.yte'], it] }.collect()    : Channel.fromPath("${projectDir}/assets/varlociraptor_germline.yte.yaml").collect(),
        params.varlociraptor_scenario_somatic    ? Channel.fromPath(params.varlociraptor_scenario_somatic).map { it -> [[id: it.baseName - '.yte'], it] }.collect()     : Channel.fromPath("${projectDir}/assets/varlociraptor_somatic.yte.yaml").collect(),
        params.varlociraptor_scenario_tumor_only ? Channel.fromPath(params.varlociraptor_scenario_tumor_only).map { it -> [[id: it.baseName - '.yte'], it] }.collect()  : Channel.fromPath("${projectDir}/assets/varlociraptor_tumor_only.yte.yaml").collect(),
        snpeff_cache,
        params.snpeff_db,
        vep_cache,
        params.vep_cache_version,
        vep_extra_files,
        PREPARE_GENOME.out.vep_fasta,
        params.vep_genome,
        params.vep_species,
        versions,
    )

    emit:
    multiqc_report      = SAREK.out.multiqc_report      // channel: /path/to/multiqc_report.html
    multiqc_data        = SAREK.out.multiqc_data        // channel: /path/to/multiqc_data
    multiqc_plots       = SAREK.out.multiqc_plots       // channel: /path/to/multiqc_plots
    versions            = SAREK.out.versions            // channel: [ path(versions.yml) ]
    reports             = SAREK.out.reports             // channel: [ reports ] - aggregate QC reports
    // Individual QC report channels with meta
    fastqc_zip          = SAREK.out.fastqc_zip          // channel: [ meta, zip ]
    fastqc_html         = SAREK.out.fastqc_html         // channel: [ meta, html ]
    samtools_stats      = SAREK.out.samtools_stats      // channel: [ meta, stats ]
    mosdepth_global     = SAREK.out.mosdepth_global     // channel: [ meta, txt ]
    mosdepth_region     = SAREK.out.mosdepth_region     // channel: [ meta, txt ]
    mosdepth_summary    = SAREK.out.mosdepth_summary    // channel: [ meta, txt ]
    mosdepth_regions_bed = SAREK.out.mosdepth_regions_bed // channel: [ meta, bed.gz ]
    mosdepth_regions_csi = SAREK.out.mosdepth_regions_csi // channel: [ meta, csi ]
    // Markduplicates-stage QC channels
    md_samtools_stats      = SAREK.out.md_samtools_stats      // channel: [ meta, stats ]
    md_mosdepth_global     = SAREK.out.md_mosdepth_global     // channel: [ meta, txt ]
    md_mosdepth_region     = SAREK.out.md_mosdepth_region     // channel: [ meta, txt ]
    md_mosdepth_summary    = SAREK.out.md_mosdepth_summary    // channel: [ meta, txt ]
    md_mosdepth_regions_bed = SAREK.out.md_mosdepth_regions_bed // channel: [ meta, bed.gz ]
    md_mosdepth_regions_csi = SAREK.out.md_mosdepth_regions_csi // channel: [ meta, csi ]
    bcftools_stats      = SAREK.out.bcftools_stats      // channel: [ meta, stats ]
    vcftools_tstv_counts = SAREK.out.vcftools_tstv_counts // channel: [ meta, counts ]
    vcftools_tstv_qual  = SAREK.out.vcftools_tstv_qual  // channel: [ meta, qual ]
    vcftools_filter_summary = SAREK.out.vcftools_filter_summary // channel: [ meta, summary ]
    // Preprocessing outputs
    cram_mapped         = SAREK.out.cram_mapped         // channel: [ meta, cram, crai ]
    bam_mapped          = SAREK.out.bam_mapped          // channel: [ meta, bam, bai ]
    cram_markduplicates = SAREK.out.cram_markduplicates // channel: [ meta, cram, crai ]
    bam_markduplicates  = SAREK.out.bam_markduplicates  // channel: [ meta, bam, bai ]
    cram_recalibrated   = SAREK.out.cram_recalibrated   // channel: [ meta, cram, crai ]
    recal_table         = SAREK.out.recal_table         // channel: [ meta, table ]
    markduplicates_metrics = SAREK.out.markduplicates_metrics // channel: [ meta, metrics ]
    // Variant calling outputs - VCF
    vcf_germline        = SAREK.out.vcf_germline        // channel: [ meta, vcf ]
    vcf_somatic         = SAREK.out.vcf_somatic         // channel: [ meta, vcf ]
    vcf_tumor_only      = SAREK.out.vcf_tumor_only      // channel: [ meta, vcf ]
    // Variant calling outputs - TBI (separate)
    tbi_germline        = SAREK.out.tbi_germline        // channel: [ meta, tbi ]
    tbi_somatic         = SAREK.out.tbi_somatic         // channel: [ meta, tbi ]
    tbi_tumor_only      = SAREK.out.tbi_tumor_only      // channel: [ meta, tbi ]
    // Strelka genome VCF/TBI (separate from variant VCF)
    vcf_strelka_genome  = SAREK.out.vcf_strelka_genome  // channel: [ meta, vcf ]
    tbi_strelka_genome  = SAREK.out.tbi_strelka_genome  // channel: [ meta, tbi ]
    // Annotation outputs
    vcf_annotated       = SAREK.out.vcf_annotated       // channel: [ meta, vcf ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SAREK(PIPELINE_INITIALISATION.out.samplesheet)

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_SAREK.out.multiqc_report,
    )

    publish:
    // QC Reports - MultiQC
    multiqc_report = NFCORE_SAREK.out.multiqc_report
    multiqc_data = NFCORE_SAREK.out.multiqc_data
    multiqc_plots = NFCORE_SAREK.out.multiqc_plots
    // Versions
    versions = NFCORE_SAREK.out.versions
    // Individual QC reports with meta
    fastqc_zip = NFCORE_SAREK.out.fastqc_zip
    fastqc_html = NFCORE_SAREK.out.fastqc_html
    samtools_stats = NFCORE_SAREK.out.samtools_stats
    mosdepth_global = NFCORE_SAREK.out.mosdepth_global
    mosdepth_region = NFCORE_SAREK.out.mosdepth_region
    mosdepth_summary = NFCORE_SAREK.out.mosdepth_summary
    mosdepth_regions_bed = NFCORE_SAREK.out.mosdepth_regions_bed
    mosdepth_regions_csi = NFCORE_SAREK.out.mosdepth_regions_csi
    // Markduplicates-stage QC
    md_samtools_stats = NFCORE_SAREK.out.md_samtools_stats
    md_mosdepth_global = NFCORE_SAREK.out.md_mosdepth_global
    md_mosdepth_region = NFCORE_SAREK.out.md_mosdepth_region
    md_mosdepth_summary = NFCORE_SAREK.out.md_mosdepth_summary
    md_mosdepth_regions_bed = NFCORE_SAREK.out.md_mosdepth_regions_bed
    md_mosdepth_regions_csi = NFCORE_SAREK.out.md_mosdepth_regions_csi
    bcftools_stats = NFCORE_SAREK.out.bcftools_stats
    vcftools_tstv_counts = NFCORE_SAREK.out.vcftools_tstv_counts
    vcftools_tstv_qual = NFCORE_SAREK.out.vcftools_tstv_qual
    vcftools_filter_summary = NFCORE_SAREK.out.vcftools_filter_summary
    // Preprocessing outputs
    cram_mapped = NFCORE_SAREK.out.cram_mapped
    bam_mapped = NFCORE_SAREK.out.bam_mapped
    cram_markduplicates = NFCORE_SAREK.out.cram_markduplicates
    bam_markduplicates = NFCORE_SAREK.out.bam_markduplicates
    cram_recalibrated = NFCORE_SAREK.out.cram_recalibrated
    recal_table = NFCORE_SAREK.out.recal_table
    markduplicates_metrics = NFCORE_SAREK.out.markduplicates_metrics
    // Variant calling outputs - VCF
    vcf_germline = NFCORE_SAREK.out.vcf_germline
    vcf_somatic = NFCORE_SAREK.out.vcf_somatic
    vcf_tumor_only = NFCORE_SAREK.out.vcf_tumor_only
    // Variant calling outputs - TBI (separate)
    tbi_germline = NFCORE_SAREK.out.tbi_germline
    tbi_somatic = NFCORE_SAREK.out.tbi_somatic
    tbi_tumor_only = NFCORE_SAREK.out.tbi_tumor_only
    // Strelka genome VCF/TBI (separate from variant VCF)
    vcf_strelka_genome = NFCORE_SAREK.out.vcf_strelka_genome
    tbi_strelka_genome = NFCORE_SAREK.out.tbi_strelka_genome
    // Annotation outputs
    vcf_annotated = NFCORE_SAREK.out.vcf_annotated
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    // QC Reports - MultiQC
    multiqc_report {
        path 'multiqc'
    }
    multiqc_data {
        path 'multiqc'
    }
    multiqc_plots {
        path 'multiqc'
    }

    // Versions
    versions {
        path 'pipeline_info'
    }

    // Individual QC Reports
    fastqc_zip {
        path { meta, zip -> "reports/fastqc/${meta.id}" }
    }
    fastqc_html {
        path { meta, html -> "reports/fastqc/${meta.id}" }
    }
    samtools_stats {
        path { meta, stats -> "reports/samtools/${meta.id}" }
    }
    mosdepth_global {
        path { meta, txt -> "reports/mosdepth/${meta.id}" }
    }
    mosdepth_region {
        path { meta, txt -> "reports/mosdepth/${meta.id}" }
    }
    mosdepth_summary {
        path { meta, txt -> "reports/mosdepth/${meta.id}" }
    }
    mosdepth_regions_bed {
        path { meta, bed -> "reports/mosdepth/${meta.id}" }
    }
    mosdepth_regions_csi {
        path { meta, csi -> "reports/mosdepth/${meta.id}" }
    }

    // Markduplicates-stage QC Reports
    md_samtools_stats {
        path { meta, stats -> "reports/samtools/${meta.id}" }
    }
    md_mosdepth_global {
        path { meta, txt -> "reports/mosdepth/${meta.id}" }
    }
    md_mosdepth_region {
        path { meta, txt -> "reports/mosdepth/${meta.id}" }
    }
    md_mosdepth_summary {
        path { meta, txt -> "reports/mosdepth/${meta.id}" }
    }
    md_mosdepth_regions_bed {
        path { meta, bed -> "reports/mosdepth/${meta.id}" }
    }
    md_mosdepth_regions_csi {
        path { meta, csi -> "reports/mosdepth/${meta.id}" }
    }

    // TODO: Migrate to reports/bcftools/${meta.variantcaller}/${meta.id} to match vcftools pattern
    bcftools_stats {
        path { meta, stats -> "reports/bcftools/${meta.id}" }
    }
    vcftools_tstv_counts {
        path { meta, counts -> "reports/vcftools/${meta.variantcaller}/${meta.id}" }
    }
    vcftools_tstv_qual {
        path { meta, qual -> "reports/vcftools/${meta.variantcaller}/${meta.id}" }
    }
    vcftools_filter_summary {
        path { meta, summary -> "reports/vcftools/${meta.variantcaller}/${meta.id}" }
    }

    // Preprocessing - Mapped
    cram_mapped {
        path { meta, cram, crai -> "preprocessing/mapped/${meta.id}" }
    }
    bam_mapped {
        path { meta, bam, bai -> "preprocessing/mapped/${meta.id}" }
    }

    // Preprocessing - Markduplicates
    cram_markduplicates {
        path { meta, cram, crai -> "preprocessing/markduplicates/${meta.id}" }
    }
    bam_markduplicates {
        path { meta, bam, bai -> "preprocessing/markduplicates/${meta.id}" }
    }

    // Preprocessing - Recalibrated
    cram_recalibrated {
        path { meta, cram, crai -> "preprocessing/recalibrated/${meta.id}" }
    }

    // Preprocessing - Recal Table
    recal_table {
        path { meta, table -> "preprocessing/recal_table/${meta.id}" }
    }

    // Preprocessing - Markduplicates Metrics
    markduplicates_metrics {
        path { meta, metrics -> "reports/markduplicates/${meta.id}" }
    }

    // Variant Calling - VCF
    vcf_germline {
        path { meta, vcf -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }
    vcf_somatic {
        path { meta, vcf -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }
    vcf_tumor_only {
        path { meta, vcf -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }

    // Variant Calling - TBI (published to same directory as VCF)
    tbi_germline {
        path { meta, tbi -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }
    tbi_somatic {
        path { meta, tbi -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }
    tbi_tumor_only {
        path { meta, tbi -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }

    // Strelka genome VCF/TBI (separate from variant VCF)
    vcf_strelka_genome {
        path { meta, vcf -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }
    tbi_strelka_genome {
        path { meta, tbi -> "variant_calling/${meta.variantcaller}/${meta.id}" }
    }

    // Annotation
    vcf_annotated {
        path { meta, vcf -> "annotation/${meta.variantcaller}/${meta.id}" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[params.genome].containsKey(attribute)) {
            return params.genomes[params.genome][attribute]
        }
    }
    return null
}
