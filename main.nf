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
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK                           } from './workflows/sarek'
include { ANNOTATION_CACHE_INITIALISATION } from './subworkflows/local/annotation_cache_initialisation'
include { DOWNLOAD_CACHE_SNPEFF_VEP       } from './subworkflows/local/download_cache_snpeff_vep'
include { PIPELINE_COMPLETION             } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PIPELINE_INITIALISATION         } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PREPARE_INTERVALS               } from './subworkflows/local/prepare_intervals'
include { PREPARE_REFERENCE_CNVKIT        } from './subworkflows/local/prepare_reference_cnvkit'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.references,
    )

    // WORKFLOW: Run main workflow
    NFCORE_SAREK(PIPELINE_INITIALISATION.out.samplesheet, PIPELINE_INITIALISATION.out.references, params.aligner)

    // SUBWORKFLOW: Run completion tasks
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_SAREK.out.multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Run main nf-core/sarek analysis pipeline
workflow NFCORE_SAREK {
    take:
    samplesheet
    references
    aligner

    main:
    versions = Channel.empty()

    // reduce the references map to keep the minimal set of keys
    def reduce_references_meta = { meta -> meta.subMap(['id']) }

    fasta = references.map { meta, fasta -> [reduce_references_meta(meta), file(fasta)] }.collect()
    fasta_dict = references.map { meta, _fasta -> meta.fasta_dict ? [reduce_references_meta(meta), file(meta.fasta_dict)] : null }.collect()
    fasta_fai = references.map { meta, _fasta -> meta.fasta_fai ? [reduce_references_meta(meta), file(meta.fasta_fai)] : null }.collect()
    bwamem1_index = references.map { meta, _fasta -> meta.bwamem1_index ? [reduce_references_meta(meta), file(meta.bwamem1_index)] : null }.collect()
    bwamem2_index = references.map { meta, _fasta -> meta.bwamem2_index ? [reduce_references_meta(meta), file(meta.bwamem2_index)] : null }.collect()
    dragmap_hashtable = references.map { meta, _fasta -> meta.dragmap_hashtable ? [reduce_references_meta(meta), file(meta.dragmap_hashtable)] : null }.collect()
    intervals_bed = references.map { meta, _fasta -> meta.intervals_bed ? [reduce_references_meta(meta), file(meta.intervals_bed)] : null }.collect()
    dbsnp = references.map { meta, _fasta -> meta.vcf.dbsnp.vcf ? [reduce_references_meta(meta), file(meta.vcf.dbsnp.vcf)] : null }.transpose().collect()
    dbsnp_tbi = references.map { meta, _fasta -> meta.vcf.dbsnp.vcf_tbi ? [reduce_references_meta(meta), file(meta.vcf.dbsnp.vcf_tbi)] : null }.transpose().collect()
    germline_resource = references.map { meta, _fasta -> meta.vcf.germline_resource.vcf ? [reduce_references_meta(meta), file(meta.vcf.germline_resource.vcf)] : null }.collect()
    germline_resource_tbi = references.map { meta, _fasta -> meta.vcf.germline_resource.vcf_tbi ? [reduce_references_meta(meta), file(meta.vcf.germline_resource.vcf_tbi)] : null }.collect()
    known_indels = references.map { meta, _fasta -> meta.vcf.known_indels.vcf ? [reduce_references_meta(meta), file(meta.vcf.known_indels.vcf)] : null }.collect()
    known_indels_tbi = references.map { meta, _fasta -> meta.vcf.known_indels.vcf_tbi ? [reduce_references_meta(meta), file(meta.vcf.known_indels.vcf_tbi)] : null }.collect()
    known_snps = references.map { meta, _fasta -> meta.vcf.known_snps.vcf ? [reduce_references_meta(meta), file(meta.vcf.known_snps.vcf)] : null }.collect()
    known_snps_tbi = references.map { meta, _fasta -> meta.vcf.known_snps.vcf_tbi ? [reduce_references_meta(meta), file(meta.vcf.known_snps.vcf_tbi)] : null }.collect()
    pon = references.map { meta, _fasta -> meta.vcf.pon.vcf ? [reduce_references_meta(meta), file(meta.vcf.pon.vcf)] : null }.collect()
    pon_tbi = references.map { meta, _fasta -> meta.vcf.pon.vcf_tbi ? [reduce_references_meta(meta), file(meta.vcf.pon.vcf_tbi)] : null }.collect()

    // Gather index for mapping given the chosen aligner
    index_alignment = aligner == "bwa-mem" || aligner == "sentieon-bwamem"
        ? bwamem1_index
        : aligner == "bwa-mem2"
            ? bwamem2_index
            : dragmap_hashtable

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels = dbsnp.mix(known_indels).groupTuple().collect()
    known_sites_indels_tbi = dbsnp_tbi.mix(known_indels_tbi).groupTuple().collect()
    known_sites_snps = dbsnp.mix(known_snps).groupTuple().collect()
    known_sites_snps_tbi = dbsnp_tbi.mix(known_snps_tbi).groupTuple().collect()

    // Initialize value channels based on params for annotation

    vep_extra_files = []

    if (params.dbnsfp && params.dbnsfp_tbi) {
        vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
        vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
    }

    if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
        vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
    }

    // Prepare intervals for spread/gather (split and tabix index)
    // This depends on params.nucleotides_per_second so it needs to be done on runtime
    PREPARE_INTERVALS(intervals_bed, params.no_intervals, params.nucleotides_per_second, params.outdir, params.step)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined
    intervals_bed_combined_for_variant_calling = PREPARE_INTERVALS.out.intervals_bed_combined

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes
        ? intervals_bed_combined.map { it -> [[id: it.baseName], it] }.collect()
        : Channel.value([[id: 'null'], []])
    intervals = PREPARE_INTERVALS.out.intervals_bed
    // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi
    // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather
    intervals_and_num_intervals = intervals.map { intervals_, num_intervals ->
        if (num_intervals < 1) {
            [[], num_intervals]
        }
        else {
            [intervals_, num_intervals]
        }
    }
    intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map { intervals_, num_intervals ->
        if (num_intervals < 1) {
            [[], [], num_intervals]
        }
        else {
            [intervals_[0], intervals_[1], num_intervals]
        }
    }
    if (params.tools && params.tools.split(',').contains('cnvkit')) {
        if (params.cnvkit_reference) {
            cnvkit_reference = Channel.fromPath(params.cnvkit_reference).collect()
        }
        else {
            PREPARE_REFERENCE_CNVKIT(fasta, intervals_bed_combined)
            cnvkit_reference = PREPARE_REFERENCE_CNVKIT.out.cnvkit_reference
            versions = versions.mix(PREPARE_REFERENCE_CNVKIT.out.versions)
        }
    }
    else {
        cnvkit_reference = Channel.value([])
    }

    // Gather used softwares versions
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    vep_fasta = params.vep_include_fasta ? fasta.map { fasta_ -> [[id: fasta_.baseName], fasta_] } : [[id: 'null'], []]

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

    //
    // WORKFLOW: Run pipeline
    //
    SAREK(
        samplesheet,
        [],
        [],
        [],
        [],
        [],
        [],
        cnvkit_reference,
        dbsnp,
        dbsnp_tbi,
        "",
        fasta_dict,
        fasta,
        fasta_fai,
        [],
        germline_resource,
        germline_resource_tbi,
        index_alignment,
        intervals_and_num_intervals,
        intervals_bed_combined,
        intervals_bed_combined_for_variant_calling,
        intervals_bed_gz_tbi_and_num_intervals,
        intervals_bed_gz_tbi_combined,
        intervals_for_preprocessing,
        "",
        known_sites_indels,
        known_sites_indels_tbi,
        known_sites_snps,
        known_sites_snps_tbi,
        "",
        [],
        [],
        [],
        [],
        pon,
        pon_tbi,
        [],
        [],
        snpeff_cache,
        vep_cache,
        "",
        vep_extra_files,
        vep_fasta,
        "",
        "",
    )

    emit:
    multiqc_report = SAREK.out.multiqc_report // channel: /path/to/multiqc_report.html
}
