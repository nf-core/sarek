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
        params.pon,
        params.pon_tbi,
        params.aligner,
        params.step,
        params.tools,
        params.vep_include_fasta,
    )

    // TODO: add a params for msisensorpro_scan
    msisensorpro_scan = PREPARE_GENOME.out.msisensorpro_scan

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
        PREPARE_GENOME.out.ascat_alleles,
        PREPARE_GENOME.out.ascat_loci,
        PREPARE_GENOME.out.ascat_loci_gc,
        PREPARE_GENOME.out.ascat_loci_rt,
        params.aligner,
        PREPARE_GENOME.out.bcftools_annotations,
        PREPARE_GENOME.out.bcftools_annotations_tbi,
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
        msisensorpro_scan,
        params.ngscheckmate_bed ? Channel.value(params.ngscheckmate_bed) : Channel.empty(),
        PREPARE_GENOME.out.pon,
        PREPARE_GENOME.out.pon_tbi,
        params.sentieon_dnascope_model ? Channel.fromPath(params.sentieon_dnascope_model).collect() : Channel.value([]),
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
    multiqc_report = SAREK.out.multiqc_report // channel: /path/to/multiqc_report.html
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
