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
        args,
        params.outdir,
        params.input,
        params.references,
        params.step,
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

    // References' files from the references yaml or params
    ascat_alleles           = extract_references_file(references, params.ascat_alleles, 'ascat_alleles').collect()
    ascat_loci              = extract_references_file(references, params.ascat_loci, 'ascat_loci').collect()
    ascat_loci_gc           = extract_references_file(references, params.ascat_loci_gc, 'ascat_loci_gc').collect()
    ascat_loci_rt           = extract_references_file(references, params.ascat_loci_rt, 'ascat_loci_rt').collect()
    bwamem1_index           = extract_references_file(references, params.bwa, 'bwamem1_index').collect()
    bwamem2_index           = extract_references_file(references, params.bwamem2, 'bwamem2_index').collect()
    cf_chrom_len            = extract_references_file(references, params.cf_chrom_len, 'cf_chrom_len').collect()
    chr_dir                 = extract_references_file(references, params.chr_dir, 'chr_dir').collect()
    dragmap_hashtable       = extract_references_file(references, params.dragmap, 'dragmap_hashtable').collect()
    fasta                   = extract_references_file(references, params.fasta, 'fasta').collect()
    fasta_dict              = extract_references_file(references, params.dict, 'fasta_dict').collect()
    fasta_fai               = extract_references_file(references, params.fasta_fai, 'fasta_fai').collect()
    intervals_bed           = extract_references_file(references, params.intervals, 'intervals_bed').collect()
    mappability             = extract_references_file(references, params.mappability, 'mappability').collect()
    msisensorpro_scan       = extract_references_file(references, params.msisensorpro_scan, 'msisensorpro_scan').collect()
    ngscheckmate_bed        = extract_references_file(references, params.ngscheckmate_bed, 'ngscheckmate_bed').collect()
    sentieon_dnascope_model = extract_references_file(references, params.sentieon_dnascope_model, 'sentieon_dnascope_model').collect()

    // References' VCFs and related from the references yaml or params
    dbsnp                   = extract_references_vcf(references, params.dbsnp, 'dbsnp', 'vcf').collect()
    dbsnp_tbi               = extract_references_vcf(references, params.dbsnp_tbi, 'dbsnp', 'vcf_tbi').collect()
    dbsnp_vqsr              = extract_references_vcf(references, params.dbsnp_vqsr, 'dbsnp', 'vcf_vqsr').collect()
    germline_resource       = extract_references_vcf(references, params.germline_resource, 'germline_resource', 'vcf').collect()
    germline_resource_tbi   = extract_references_vcf(references, params.germline_resource_tbi, 'germline_resource', 'vcf_tbi').collect()
    known_indels            = extract_references_vcf(references, params.known_indels, 'known_indels', 'vcf').collect()
    known_indels_tbi        = extract_references_vcf(references, params.known_indels_tbi, 'known_indels', 'vcf_tbi').collect()
    known_indels_vqsr       = extract_references_vcf(references, params.known_indels_vqsr, 'known_indels', 'vcf_vqsr').collect()
    known_snps              = extract_references_vcf(references, params.known_snps, 'known_snps', 'vcf').collect()
    known_snps_tbi          = extract_references_vcf(references, params.known_snps_tbi, 'known_snps', 'vcf_tbi').collect()
    known_snps_vqsr         = extract_references_vcf(references, params.known_snps_vqsr, 'known_snps', 'vcf_vqsr').collect()
    pon                     = extract_references_vcf(references, params.pon, 'pon', 'vcf').collect()
    pon_tbi                 = extract_references_vcf(references, params.pon_tbi, 'pon', 'vcf_tbi').collect()

    // References' values from the references yaml or params
    ascat_genome            = extract_references_value(references, params.ascat_genome, 'ascat_genome').collect()
    snpeff_db               = extract_references_value(references, params.snpeff_db, 'snpeff_db').collect()
    vep_cache_version       = extract_references_value(references, params.vep_cache_version, 'vep_cache_version').collect()
    vep_genome              = extract_references_value(references, params.vep_genome, 'vep_genome').collect()
    vep_species             = extract_references_value(references, params.vep_species, 'vep_species').collect()

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels      = dbsnp.mix(known_indels).groupTuple().collect()
    known_sites_indels_tbi  = dbsnp_tbi.mix(known_indels_tbi).groupTuple().collect()
    known_sites_snps        = dbsnp.mix(known_snps).groupTuple().collect()
    known_sites_snps_tbi    = dbsnp_tbi.mix(known_snps_tbi).groupTuple().collect()

    // Gather index for mapping given the chosen aligner
    index_alignment = aligner == "bwa-mem" || aligner == "sentieon-bwamem"
        ? bwamem1_index
        : aligner == "bwa-mem2"
            ? bwamem2_index
            : dragmap_hashtable

    bcftools_annotations        = params.bcftools_annotations ? Channel.fromPath(params.bcftools_annotations).collect() : Channel.value([])
    bcftools_annotations_tbi    = params.bcftools_annotations ? params.bcftools_annotations_tbi ? Channel.fromPath(params.bcftools_annotations_tbi).collect() : Channel.value([]) : Channel.value([])
    bcftools_header_lines       = params.bcftools_header_lines ?: Channel.value([])

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

    // Build intervals if needed
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
    intervals_and_num_intervals = intervals.map { interval, num_intervals ->
        if (num_intervals < 1) {
            [[], num_intervals]
        }
        else {
            [interval, num_intervals]
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

    // Looks for cache information either locally or on the cloud
    ANNOTATION_CACHE_INITIALISATION(
        (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
        params.snpeff_cache,
        snpeff_db,
        (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
        params.vep_cache,
        vep_species,
        vep_cache_version,
        vep_genome,
        params.vep_custom_args,
        "Please refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.",
    )

    snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
    vep_cache = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache

    //
    // WORKFLOW: Run pipeline
    //
    SAREK(
        samplesheet,
        ascat_alleles,
        bcftools_annotations,
        bcftools_annotations_tbi,
        bcftools_header_lines,
        cf_chrom_len,
        chr_dir,
        cnvkit_reference,
        dbsnp,
        dbsnp_tbi,
        dbsnp_vqsr,
        fasta_dict,
        fasta,
        fasta_fai,
        ascat_loci_gc,
        germline_resource,
        germline_resource_tbi,
        index_alignment,
        intervals_and_num_intervals,
        intervals_bed_combined,
        intervals_bed_combined_for_variant_calling,
        intervals_bed_gz_tbi_and_num_intervals,
        intervals_bed_gz_tbi_combined,
        intervals_for_preprocessing,
        known_indels_vqsr,
        known_sites_indels,
        known_sites_indels_tbi,
        known_sites_snps,
        known_sites_snps_tbi,
        known_snps_vqsr,
        ascat_loci,
        mappability,
        msisensorpro_scan,
        ngscheckmate_bed,
        pon,
        pon_tbi,
        ascat_loci_rt,
        sentieon_dnascope_model,
        snpeff_cache,
        snpeff_db,
        vep_cache,
        vep_cache_version,
        vep_extra_files,
        vep_fasta,
        vep_genome,
        vep_species,
    )

    emit:
    multiqc_report = SAREK.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS TO EXTRACT REFERENCES FILES OR VALUES FROM THE REFERENCES YAML OR PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def extract_references_file(references, param, attribute) {
    references.map { meta, _readme ->
        if (param || meta[attribute]) {
            return [meta.subMap(['id']), file(param ?: meta[attribute] ?: [])]
        }
        else {
            return null
        }
    }
}

def extract_references_vcf(references, param, vcf_attribute, attribute) {
    references.map { meta, _readme ->
        if (param || meta.vcf[vcf_attribute][attribute]) {
            return [meta.subMap(['id']), file(param ?: meta.vcf[vcf_attribute][attribute] ?: [])]
        }
        else {
            return null
        }
    }
}

def extract_references_value(references, param, attribute) {
    references.map { meta, _readme ->
        if (param || meta[attribute]) {
            return [meta.subMap(['id']), param ?: meta[attribute] ?: []]
        }
        else {
            return null
        }
    }
}
