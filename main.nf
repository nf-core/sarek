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

    // References' files from the references asset yaml or params
    ascat_alleles = references.map { meta, _fasta -> meta.ascat_alleles ? [meta.subMap(['id']), file(params.ascat_alleles ?: meta.ascat_alleles)] : null }.collect()
    ascat_loci = references.map { meta, _fasta -> meta.ascat_loci ? [meta.subMap(['id']), file(params.ascat_loci ?: meta.ascat_loci)] : null }.collect()
    ascat_loci_gc = references.map { meta, _fasta -> meta.ascat_loci_gc ? [meta.subMap(['id']), file(params.ascat_loci_gc ?: meta.ascat_loci_gc)] : null }.collect()
    ascat_loci_rt = references.map { meta, _fasta -> meta.ascat_loci_rt ? [meta.subMap(['id']), file(params.ascat_loci_rt ?: meta.ascat_loci_rt)] : null }.collect()
    bwamem1_index = references.map { meta, _fasta -> meta.bwamem1_index ? [meta.subMap(['id']), file(params.bwa ?: meta.bwamem1_index)] : null }.collect()
    bwamem2_index = references.map { meta, _fasta -> meta.bwamem2_index ? [meta.subMap(['id']), file(params.bwamem2 ?: meta.bwamem2_index)] : null }.collect()
    cf_chrom_len = references.map { meta, _fasta -> meta.cf_chrom_len ? [meta.subMap(['id']), file(params.cf_chrom_len ?: meta.cf_chrom_len)] : null }.collect()
    chr_dir = references.map { meta, _fasta -> meta.chr_dir ? [meta.subMap(['id']), file(params.chr_dir ?: meta.chr_dir)] : null }.collect()
    dbsnp = references.map { meta, _fasta -> meta.vcf.dbsnp.vcf ? [meta.subMap(['id']), file(meta.vcf.dbsnp.vcf)] : null }.transpose().collect()
    dbsnp_tbi = references.map { meta, _fasta -> meta.vcf.dbsnp.vcf_tbi ? [meta.subMap(['id']), file(meta.vcf.dbsnp.vcf_tbi)] : null }.transpose().collect()
    dragmap_hashtable = references.map { meta, _fasta -> meta.dragmap_hashtable ? [meta.subMap(['id']), file(meta.dragmap_hashtable)] : null }.collect()
    fasta = references.map { meta, fasta -> [meta.subMap(['id']), file(params.fasta ?: fasta)] }.collect()
    fasta_dict = references.map { meta, _fasta -> [meta.subMap(['id']), file(params.dict ?: meta.fasta_dict)] }.collect()
    fasta_fai = references.map { meta, _fasta -> [meta.subMap(['id']), file(params.fasta_fai ?: meta.fasta_fai)] }.collect()
    germline_resource = references.map { meta, _fasta -> meta.vcf.germline_resource.vcf ? [meta.subMap(['id']), file(meta.vcf.germline_resource.vcf)] : null }.collect()
    germline_resource_tbi = references.map { meta, _fasta -> meta.vcf.germline_resource.vcf_tbi ? [meta.subMap(['id']), file(meta.vcf.germline_resource.vcf_tbi)] : null }.collect()
    intervals_bed = references.map { meta, _fasta -> meta.intervals_bed ? [meta.subMap(['id']), file(params.intervals ?: meta.intervals_bed)] : null }.collect()
    known_indels = references.map { meta, _fasta -> meta.vcf.known_indels.vcf ? [meta.subMap(['id']), file(meta.vcf.known_indels.vcf)] : null }.collect()
    known_indels_tbi = references.map { meta, _fasta -> meta.vcf.known_indels.vcf_tbi ? [meta.subMap(['id']), file(meta.vcf.known_indels.vcf_tbi)] : null }.collect()
    known_snps = references.map { meta, _fasta -> meta.vcf.known_snps.vcf ? [meta.subMap(['id']), file(meta.vcf.known_snps.vcf)] : null }.collect()
    known_snps_tbi = references.map { meta, _fasta -> meta.vcf.known_snps.vcf_tbi ? [meta.subMap(['id']), file(meta.vcf.known_snps.vcf_tbi)] : null }.collect()
    mappability = []
    // mappability = references.map { meta, _fasta -> [meta.subMap(['id']), file(params.mappability ?: meta.mappability)] }.collect()
    msisensorpro_scan = []
    // msisensorpro_scan = references.map { meta, _fasta -> [meta.subMap(['id']), file(params.msisensorpro_scan ?: meta.msisensorpro_scan)] }.collect()
    ngscheckmate_bed = []
    // ngscheckmate_bed = references.map { meta, _fasta -> [meta.subMap(['id']), file(params.ngscheckmate_bed ?: meta.ngscheckmate_bed)] }.collect()
    pon = references.map { meta, _fasta -> meta.vcf.pon.vcf ? [meta.subMap(['id']), file(params.pon ?: meta.vcf.pon.vcf)] : null }.collect()
    pon_tbi = references.map { meta, _fasta -> meta.vcf.pon.vcf_tbi ? [meta.subMap(['id']), file(params.pon_tbi ?: meta.vcf.pon.vcf_tbi)] : null }.collect()
    sentieon_dnascope_model = references.map { meta, _fasta -> [meta.subMap(['id']), file(params.sentieon_dnascope_model ?: meta.sentieon_dnascope_model)] }.collect()

    // References' values from the references asset yaml or params
    ascat_genome = references.map { meta, _fasta -> meta.ascat_genome ? [meta.subMap(['id']), params.ascat_genome ?: meta.ascat_genome] : null }.collect()
    snpeff_db = references.map { meta, _fasta -> meta.snpeff_db ? [meta.subMap(['id']), params.snpeff_db ?: meta.snpeff_db] : null }.collect()
    vep_cache_version = references.map { meta, _fasta -> meta.vep_cache_version ? [meta.subMap(['id']), params.vep_cache_version ?: meta.vep_cache_version] : null }.collect()
    vep_genome = references.map { meta, _fasta -> meta.vep_genome ? [meta.subMap(['id']), params.vep_genome ?: meta.vep_genome] : null }.collect()
    vep_species = references.map { meta, _fasta -> meta.vep_species ? [meta.subMap(['id']), params.vep_species ?: meta.vep_species] : null }.collect()

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels = dbsnp.mix(known_indels).groupTuple().collect()
    known_sites_indels_tbi = dbsnp_tbi.mix(known_indels_tbi).groupTuple().collect()
    known_sites_snps = dbsnp.mix(known_snps).groupTuple().collect()
    known_sites_snps_tbi = dbsnp_tbi.mix(known_snps_tbi).groupTuple().collect()

    // Gather index for mapping given the chosen aligner
    index_alignment = aligner == "bwa-mem" || aligner == "sentieon-bwamem"
        ? bwamem1_index
        : aligner == "bwa-mem2"
            ? bwamem2_index
            : dragmap_hashtable

    bcftools_annotations = params.bcftools_annotations ? Channel.fromPath(params.bcftools_annotations).collect() : Channel.value([])
    bcftools_annotations_tbi = params.bcftools_annotations ? params.bcftools_annotations_tbi ? Channel.fromPath(params.bcftools_annotations_tbi).collect() : Channel.value([]) : Channel.value([])
    bcftools_header_lines = params.bcftools_header_lines ?: Channel.value([])

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
        ascat_alleles,
        bcftools_annotations,
        bcftools_annotations_tbi,
        bcftools_header_lines,
        cf_chrom_len,
        chr_dir,
        cnvkit_reference,
        dbsnp,
        dbsnp_tbi,
        "dbsnp_vqsr",
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
        "known_indels_vqsr",
        known_sites_indels,
        known_sites_indels_tbi,
        known_sites_snps,
        known_sites_snps_tbi,
        "known_snps_vqsr",
        ascat_loci,
        mappability,
        msisensorpro_scan,
        ngscheckmate_bed,
        pon,
        pon_tbi,
        ascat_loci_rt,
        sentieon_dnascope_model,
        snpeff_cache,
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
