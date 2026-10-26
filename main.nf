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

include { getGenomeAttributeStatic } from 'plugin/nf-core-utils'
include { getReferenceAttribute } from 'plugin/nf-core-utils'
include { referenceMeta } from 'plugin/nf-core-utils'
include { selectReferenceRow } from 'plugin/nf-core-utils'
include { samplesheetToList } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPORT MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

include { SAREK } from './workflows/sarek'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PREPARE_GENOME } from './subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS } from './subworkflows/local/prepare_intervals'
include { PREPARE_REFERENCE_CNVKIT } from './subworkflows/local/prepare_reference_cnvkit'
include { PREPARE_SNPSIFT_DATABASES } from './subworkflows/local/prepare_snpsift_databases'
include { CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF } from './subworkflows/nf-core/cache_download_ensemblvep_snpeff'
include { UTILS_ANNOTATION_CACHE } from './subworkflows/nf-core/utils_annotation_cache'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// WORKFLOW: Run main nf-core/sarek analysis pipeline
workflow NFCORE_SAREK {
    take:
    samplesheet
    datasheet

    main:
    // Resolve reference attributes with fallback to igenomes
    params.ascat_alleles = getReferenceAttribute(param: 'ascat_alleles', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.ascat_genome = getReferenceAttribute(param: 'ascat_genome', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.ascat_loci = getReferenceAttribute(param: 'ascat_loci', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.ascat_loci_gc = getReferenceAttribute(param: 'ascat_loci_gc', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.ascat_loci_rt = getReferenceAttribute(param: 'ascat_loci_rt', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.bwa = getReferenceAttribute(param: 'bwa', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.bwamem2 = getReferenceAttribute(param: 'bwamem2', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.cf_chrom_len = getReferenceAttribute(param: 'cf_chrom_len', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.chr_dir = getReferenceAttribute(param: 'chr_dir', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.dbsnp = getReferenceAttribute(param: 'dbsnp', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.dbsnp_tbi = getReferenceAttribute(param: 'dbsnp_tbi', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.dbsnp_vqsr = getReferenceAttribute(param: 'dbsnp_vqsr', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.dict = getReferenceAttribute(param: 'dict', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.dragmap = getReferenceAttribute(param: 'dragmap', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.fasta = getReferenceAttribute(param: 'fasta', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.fasta_fai = getReferenceAttribute(param: 'fasta_fai', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.germline_resource = getReferenceAttribute(param: 'germline_resource', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.germline_resource_tbi = getReferenceAttribute(param: 'germline_resource_tbi', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.intervals = getReferenceAttribute(param: 'intervals', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.known_indels = getReferenceAttribute(param: 'known_indels', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.known_indels_tbi = getReferenceAttribute(param: 'known_indels_tbi', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.known_indels_vqsr = getReferenceAttribute(param: 'known_indels_vqsr', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.known_snps = getReferenceAttribute(param: 'known_snps', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.known_snps_tbi = getReferenceAttribute(param: 'known_snps_tbi', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.known_snps_vqsr = getReferenceAttribute(param: 'known_snps_vqsr', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.mappability = getReferenceAttribute(param: 'mappability', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.msisensor2_models = getReferenceAttribute(param: 'msisensor2_models', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.msisensorpro_scan = getReferenceAttribute(param: 'msisensorpro_scan', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.ngscheckmate_bed = getReferenceAttribute(param: 'ngscheckmate_bed', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.pon = getReferenceAttribute(param: 'pon', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.pon_tbi = getReferenceAttribute(param: 'pon_tbi', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.sentieon_dnascope_model = getReferenceAttribute(param: 'sentieon_dnascope_model', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.snpeff_db = getReferenceAttribute(param: 'snpeff_db', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.vep_cache_version = getReferenceAttribute(param: 'vep_cache_version', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.vep_genome = getReferenceAttribute(param: 'vep_genome', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)
    params.vep_species = getReferenceAttribute(param: 'vep_species', datasheet: datasheet, genomeMap: params.genomes, genomeKey: params.genome)

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
    intervals_bed_combined = params.no_intervals ? channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined = params.no_intervals ? channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined
    intervals_bed_combined_for_variant_calling = PREPARE_INTERVALS.out.intervals_bed_combined

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes
        ? intervals_bed_combined.map { bed -> [[id: bed.baseName], bed] }.collect()
        : channel.value([[id: 'null'], []])
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
            cnvkit_reference = channel.fromPath(params.cnvkit_reference).collect()
        }
        else {
            PREPARE_REFERENCE_CNVKIT(PREPARE_GENOME.out.fasta, intervals_bed_combined)
            cnvkit_reference = PREPARE_REFERENCE_CNVKIT.out.cnvkit_reference
        }
    }
    else {
        cnvkit_reference = channel.value([])
    }

    // Fails when consensus calling is specified without normalization
    if (params.snv_consensus_calling && !params.normalize_vcfs) {
        error("Consensus calling was specified without normalization. Set --normalize_vcfs in addition. See: https://www.biostars.org/p/307035/")
    }


    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, sarek will download the cache
        ensemblvep_info = channel.of([[id: "${params.vep_cache_version}_${params.vep_genome}"], params.vep_genome, params.vep_species, params.vep_cache_version])
        snpeff_info = channel.of([[id: "${params.snpeff_db}"], params.snpeff_db])
        CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF(ensemblvep_info, snpeff_info, params.vep_cache_preflight_check)
        snpeff_cache = CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF.out.snpeff_cache
        vep_cache = CACHE_DOWNLOAD_ENSEMBLVEP_SNPEFF.out.ensemblvep_cache
    }
    else {
        // Looks for cache information either locally or on the cloud
        UTILS_ANNOTATION_CACHE(
            params.vep_cache,
            params.vep_cache_version,
            params.vep_custom_args,
            params.vep_genome,
            params.vep_species,
            (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            params.snpeff_db,
            (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            "Please refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.",
        )

        snpeff_cache = UTILS_ANNOTATION_CACHE.out.snpeff_cache
        vep_cache = UTILS_ANNOTATION_CACHE.out.ensemblvep_cache
    }

    vep_extra_files = []

    if (params.dbnsfp && params.dbnsfp_tbi) {
        vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
        vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
    }
    else if (params.dbnsfp && !params.dbnsfp_tbi) {
        System.err.println("DBNSFP: ${params.dbnsfp} has been provided with `--dbnsfp, but no dbnsfp_tbi has")
        System.err.println("cf: https://nf-co.re/sarek/parameters#dbnsfp")
        error("Execution halted due to dbnsfp inconsistency.")
    }

    if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
        vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
    }

    if (params.phenotypes_file) {
        vep_extra_files.add(file(params.phenotypes_file, checkIfExists: true))
        if (params.phenotypes_file_tbi) {
            vep_extra_files.add(file(params.phenotypes_file_tbi, checkIfExists: true))
        }
    }

    // Build SnpSift annotation databases configuration from CSV samplesheet
    // CSV format: vcf,tbi,fields,prefix,vardb
    // - vcf: Path to annotation VCF (required)
    // - tbi: Path to tabix index (optional, defaults to ${vcf}.tbi)
    // - fields: Semicolon-separated INFO fields to extract (optional)
    // - prefix: Prefix for annotated field names (optional)
    // - vardb: Path to pre-built .snpsift.vardb directory (optional)
    snpsift_db_configs = []

    if (params.snpsift_databases) {
        // Parse and validate CSV using nf-schema
        // Returns list of tuples: [vcf, tbi, fields, prefix, vardb]
        def db_list = samplesheetToList(params.snpsift_databases, "${projectDir}/assets/schema_snpsift_databases.json")

        db_list.each { vcf, tbi, fields, prefix, vardb ->
            // Fields are required when vardb is not provided (needed to build the database)
            if (!vardb && !fields) {
                error("SnpSift database '${vcf}': 'fields' column is required when 'vardb' is not provided (needed for database creation)")
            }

            def vcf_file = file(vcf, checkIfExists: true)
            def tbi_file = tbi ? file(tbi, checkIfExists: true) : file("${vcf}.tbi", checkIfExists: true)
            def vardb_file = vardb ? file(vardb, checkIfExists: true) : null

            snpsift_db_configs.add(
                [
                    vcf: vcf_file,
                    tbi: tbi_file,
                    fields: fields ?: '',
                    prefix: prefix ?: '',
                    vardb: vardb_file,
                ]
            )
        }
    }

    // Prepare SnpSift databases (build if vardb not provided, returns tuple for SNPSIFT_ANNMEM)
    ch_snpsift_db = channel.value([[], [], [], [], []])
    if (params.tools && params.tools.split(',').contains('snpsift') && snpsift_db_configs) {
        PREPARE_SNPSIFT_DATABASES(snpsift_db_configs)
        ch_snpsift_db = PREPARE_SNPSIFT_DATABASES.out.db_tuple
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
        params.bcftools_columns ? channel.fromPath(params.bcftools_columns).collect() : channel.value([[]]),
        params.bcftools_header_lines ? channel.fromPath(params.bcftools_header_lines).collect() : channel.empty(),
        params.cf_chrom_len ? channel.fromPath(params.cf_chrom_len).collect() : [],
        PREPARE_GENOME.out.chr_dir,
        cnvkit_reference,
        PREPARE_GENOME.out.dbsnp,
        PREPARE_GENOME.out.dbsnp_tbi,
        params.dbsnp_vqsr ? channel.value(params.dbsnp_vqsr) : channel.empty(),
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
        params.known_indels_vqsr ? channel.value(params.known_indels_vqsr) : channel.empty(),
        PREPARE_GENOME.out.known_sites_indels,
        PREPARE_GENOME.out.known_sites_indels_tbi,
        PREPARE_GENOME.out.known_sites_snps,
        PREPARE_GENOME.out.known_sites_snps_tbi,
        params.known_snps_vqsr ? channel.value(params.known_snps_vqsr) : channel.empty(),
        params.mappability ? channel.fromPath(params.mappability).collect() : channel.value([]),
        PREPARE_GENOME.out.msisensor2_models,
        PREPARE_GENOME.out.msisensorpro_scan,
        params.ngscheckmate_bed ? channel.value(params.ngscheckmate_bed) : channel.empty(),
        PREPARE_GENOME.out.pon,
        PREPARE_GENOME.out.pon_tbi,
        params.sentieon_dnascope_model ? channel.fromPath(params.sentieon_dnascope_model).collect() : channel.value([]),
        params.varlociraptor_chunk_size,
        params.varlociraptor_events_germline,
        params.varlociraptor_events_somatic,
        params.varlociraptor_events_tumor_only,
        params.varlociraptor_fdr,
        params.varlociraptor_scenario_germline ? channel.fromPath(params.varlociraptor_scenario_germline).collect() : channel.fromPath("${projectDir}/assets/varlociraptor_germline.yte.yaml").collect(),
        params.varlociraptor_scenario_somatic ? channel.fromPath(params.varlociraptor_scenario_somatic).collect() : channel.fromPath("${projectDir}/assets/varlociraptor_somatic.yte.yaml").collect(),
        params.varlociraptor_scenario_tumor_only ? channel.fromPath(params.varlociraptor_scenario_tumor_only).collect() : channel.fromPath("${projectDir}/assets/varlociraptor_tumor_only.yte.yaml").collect(),
        snpeff_cache,
        params.snpeff_db,
        vep_cache,
        params.vep_cache_version,
        vep_extra_files,
        PREPARE_GENOME.out.vep_fasta,
        params.vep_genome,
        params.vep_species,
        ch_snpsift_db,
    )

    emit:
    multiqc_publish = SAREK.out.multiqc_publish
    multiqc_report = SAREK.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SAREK(PIPELINE_INITIALISATION.out.samplesheet, PIPELINE_INITIALISATION.out.datasheet.getVal())

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_SAREK.out.multiqc_report,
    )

    publish:
    multiqc = NFCORE_SAREK.out.multiqc_publish
}

output {
    multiqc {
        path "multiqc"
        index {
            path "multiqc/index.json"
            sep ":"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
