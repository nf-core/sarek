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

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK  } from './workflows/sarek'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_sarek_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_sarek_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.ascat_alleles           = WorkflowMain.getGenomeAttribute(params, 'ascat_alleles')
params.ascat_genome            = WorkflowMain.getGenomeAttribute(params, 'ascat_genome')
params.ascat_loci              = WorkflowMain.getGenomeAttribute(params, 'ascat_loci')
params.ascat_loci_gc           = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_gc')
params.ascat_loci_rt           = WorkflowMain.getGenomeAttribute(params, 'ascat_loci_rt')
params.bwa                     = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2                 = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.cf_chrom_len            = WorkflowMain.getGenomeAttribute(params, 'cf_chrom_len')
params.chr_dir                 = WorkflowMain.getGenomeAttribute(params, 'chr_dir')
params.dbsnp                   = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi               = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.dbsnp_vqsr              = WorkflowMain.getGenomeAttribute(params, 'dbsnp_vqsr')
params.dict                    = WorkflowMain.getGenomeAttribute(params, 'dict')
params.dragmap                 = WorkflowMain.getGenomeAttribute(params, 'dragmap')
params.fasta                   = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai               = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.germline_resource       = WorkflowMain.getGenomeAttribute(params, 'germline_resource')
params.germline_resource_tbi   = WorkflowMain.getGenomeAttribute(params, 'germline_resource_tbi')
params.intervals               = WorkflowMain.getGenomeAttribute(params, 'intervals')
params.known_indels            = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi        = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')
params.known_indels_vqsr       = WorkflowMain.getGenomeAttribute(params, 'known_indels_vqsr')
params.known_snps              = WorkflowMain.getGenomeAttribute(params, 'known_snps')
params.known_snps_tbi          = WorkflowMain.getGenomeAttribute(params, 'known_snps_tbi')
params.known_snps_vqsr         = WorkflowMain.getGenomeAttribute(params, 'known_snps_vqsr')
params.mappability             = WorkflowMain.getGenomeAttribute(params, 'mappability')
params.minimap2                = WorkflowMain.getGenomeAttribute(params, 'minimap2')
params.ngscheckmate_bed        = WorkflowMain.getGenomeAttribute(params, 'ngscheckmate_bed')
params.pon                     = WorkflowMain.getGenomeAttribute(params, 'pon')
params.pon_tbi                 = WorkflowMain.getGenomeAttribute(params, 'pon_tbi')
params.sentieon_dnascope_model = WorkflowMain.getGenomeAttribute(params, 'sentieon_dnascope_model')
params.snpeff_db               = WorkflowMain.getGenomeAttribute(params, 'snpeff_db')
params.snpeff_genome           = WorkflowMain.getGenomeAttribute(params, 'snpeff_genome')
params.vep_cache_version       = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.vep_genome              = WorkflowMain.getGenomeAttribute(params, 'vep_genome')
params.vep_species             = WorkflowMain.getGenomeAttribute(params, 'vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK                            } from './workflows/sarek'
include { ANNOTATION_CACHE_INITIALISATION  } from './subworkflows/local/annotation_cache_initialisation'
include { DOWNLOAD_CACHE_SNPEFF_VEP        } from './subworkflows/local/download_cache_snpeff_vep'
include { PIPELINE_COMPLETION              } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PIPELINE_INITIALISATION          } from './subworkflows/local/utils_nfcore_sarek_pipeline'
include { PREPARE_GENOME                   } from './subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS                } from './subworkflows/local/prepare_intervals'
include { PREPARE_REFERENCE_CNVKIT         } from './subworkflows/local/prepare_reference_cnvkit'

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
bcftools_annotations    = params.bcftools_annotations    ? Channel.fromPath(params.bcftools_annotations).collect()  : Channel.empty()
bcftools_header_lines   = params.bcftools_header_lines   ? Channel.fromPath(params.bcftools_header_lines).collect() : Channel.empty()
cf_chrom_len            = params.cf_chrom_len            ? Channel.fromPath(params.cf_chrom_len)                    : []
dbsnp                   = params.dbsnp                   ? Channel.fromPath(params.dbsnp).collect()                 : Channel.value([])
fasta                   = params.fasta                   ? Channel.fromPath(params.fasta).collect()                 : Channel.empty()
fasta_fai               = params.fasta_fai               ? Channel.fromPath(params.fasta_fai).collect()             : Channel.empty()
germline_resource       = params.germline_resource       ? Channel.fromPath(params.germline_resource)               : Channel.value([]) // Mutect2 does not require a germline resource, so set to optional input
known_indels            = params.known_indels            ? Channel.fromPath(params.known_indels).collect()          : Channel.value([])
known_snps              = params.known_snps              ? Channel.fromPath(params.known_snps).collect()            : Channel.value([])
mappability             = params.mappability             ? Channel.fromPath(params.mappability)                     : Channel.value([])
pon                     = params.pon                     ? Channel.fromPath(params.pon)                             : Channel.value([]) // PON is optional for Mutect2 (but highly recommended)
sentieon_dnascope_model = params.sentieon_dnascope_model ? Channel.fromPath(params.sentieon_dnascope_model)         : Channel.value([])

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
ascat_genome                = params.ascat_genome             ?: Channel.empty()
dbsnp_vqsr                  = params.dbsnp_vqsr               ? Channel.value(params.dbsnp_vqsr) : Channel.empty()
known_indels_vqsr           = params.known_indels_vqsr        ? Channel.value(params.known_indels_vqsr) : Channel.empty()
known_snps_vqsr             = params.known_snps_vqsr          ? Channel.value(params.known_snps_vqsr) : Channel.empty()
ngscheckmate_bed            = params.ngscheckmate_bed         ? Channel.value(params.ngscheckmate_bed) : Channel.empty()
snpeff_db                   = params.snpeff_db                ?: Channel.empty()
vep_cache_version           = params.vep_cache_version        ?: Channel.empty()
vep_genome                  = params.vep_genome               ?: Channel.empty()
vep_species                 = params.vep_species              ?: Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GATK.GRCh38 -profile docker --outdir results"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAREK } from './workflows/sarek'

// WORKFLOW: Run main nf-core/sarek analysis pipeline
workflow NFCORE_SAREK {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    versions = Channel.empty()

    // build indexes if needed
    PREPARE_GENOME(
        params.ascat_alleles,
        params.ascat_loci,
        params.ascat_loci_gc,
        params.ascat_loci_rt,
        bcftools_annotations,
        params.chr_dir,
        dbsnp,
        fasta,
        fasta_fai,
        germline_resource,
        known_indels,
        known_snps,
        pon)

    // Gather built indices or get them from the params
    // Built from the fasta file:
    dict       = params.dict        ? Channel.fromPath(params.dict).map{ it -> [ [id:'dict'], it ] }
                                    : PREPARE_GENOME.out.dict
    fasta_fai  = params.fasta_fai   ? Channel.fromPath(params.fasta_fai)
                                    : PREPARE_GENOME.out.fasta_fai
    bwa        = params.bwa         ? Channel.fromPath(params.bwa).collect()
                                    : PREPARE_GENOME.out.bwa
    bwamem2    = params.bwamem2     ? Channel.fromPath(params.bwamem2).collect()
                                    : PREPARE_GENOME.out.bwamem2
    dragmap    = params.dragmap     ? Channel.fromPath(params.dragmap).collect()
                                    : PREPARE_GENOME.out.hashtable
    minimap2   = params.dragmap     ? Channel.fromPath(params.minimap2).collect()
                                    : PREPARE_GENOME.out.hashtable

    // Gather index for mapping given the chosen aligner
    index_alignement = (params.aligner == "bwa-mem" || params.aligner == "sentieon-bwamem") ? bwa :
        params.aligner == "bwa-mem2" ? bwamem2 : params.aligner == "dragmap" ? dragmap : minimap2

    // TODO: add a params for msisensorpro_scan
    msisensorpro_scan      = PREPARE_GENOME.out.msisensorpro_scan

    // For ASCAT, extracted from zip or tar.gz files
    allele_files           = PREPARE_GENOME.out.allele_files
    chr_files              = PREPARE_GENOME.out.chr_files
    gc_file                = PREPARE_GENOME.out.gc_file
    loci_files             = PREPARE_GENOME.out.loci_files
    rt_file                = PREPARE_GENOME.out.rt_file

    // Tabix indexed vcf files
    bcftools_annotations_tbi  = params.bcftools_annotations    ? params.bcftools_annotations_tbi ? Channel.fromPath(params.bcftools_annotations_tbi)   : PREPARE_GENOME.out.bcftools_annotations_tbi : Channel.empty([])
    dbsnp_tbi                 = params.dbsnp                   ? params.dbsnp_tbi                ? Channel.fromPath(params.dbsnp_tbi)                  : PREPARE_GENOME.out.dbsnp_tbi                : Channel.value([])
    germline_resource_tbi     = params.germline_resource       ? params.germline_resource_tbi    ? Channel.fromPath(params.germline_resource_tbi)      : PREPARE_GENOME.out.germline_resource_tbi    : [] //do not change to Channel.value([]), the check for its existence then fails for Getpileupsumamries
    known_indels_tbi          = params.known_indels            ? params.known_indels_tbi         ? Channel.fromPath(params.known_indels_tbi).collect() : PREPARE_GENOME.out.known_indels_tbi         : Channel.value([])
    known_snps_tbi            = params.known_snps              ? params.known_snps_tbi           ? Channel.fromPath(params.known_snps_tbi)             : PREPARE_GENOME.out.known_snps_tbi           : Channel.value([])
    pon_tbi                   = params.pon                     ? params.pon_tbi                  ? Channel.fromPath(params.pon_tbi)                    : PREPARE_GENOME.out.pon_tbi                  : Channel.value([])

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()
    known_sites_snps       = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi   = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai, params.intervals, params.no_intervals, params.nucleotides_per_second, params.outdir, params.step)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined         = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined  = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined
    intervals_bed_combined_for_variant_calling = PREPARE_INTERVALS.out.intervals_bed_combined

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes ?
        intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect() :
        Channel.value([ [ id:'null' ], [] ])
    intervals            = PREPARE_INTERVALS.out.intervals_bed        // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather
    intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }
    intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map{ intervals, num_intervals ->
        if ( num_intervals < 1 ) [ [], [], num_intervals ]
        else [ intervals[0], intervals[1], num_intervals ]
    }
    if (params.tools && params.tools.split(',').contains('cnvkit')) {
        if (params.cnvkit_reference) {
            cnvkit_reference = Channel.fromPath(params.cnvkit_reference).collect()
        } else {
            PREPARE_REFERENCE_CNVKIT(fasta, intervals_bed_combined)
            cnvkit_reference = PREPARE_REFERENCE_CNVKIT.out.cnvkit_reference
            versions = versions.mix(PREPARE_REFERENCE_CNVKIT.out.versions)
        }
    } else {
        cnvkit_reference = Channel.value([])
    }
    // Gather used softwares versions
    versions = versions.mix(PREPARE_GENOME.out.versions)
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    vep_fasta = (params.vep_include_fasta) ? fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] } : [[id: 'null'], []]

    // Download cache
    if (params.download_cache) {
        // Assuming that even if the cache is provided, if the user specify download_cache, sarek will download the cache
        ensemblvep_info = Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
        snpeff_info     = Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])
        DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)
        snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
        vep_cache    = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        versions = versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    } else {
        // Looks for cache information either locally or on the cloud
        ANNOTATION_CACHE_INITIALISATION(
            (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
            params.snpeff_cache,
            params.snpeff_genome,
            params.snpeff_db,
            (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
            params.vep_cache,
            params.vep_species,
            params.vep_cache_version,
            params.vep_genome,
            "Please refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")

            snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
            vep_cache    = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    }

    //
    // WORKFLOW: Run pipeline
    //
    SAREK (
        samplesheet
    )

    emit:
    multiqc_report = SAREK.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_SAREK (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_SAREK.out.multiqc_report
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
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
