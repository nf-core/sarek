#!/usr/bin/env nextflow

/*
--------------------------------------------------------------------------------
                                  nf-core/sarek
--------------------------------------------------------------------------------
Started March 2016.
Ported to nf-core May 2019.
Ported to DSL 2 July 2020.
--------------------------------------------------------------------------------
nf-core/sarek:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://nf-co.re/sarek
--------------------------------------------------------------------------------
 @Documentation
 https://nf-co.re/sarek/latest/usage
--------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/sarek -profile docker --input sample.tsv --genome GRCh37"
    log.info Schema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --          INCLUDE SAREK FUNCTIONS         -- */
////////////////////////////////////////////////////

include {
    check_parameter_existence;
    check_parameter_list;
    define_anno_list;
    define_skip_qc_list;
    define_step_list;
    define_tool_list;
    extract_bam;
    extract_fastq;
    extract_fastq_from_dir;
    extract_recal;
    has_extension
} from './modules/local/functions'

////////////////////////////////////////////////////
/* --      SET UP CONFIGURATION VARIABLES      -- */
////////////////////////////////////////////////////

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

// Check that conda channels are set-up correctly
if (params.enable_conda) {
    Checks.check_conda_channels(log)
}

// Check AWS batch settings
Checks.aws_batch(workflow, params)

// Check the hostnames against configured profiles
Checks.hostname(workflow, params, log)

// Check genome key exists if provided
Checks.genome_exists(params, log)

step_list = define_step_list()
step = params.step ? params.step.toLowerCase().replaceAll('-', '').replaceAll('_', '') : ''

if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!check_parameter_existence(step, step_list)) exit 1, "Unknown step ${step}, see --help for more information"

tool_list = define_tool_list()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
if (step == 'controlfreec') tools = ['controlfreec']
if (!check_parameter_list(tools, tool_list)) exit 1, 'Unknown tool(s), see --help for more information'

skip_qc_list = define_skip_qc_list()
skip_qc = params.skip_qc ? params.skip_qc == 'all' ? skip_qc_list : params.skip_qc.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
if (!check_parameter_list(skip_qc, skip_qc_list)) exit 1, 'Unknown QC tool(s), see --help for more information'

anno_list = define_anno_list()
annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '')} : []
if (!check_parameter_list(annotate_tools,anno_list)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'

if (!(params.aligner in ['bwa-mem', 'bwa-mem2'])) exit 1, 'Unknown aligner, see --help for more information'

// // Check parameters
if ((params.ascat_ploidy && !params.ascat_purity) || (!params.ascat_ploidy && params.ascat_purity)) exit 1, 'Please specify both --ascat_purity and --ascat_ploidy, or none of them'
if (params.cf_window && params.cf_coeff) exit 1, 'Please specify either --cf_window OR --cf_coeff, but not both of them'
if (params.umi && !(params.read_structure1 && params.read_structure2)) exit 1, 'Please specify both --read_structure1 and --read_structure2, when using --umi'

// Handle input
tsv_path = null
if (params.input && (has_extension(params.input, "tsv") || has_extension(params.input, "vcf") || has_extension(params.input, "vcf.gz"))) tsv_path = params.input
if (params.input && (has_extension(params.input, "vcf") || has_extension(params.input, "vcf.gz"))) step = "annotate"

save_bam_mapped = params.skip_markduplicates ? true : params.save_bam_mapped ? true : false

// If no input file specified, trying to get TSV files corresponding to step in the TSV directory
// only for steps preparerecalibration, recalibrate, variantcalling and controlfreec
if (!params.input && params.sentieon) {
    switch (step) {
        case 'mapping': break
        case 'recalibrate': tsv_path = "${params.outdir}/preprocessing/tsv/sentieon_deduped.tsv"; break
        case 'variantcalling': tsv_path = "${params.outdir}/preprocessing/tsv/sentieon_recalibrated.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && !params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsv_path = "${params.outdir}/preprocessing/tsv/markduplicates_no_table.tsv"; break
        case 'recalibrate': tsv_path = "${params.outdir}/preprocessing/tsv/markduplicates.tsv"; break
        case 'variantcalling': tsv_path = "${params.outdir}/preprocessing/tsv/recalibrated.tsv"; break
        case 'controlfreec': tsv_path = "${params.outdir}/variant_calling/tsv/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsv_path = "${params.outdir}/preprocessing/tsv/mapped.tsv"; break
        case 'recalibrate': tsv_path = "${params.outdir}/preprocessing/tsv/mapped_no_markduplicates.tsv"; break
        case 'variantcalling': tsv_path = "${params.outdir}/preprocessing/tsv/recalibrated.tsv"; break
        case 'controlfreec': tsv_path = "${params.outdir}/variant_calling/tsv/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
}

input_sample = Channel.empty()
if (tsv_path) {
    tsv_file = file(tsv_path)
    switch (step) {
        case 'mapping': input_sample = extract_fastq(tsv_file); break
        case 'preparerecalibration': input_sample = extract_bam(tsv_file); break
        case 'recalibrate': input_sample = extract_recal(tsv_file); break
        case 'variantcalling': input_sample = extract_bam(tsv_file); break
        case 'controlfreec': input_sample = extract_pileup(tsv_file); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (params.input && !has_extension(params.input, "tsv")) {
    log.info "No TSV file"
    if (step != 'mapping') exit 1, 'No step other than "mapping" supports a directory as an input'
    log.info "Reading ${params.input} directory"
    log.warn "[nf-core/sarek] in ${params.input} directory, all fastqs are assuming to be from the same sample, which is assumed to be a germline one"
    input_sample = extract_fastq_from_dir(params.input)
    tsv_file = params.input  // used in the reports
} else if (tsv_path && step == 'annotate') {
    log.info "Annotating ${tsv_path}"
} else if (step == 'annotate') {
    log.info "Trying automatic annotation on files in the VariantCalling/ directory"
} else exit 1, 'No sample were defined, see --help'

////////////////////////////////////////////////////
/* --  UPDATE MODULES OPTIONS BASED ON PARAMS  -- */
////////////////////////////////////////////////////

modules = params.modules

if (params.save_reference)      modules['build_intervals'].publish_files         = ['bed':'intervals']
if (params.save_reference)      modules['bwa_index'].publish_files               = ['amb':'bwa', 'ann':'bwa', 'bwt':'bwa', 'pac':'bwa', 'sa':'bwa']
if (params.save_reference)      modules['bwamem2_index'].publish_files           = ['0123':'bwamem2', 'amb':'bwamem2', 'ann':'bwamem2', 'bwt.2bit.64':'bwamem2', 'bwt.8bit.32':'bwamem2', 'pac':'bwamem2']
if (params.save_reference)      modules['create_intervals_bed'].publish_files    = ['bed':'intervals']
if (params.save_reference)      modules['dict'].publish_files                    = ['dict':'dict']
if (params.save_reference)      modules['index_target_bed'].publish_files        = ['bed.gz':'target', 'bed.gz.tbi':'target']
if (params.save_reference)      modules['msisensor_scan'].publish_files          = ['list':'msi']
if (params.save_reference)      modules['samtools_faidx'].publish_files          = ['fai':'fai']
if (params.save_reference)      modules['tabix_dbsnp'].publish_files             = ['vcf.gz.tbi':'dbsnp']
if (params.save_reference)      modules['tabix_germline_resource'].publish_files = ['vcf.gz.tbi':'germline_resource']
if (params.save_reference)      modules['tabix_known_indels'].publish_files      = ['vcf.gz.tbi':'known_indels']
if (params.save_reference)      modules['tabix_pon'].publish_files               = ['vcf.gz.tbi':'pon']
if (save_bam_mapped)            modules['samtools_index_mapping'].publish_files  = ['bam':'mapped', 'bai':'mapped']
if (params.skip_markduplicates) modules['baserecalibrator'].publish_files        = ['recal.table':'mapped']
if (params.skip_markduplicates) modules['gatherbqsrreports'].publish_files       = ['recal.table':'mapped']

////////////////////////////////////////////////////
/* --        REFERENCES PARAMETER VALUES       -- */
////////////////////////////////////////////////////
/* -- Initialize each params in params.genomes -- */
/* --  catch the command line first if defined -- */
////////////////////////////////////////////////////

params.ac_loci                 = Checks.get_genome_attribute(params, 'ac_loci')
params.ac_loci_gc              = Checks.get_genome_attribute(params, 'ac_loci_gc')
params.bwa                     = Checks.get_genome_attribute(params, 'bwa')
params.chr_dir                 = Checks.get_genome_attribute(params, 'chr_dir')
params.chr_length              = Checks.get_genome_attribute(params, 'chr_length')
params.dbsnp                   = Checks.get_genome_attribute(params, 'dbsnp')
params.dbsnp_index             = Checks.get_genome_attribute(params, 'dbsnp_index')
params.dict                    = Checks.get_genome_attribute(params, 'dict')
params.fasta                   = Checks.get_genome_attribute(params, 'fasta')
params.fasta_fai               = Checks.get_genome_attribute(params, 'fasta_fai')
params.germline_resource       = Checks.get_genome_attribute(params, 'germline_resource')
params.germline_resource_index = Checks.get_genome_attribute(params, 'germline_resource_index')
params.intervals               = Checks.get_genome_attribute(params, 'intervals')
params.known_indels            = Checks.get_genome_attribute(params, 'known_indels')
params.known_indels_index      = Checks.get_genome_attribute(params, 'known_indels_index')
params.mappability             = Checks.get_genome_attribute(params, 'mappability')
params.snpeff_db               = Checks.get_genome_attribute(params, 'snpeff_db')
params.species                 = Checks.get_genome_attribute(params, 'species')
params.vep_cache_version       = Checks.get_genome_attribute(params, 'vep_cache_version')

file("${params.outdir}/no_file").text = "no_file\n"

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
chr_dir           = params.chr_dir           ? file(params.chr_dir)           : file("${params.outdir}/no_file")
chr_length        = params.chr_length        ? file(params.chr_length)        : file("${params.outdir}/no_file")
dbsnp             = params.dbsnp             ? file(params.dbsnp)             : file("${params.outdir}/no_file")
fasta             = params.fasta             ? file(params.fasta)             : file("${params.outdir}/no_file")
germline_resource = params.germline_resource ? file(params.germline_resource) : file("${params.outdir}/no_file")
known_indels      = params.known_indels      ? file(params.known_indels)      : file("${params.outdir}/no_file")
loci              = params.ac_loci           ? file(params.ac_loci)           : file("${params.outdir}/no_file")
loci_gc           = params.ac_loci_gc        ? file(params.ac_loci_gc)        : file("${params.outdir}/no_file")
mappability       = params.mappability       ? file(params.mappability)       : file("${params.outdir}/no_file")

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
snpeff_db         = params.snpeff_db         ?: Channel.empty()
snpeff_species    = params.species           ?: Channel.empty()
vep_cache_version = params.vep_cache_version ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
cadd_indels       = params.cadd_indels       ? file(params.cadd_indels)      : file("${params.outdir}/no_file")
cadd_indels_tbi   = params.cadd_indels_tbi   ? file(params.cadd_indels_tbi)  : file("${params.outdir}/no_file")
cadd_wg_snvs      = params.cadd_wg_snvs      ? file(params.cadd_wg_snvs)     : file("${params.outdir}/no_file")
cadd_wg_snvs_tbi  = params.cadd_wg_snvs_tbi  ? file(params.cadd_wg_snvs_tbi) : file("${params.outdir}/no_file")
pon               = params.pon               ? file(params.pon)              : file("${params.outdir}/no_file")
snpeff_cache      = params.snpeff_cache      ? file(params.snpeff_cache)     : file("${params.outdir}/no_file")
target_bed        = params.target_bed        ? file(params.target_bed)       : file("${params.outdir}/no_file")
vep_cache         = params.vep_cache         ? file(params.vep_cache)        : file("${params.outdir}/no_file")

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
read_structure1   = params.read_structure1   ?: Channel.empty()
read_structure2   = params.read_structure2   ?: Channel.empty()

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

if ('mutect2' in tools && !(params.pon)) log.warn "[nf-core/sarek] Mutect2 was requested, but as no panel of normals were given, results will not be optimal"
if (params.sentieon) log.warn "[nf-core/sarek] Sentieon will be used, only works if Sentieon is available where nf-core/sarek is run"

////////////////////////////////////////////////////
/* --        INCLUDE LOCAL SUBWORKFLOWS        -- */
////////////////////////////////////////////////////

include { BUILD_INDICES } from './modules/local/subworkflow/build_indices' addParams(
    build_intervals_options:         modules['build_intervals'],
    bwa_index_options:               modules['bwa_index'],
    bwamem2_index_options:           modules['bwamem2_index'],
    create_intervals_bed_options:    modules['create_intervals_bed'],
    gatk_dict_options:               modules['dict'],
    index_target_bed_options:        modules['index_target_bed'],
    msisensor_scan_options:          modules['msisensor_scan'],
    samtools_faidx_options:          modules['samtools_faidx'],
    tabix_dbsnp_options:             modules['tabix_dbsnp'],
    tabix_germline_resource_options: modules['tabix_germline_resource'],
    tabix_known_indels_options:      modules['tabix_known_indels'],
    tabix_pon_options:               modules['tabix_pon']
)
include { MAPPING } from './modules/local/subworkflow/mapping' addParams(
    bwamem1_mem_options:             modules['bwa_mem1_mem'],
    bwamem2_mem_options:             modules['bwa_mem2_mem'],
    merge_bam_options:               modules['merge_bam_mapping'],
    qualimap_bamqc_options:          modules['qualimap_bamqc_mapping'],
    samtools_index_options:          modules['samtools_index_mapping'],
    samtools_stats_options:          modules['samtools_stats_mapping']
)
include { MARKDUPLICATES } from './modules/local/subworkflow/markduplicates' addParams(
    markduplicates_options:          modules['markduplicates']
)
include { PREPARE_RECALIBRATION } from './modules/local/subworkflow/prepare_recalibration' addParams(
    baserecalibrator_options:        modules['baserecalibrator'],
    gatherbqsrreports_options:       modules['gatherbqsrreports']
)
include { RECALIBRATE } from './modules/local/subworkflow/recalibrate' addParams(
    applybqsr_options:               modules['applybqsr'],
    merge_bam_options:               modules['merge_bam_recalibrate'],
    qualimap_bamqc_options:          modules['qualimap_bamqc_recalibrate'],
    samtools_index_options:          modules['samtools_index_recalibrate'],
    samtools_stats_options:          modules['samtools_stats_recalibrate']
)
include { GERMLINE_VARIANT_CALLING } from './modules/local/subworkflow/germline_variant_calling' addParams(
    concat_gvcf_options:             modules['concat_gvcf'],
    concat_haplotypecaller_options:  modules['concat_haplotypecaller'],
    genotypegvcf_options:            modules['genotypegvcf'],
    haplotypecaller_options:         modules['haplotypecaller'],
    strelka_options:                 modules['strelka_germline']
)
// include { TUMOR_VARIANT_CALLING } from './modules/local/subworkflow/tumor_variant_calling' addParams(
// )
include { PAIR_VARIANT_CALLING } from './modules/local/subworkflow/pair_variant_calling' addParams(
    manta_options:                   modules['manta_somatic'],
    msisensor_msi_options:           modules['msisensor_msi'],
    strelka_bp_options:              modules['strelka_somatic_bp'],
    strelka_options:                 modules['strelka_somatic']
)

////////////////////////////////////////////////////
/* --          INCLUDE NF-CORE MODULES         -- */
////////////////////////////////////////////////////

include { MULTIQC }                       from './modules/nf-core/software/multiqc'

////////////////////////////////////////////////////
/* --       INCLUDE NF-CORE SUBWORKFLOWS       -- */
////////////////////////////////////////////////////

include { QC_TRIM }                       from './modules/nf-core/subworkflow/qc_trim' addParams(
    fastqc_options:                  modules['fastqc'],
    trimgalore_options:              modules['trimgalore']
)

////////////////////////////////////////////////////
/* --             RUN THE WORKFLOW             -- */
////////////////////////////////////////////////////

workflow {


    ////////////////////////////////////////////////////
    /* --               BUILD INDICES              -- */
    ////////////////////////////////////////////////////

    BUILD_INDICES(
        dbsnp,
        fasta,
        germline_resource,
        known_indels,
        pon,
        step,
        target_bed,
        tools)

    intervals = BUILD_INDICES.out.intervals

    bwa  = params.bwa       ? file(params.bwa)       : BUILD_INDICES.out.bwa
    dict = params.dict      ? file(params.dict)      : BUILD_INDICES.out.dict
    fai  = params.fasta_fai ? file(params.fasta_fai) : BUILD_INDICES.out.fai

    dbsnp_tbi             = params.dbsnp             ? params.dbsnp_index             ? file(params.dbsnp_index)             : BUILD_INDICES.out.dbsnp_tbi                  : file("${params.outdir}/no_file")
    germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? file(params.germline_resource_index) : BUILD_INDICES.out.germline_resource_tbi      : file("${params.outdir}/no_file")
    known_indels_tbi      = params.known_indels      ? params.known_indels_index      ? file(params.known_indels_index)      : BUILD_INDICES.out.known_indels_tbi.collect() : file("${params.outdir}/no_file")
    pon_tbi               = params.pon               ? params.pon_index               ? file(params.pon_index)               : BUILD_INDICES.out.pon_tbi                    : file("${params.outdir}/no_file")

    msisensor_scan    = BUILD_INDICES.out.msisensor_scan
    target_bed_gz_tbi = BUILD_INDICES.out.target_bed_gz_tbi

    ////////////////////////////////////////////////////
    /* --               PREPROCESSING              -- */
    ////////////////////////////////////////////////////

    bam_mapped          = Channel.empty()
    bam_mapped_qc       = Channel.empty()
    bam_recalibrated_qc = Channel.empty()
    input_reads         = Channel.empty()
    qc_reports          = Channel.empty()

    // STEP 0: QC & TRIM
    // `--skip_qc fastqc` to skip fastqc
    // trim only with `--trim_fastq`
    // additional options to be set up

    QC_TRIM(
        input_sample,
        ('fastqc' in skip_qc || step != "mapping"),
        !(params.trim_fastq))

    reads_input = QC_TRIM.out.reads

    qc_reports = qc_reports.mix(
        QC_TRIM.out.fastqc_html,
        QC_TRIM.out.fastqc_zip,
        QC_TRIM.out.trimgalore_html,
        QC_TRIM.out.trimgalore_log,
        QC_TRIM.out.trimgalore_zip)

    // STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA-MEM

    MAPPING(
        ('bamqc' in skip_qc),
        ('samtools' in skip_qc),
        bwa,
        fai,
        fasta,
        reads_input,
        save_bam_mapped,
        step,
        target_bed)

    bam_mapped    = MAPPING.out.bam
    bam_mapped_qc = MAPPING.out.qc

    qc_reports = qc_reports.mix(bam_mapped_qc)

    // STEP 2: MARKING DUPLICATES

    MARKDUPLICATES(
        bam_mapped,
        step)

    bam_markduplicates = MARKDUPLICATES.out.bam

    if (step == 'preparerecalibration') bam_markduplicates = input_sample

    // STEP 3: CREATING RECALIBRATION TABLES

    PREPARE_RECALIBRATION(
        bam_markduplicates,
        dbsnp,
        dbsnp_tbi,
        dict,
        fai,
        fasta,
        intervals,
        known_indels,
        known_indels_tbi,
        step)

    table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr

    // STEP 4: RECALIBRATING
    bam_applybqsr = bam_markduplicates.join(table_bqsr)

    if (step == 'recalibrate') bam_applybqsr = input_sample

    RECALIBRATE(
        ('bamqc' in skip_qc),
        ('samtools' in skip_qc),
        bam_applybqsr,
        dict,
        fai,
        fasta,
        intervals,
        step,
        target_bed)

    bam_recalibrated    = RECALIBRATE.out.bam
    bam_recalibrated_qc = RECALIBRATE.out.qc

    qc_reports = qc_reports.mix(bam_recalibrated_qc)

    bam_variant_calling = bam_recalibrated

    if (step == 'variantcalling') bam_variant_calling = input_sample

    ////////////////////////////////////////////////////
    /* --         GERMLINE VARIANT CALLING         -- */
    ////////////////////////////////////////////////////

    GERMLINE_VARIANT_CALLING(
        bam_variant_calling,
        dbsnp,
        dbsnp_tbi,
        dict,
        fai,
        fasta,
        intervals,
        target_bed,
        target_bed_gz_tbi,
        tools)

    ////////////////////////////////////////////////////
    /* --          SOMATIC VARIANT CALLING         -- */
    ////////////////////////////////////////////////////

    // TUMOR_VARIANT_CALLING(
    //     bam_variant_calling,
    //     dbsnp,
    //     dbsnp_tbi,
    //     dict,
    //     fai,
    //     fasta,
    //     intervals,
    //     target_bed,
    //     target_bed_gz_tbi,
    //     tools)

    PAIR_VARIANT_CALLING(
        bam_variant_calling,
        dbsnp,
        dbsnp_tbi,
        dict,
        fai,
        fasta,
        intervals,
        msisensor_scan,
        target_bed,
        target_bed_gz_tbi,
        tools)

    ////////////////////////////////////////////////////
    /* --                ANNOTATION                -- */
    ////////////////////////////////////////////////////

}
