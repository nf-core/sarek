#!/usr/bin/env nextflow

/*
================================================================================
                                  nf-core/sarek
================================================================================
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

// Print help message if required

if (params.help) {
    def command = "nextflow run nf-core/sarek -profile docker --input sample.tsv"
    log.info Schema.params_help("$baseDir/nextflow_schema.json", command)
    exit 0
}

/*
================================================================================
                        INCLUDE SAREK FUNCTIONS
================================================================================
*/

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

/*
================================================================================
                         SET UP CONFIGURATION VARIABLES
================================================================================
*/

// Check parameters

Checks.aws_batch(workflow, params)     // Check AWS batch settings
Checks.hostname(workflow, params, log) // Check the hostnames against configured profiles

// MultiQC - Stage config files

multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)

// Check if genome exists in the config file
if (params.genomes && !params.genomes.containsKey(params.genome) && !params.igenomes_ignore) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
} else if (params.genomes && !params.genomes.containsKey(params.genome) && params.igenomes_ignore) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

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

/*
================================================================================
                     UPDATE MODULES OPTIONS BASED ON PARAMS
================================================================================
*/

modules = params.modules

if (params.save_reference)      modules['build_intervals'].publish_files         = ['bed':'intervals']
if (params.save_reference)      modules['bwa_index'].publish_files               = ['amb':'bwa', 'ann':'bwa', 'bwt':'bwa', 'pac':'bwa', 'sa':'bwa']
if (params.save_reference)      modules['bwamem2_index'].publish_files           = ['0123':'bwamem2', 'amb':'bwamem2', 'ann':'bwamem2', 'bwt.2bit.64':'bwamem2', 'bwt.8bit.32':'bwamem2', 'pac':'bwamem2']
if (params.save_reference)      modules['create_intervals_bed'].publish_files    = ['bed':'intervals']
if (params.save_reference)      modules['dict'].publish_files                    = ['dict':'dict']
if (params.save_reference)      modules['samtools_faidx'].publish_files          = ['fai':'fai']
if (params.save_reference)      modules['tabix_dbsnp'].publish_files             = ['vcf.gz.tbi':'dbsnp']
if (params.save_reference)      modules['tabix_germline_resource'].publish_files = ['vcf.gz.tbi':'germline_resource']
if (params.save_reference)      modules['tabix_known_indels'].publish_files      = ['vcf.gz.tbi':'known_indels']
if (params.save_reference)      modules['tabix_pon'].publish_files               = ['vcf.gz.tbi':'pon']
if (save_bam_mapped)            modules['samtools_index_mapping'].publish_files  = ['bam':'mapped', 'bai':'mapped']
if (params.skip_markduplicates) modules['baserecalibrator'].publish_files        = ['recal.table':'mapped']
if (params.skip_markduplicates) modules['gatherbqsrreports'].publish_files       = ['recal.table':'mapped']

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
params.ac_loci                 = params.genome ? params.genomes[params.genome].ac_loci                 ?: false : false
params.ac_loci_gc              = params.genome ? params.genomes[params.genome].ac_loci_gc              ?: false : false
params.bwa                     = params.genome ? params.genomes[params.genome].bwa                     ?: false : false
params.chr_dir                 = params.genome ? params.genomes[params.genome].chr_dir                 ?: false : false
params.chr_length              = params.genome ? params.genomes[params.genome].chr_length              ?: false : false
params.dbsnp                   = params.genome ? params.genomes[params.genome].dbsnp                   ?: false : false
params.dbsnp_index             = params.genome ? params.genomes[params.genome].dbsnp_index             ?: false : false
params.dict                    = params.genome ? params.genomes[params.genome].dict                    ?: false : false
params.fasta                   = params.genome ? params.genomes[params.genome].fasta                   ?: false : false
params.fasta_fai               = params.genome ? params.genomes[params.genome].fasta_fai               ?: false : false
params.germline_resource       = params.genome ? params.genomes[params.genome].germline_resource       ?: false : false
params.germline_resource_index = params.genome ? params.genomes[params.genome].germline_resource_index ?: false : false
params.intervals               = params.genome ? params.genomes[params.genome].intervals               ?: false : false
params.known_indels            = params.genome ? params.genomes[params.genome].known_indels            ?: false : false
params.known_indels_index      = params.genome ? params.genomes[params.genome].known_indels_index      ?: false : false
params.mappability             = params.genome ? params.genomes[params.genome].mappability             ?: false : false
params.snpeff_db               = params.genome ? params.genomes[params.genome].snpeff_db               ?: false : false
params.species                 = params.genome ? params.genomes[params.genome].species                 ?: false : false
params.vep_cache_version       = params.genome ? params.genomes[params.genome].vep_cache_version       ?: false : false

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

/*
================================================================================
                                PRINTING SUMMARY
================================================================================
*/

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}
summary = Schema.params_summary(workflow, params, run_name, step, tools, skip_qc, annotate_tools)
log.info Headers.nf_core(workflow, params.monochrome_logs)
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
log.info "-\033[2m----------------------------------------------------\033[0m-"

// params summary for MultiQC
workflow_summary = Schema.params_mqc_summary(summary)
workflow_summary = Channel.value(workflow_summary)

if ('mutect2' in tools && !(params.pon)) log.warn "[nf-core/sarek] Mutect2 was requested, but as no panel of normals were given, results will not be optimal"
if (params.sentieon) log.warn "[nf-core/sarek] Sentieon will be used, only works if Sentieon is available where nf-core/sarek is run"

/*
================================================================================
                              INCLUDE LOCAL MODULES
================================================================================
*/

/*
================================================================================
                           INCLUDE LOCAL SUBWORKFLOWS
================================================================================
*/

include { BUILD_INDICES } from './modules/local/subworkflow/build_indices' addParams(
    build_intervals_options:         modules['build_intervals'],
    bwa_index_options:               modules['bwa_index'],
    bwamem2_index_options:           modules['bwamem2_index'],
    create_intervals_bed_options:    modules['create_intervals_bed'],
    gatk_dict_options:               modules['dict'],
    samtools_faidx_options:          modules['samtools_faidx'],
    tabix_dbsnp_options:             modules['tabix_dbsnp'],
    tabix_germline_resource_options: modules['tabix_germline_resource'],
    tabix_known_indels_options:      modules['tabix_known_indels'],
    tabix_pon_options:               modules['tabix_pon']
)

include { SPLIT_FASTQ } from './modules/local/process/seqkit_split2' addParams( options: modules['seqkit'] )

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
    haplotypecaller_options:         modules['haplotypecaller'],
    genotypegvcf_options:            modules['genotypegvcf'],
    concat_gvcf_options:             modules['concat_gvcf'],
    concat_haplotypecaller_options:  modules['concat_haplotypecaller'],
    strelka_options:                 modules['strelka_germline']
)

/*
================================================================================
                             INCLUDE nf-core MODULES
================================================================================
*/

include { MULTIQC }                       from './modules/nf-core/software/multiqc'

/*
================================================================================
                          INCLUDE nf-core SUBWORKFLOWS
================================================================================
*/

include { QC_TRIM }                       from './modules/nf-core/subworkflow/qc_trim' addParams(
    fastqc_options:                  modules['fastqc'],
    trimgalore_options:              modules['trimgalore']
)
// PREPARING CHANNELS FOR PREPROCESSING AND QC

// input_bam = Channel.empty()
// input_pair_reads = Channel.empty()

// if (step in ['preparerecalibration', 'recalibrate', 'variantcalling', 'controlfreec', 'annotate']) {
//     input_bam.close()
//     input_pair_reads.close()
// } else input_sample.branch(input_pair_reads, input_bam) {has_extension(it[3], "bam") ? 1 : 0}

// (input_bam, input_bam_fastqc) = input_bam.into(2)

// // Removing inputFile2 which is null in case of uBAM
// input_bam_fastqc = input_bam_fastqc.map {
//     idPatient, idSample, idRun, inputFile1, inputFile2 ->
//     [idPatient, idSample, idRun, inputFile1]
// }

// if (params.split_fastq){
//    input_pair_reads = input_pair_reads
//        // newly splitfastq are named based on split, so the name is easier to catch
//        .splitFastq(by: params.split_fastq, compress:true, file:"split", pe:true)
//        .map {idPatient, idSample, idRun, reads1, reads2 ->
//            // The split fastq read1 is the 4th element (indexed 3) its name is split_3
//            // The split fastq read2's name is split_4
//            // It's followed by which split it's acutally based on the mother fastq file
//            // Index start at 1
//            // Extracting the index to get a new IdRun
//            splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
//            newIdRun = idRun + "_" + splitIndex
//            // Giving the files a new nice name
//            newReads1 = file("${idSample}_${newIdRun}_R1.fastq.gz")
//            newReads2 = file("${idSample}_${newIdRun}_R2.fastq.gz")
//            [idPatient, idSample, newIdRun, reads1, reads2]}
//}

// input_pair_reads.dump(tag:'INPUT')

// (input_pair_reads, input_pair_readstrimgalore, input_pair_readsfastqc) = input_pair_reads.into(3)


/*
================================================================================
                        RUN THE WORKFLOW
================================================================================
*/

workflow {

/*
================================================================================
                                  BUILD INDICES
================================================================================
*/

    BUILD_INDICES(
        dbsnp,
        fasta,
        germline_resource,
        known_indels,
        pon,
        step,
        tools)

    intervals = BUILD_INDICES.out.intervals

    bwa  = params.bwa       ? file(params.bwa)       : BUILD_INDICES.out.bwa
    dict = params.dict      ? file(params.dict)      : BUILD_INDICES.out.dict
    fai  = params.fasta_fai ? file(params.fasta_fai) : BUILD_INDICES.out.fai

    dbsnp_tbi             = params.dbsnp             ? params.dbsnp_index             ? file(params.dbsnp_index)             : BUILD_INDICES.out.dbsnp_tbi                  : file("${params.outdir}/no_file")
    germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? file(params.germline_resource_index) : BUILD_INDICES.out.germline_resource_tbi      : file("${params.outdir}/no_file")
    known_indels_tbi      = params.known_indels      ? params.known_indels_index      ? file(params.known_indels_index)      : BUILD_INDICES.out.known_indels_tbi.collect() : file("${params.outdir}/no_file")
    pon_tbi               = params.pon               ? params.pon_index               ? file(params.pon_index)               : BUILD_INDICES.out.pon_tbi                    : file("${params.outdir}/no_file")

/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

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

    //reads_input = QC_TRIM.out.reads

    qc_reports = qc_reports.mix(
        QC_TRIM.out.fastqc_html,
        QC_TRIM.out.fastqc_zip,
        QC_TRIM.out.trimgalore_html,
        QC_TRIM.out.trimgalore_log,
        QC_TRIM.out.trimgalore_zip)

    // STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA-MEM

    //TODO: not sure whether to split the reads before qc or after
    SPLIT_FASTQ(QC_TRIM.out.reads)
    reads_input = SPLIT_FASTQ.out.map{
            key, reads ->
                //TODO maybe this can be replaced by a regex to include part_001 etc.

                //sorts list of split fq files by :
                //[R1.part_001, R2.part_001, R1.part_002, R2.part_002,R1.part_003, R2.part_003,...]
                //TODO: determine whether it is possible to have an uneven number of parts, so remainder: true woud need to be used
                return [key, reads.sort{ a,b -> a.getName().tokenize('.')[ a.getName().tokenize('.').size() - 3] <=> b.getName().tokenize('.')[ b.getName().tokenize('.').size() - 3]}
                                        .collate(2)]
        }.transpose()

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

    /*
    ================================================================================
                                GERMLINE VARIANT CALLING
    ================================================================================
    */

    GERMLINE_VARIANT_CALLING(
        bam_variant_calling,
        dbsnp,
        dbsnp_tbi,
        dict,
        fai,
        fasta,
        intervals,
        target_bed,
        tools)

    /*
    ================================================================================
                                SOMATIC VARIANT CALLING
    ================================================================================
    */

    /*
    ================================================================================
                                    ANNOTATION
    ================================================================================
    */


    /*
    ================================================================================
                                        MultiQC
    ================================================================================
    */

    // GET_SOFTWARE_VERSIONS()

    MULTIQC(
        // GET_SOFTWARE_VERSIONS.out.yml,
        multiqc_config,
        multiqc_custom_config.ifEmpty([]),
        workflow_summary,
        qc_reports.collect())
}

/*
================================================================================
                        SEND COMPLETION EMAIL
================================================================================
 */

workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary, run_name, baseDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}

// /*
// ================================================================================
//                             GERMLINE VARIANT CALLING
// ================================================================================
// */

// // STEP MANTA.1 - SINGLE MODE

// process MantaSingle {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSample}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bamMantaSingle
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaSingle

//     when: 'manta' in tools

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     status = status_map[idPatient, idSample]
//     input_bam = status == 0 ? "--bam" : "--tumorBam"
//     vcftype = status == 0 ? "diploid" : "tumor"
//     """
//     ${beforeScript}
//     configManta.py \
//         ${input_bam} ${bam} \
//         --reference ${fasta} \
//         ${options} \
//         --runDir Manta

//     python Manta/runWorkflow.py -m local -j ${task.cpus}

//     mv Manta/results/variants/candidateSmallIndels.vcf.gz \
//         Manta_${idSample}.candidateSmallIndels.vcf.gz
//     mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
//         Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
//     mv Manta/results/variants/candidateSV.vcf.gz \
//         Manta_${idSample}.candidateSV.vcf.gz
//     mv Manta/results/variants/candidateSV.vcf.gz.tbi \
//         Manta_${idSample}.candidateSV.vcf.gz.tbi
//     mv Manta/results/variants/${vcftype}SV.vcf.gz \
//         Manta_${idSample}.${vcftype}SV.vcf.gz
//     mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
//         Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
//     """
// }

// vcfMantaSingle = vcfMantaSingle.dump(tag:'Single Manta')

// // STEP TIDDIT

// process TIDDIT {
//     tag "${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
//             else "Reports/${idSample}/TIDDIT/${it}"
//         }

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bamTIDDIT
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set val("TIDDIT"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi") into vcfTIDDIT
//         set file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig") into tidditOut

//     when: 'tiddit' in tools

//     script:
//     """
//     tiddit --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}

//     mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf

//     grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf

//     bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz

//     tabix TIDDIT_${idSample}.vcf.gz
//     """
// }

// vcfTIDDIT = vcfTIDDIT.dump(tag:'TIDDIT')

// // STEP FREEBAYES SINGLE MODE

// process FreebayesSingle {
//     tag "${idSample}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamFreebayesSingle
//         file(fasta) from fasta
//         file(fastaFai) from ch_software_versions_yaml

//     output:
//         set val("FreeBayes"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfFreebayesSingle

//     when: 'freebayes' in tools

//     script:
//     intervalsOptions = params.no_intervals ? "" : "-t ${intervalBed}"
//     """
//     freebayes \
//         -f ${fasta} \
//         --min-alternate-fraction 0.1 \
//         --min-mapping-quality 1 \
//         ${intervalsOptions} \
//         ${bam} > ${intervalBed.baseName}_${idSample}.vcf
//     """
// }

// vcfFreebayesSingle = vcfFreebayesSingle.groupTuple(by: [0,1,2])

// /*
// ================================================================================
//                              SOMATIC VARIANT CALLING
// ================================================================================
// */
// // Ascat, pileup, pileups with no intervals, recalibrated BAMs
// (bamAscat, bamMpileup, bamMpileupNoInt, bamRecalAll) = bamRecalAll.into(4)

// // separate BAM by status
// bamNormal = Channel.create()
// bamTumor = Channel.create()

// bamRecalAll
//     .choice(bamTumor, bamNormal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// // Remapping channel to remove common key idPatient
// pairBam = bamNormal.cross(bamTumor).map {
//     normal, tumor ->
//     [normal[0], normal[1], normal[2], normal[3], tumor[1], tumor[2], tumor[3]]
// }

// pairBam = pairBam.dump(tag:'BAM Somatic Pair')

// // Manta, Strelka, Mutect2, MSIsensor
// (pairBamManta, pairBamStrelka, pairBamStrelkaBP, pairBamCalculateContamination, pairBamFilterMutect2, pairBamMsisensor, pairBamCNVkit, pairBam) = pairBam.into(8)

// // Making Pair Bam for Sention

// // separate BAM by status
// bam_sention_normal = Channel.create()
// bam_sentieon_tumor = Channel.create()

// bam_sentieon_all
//     .choice(bam_sentieon_tumor, bam_sention_normal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// // Remapping channel to remove common key idPatient

// bam_pair_sentieon_TNscope = bam_sention_normal.cross(bam_sentieon_tumor).map {
//     normal, tumor ->
//     [normal[0], normal[1], normal[2], normal[3], normal[4], tumor[1], tumor[2], tumor[3], tumor[4]]
// }

// intervalPairBam = pairBam.spread(bedIntervals)

// bamMpileup = bamMpileup.spread(intMpileup)

// // intervals for Mutect2 calls, FreeBayes and pileups for Mutect2 filtering
// (pairBamMutect2, pairBamFreeBayes, pairBamPileupSummaries) = intervalPairBam.into(3)

// // STEP FREEBAYES

// process FreeBayes {
//     tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamFreeBayes
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set val("FreeBayes"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into vcfFreeBayes

//     when: 'freebayes' in tools

//     script:
//     intervalsOptions = params.no_intervals ? "" : "-t ${intervalBed}"
//     """
//     freebayes \
//         -f ${fasta} \
//         --pooled-continuous \
//         --pooled-discrete \
//         --genotype-qualities \
//         --report-genotype-likelihood-max \
//         --allele-balance-priors-off \
//         --min-alternate-fraction 0.03 \
//         --min-repeat-entropy 1 \
//         --min-alternate-count 2 \
//         ${intervalsOptions} \
//         ${bamTumor} \
//         ${bamNormal} > ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
//     """
// }

// vcfFreeBayes = vcfFreeBayes.groupTuple(by:[0,1,2])

// // STEP GATK MUTECT2.1 - RAW CALLS

// process Mutect2 {
//     tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamMutect2
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi
//         file(intervals) from intervals
//         file(pon) from pon
//         file(ponIndex) from pon_tbi

//     output:
//         set val("Mutect2"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output
//         set idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats") optional true into intervalStatsFiles
//         set idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats"), file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") optional true into mutect2Stats

//     when: 'mutect2' in tools

//     script:
//     // please make a panel-of-normals, using at least 40 samples
//     // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
//     PON = params.pon ? "--panel-of-normals ${pon}" : ""
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
//     """
//     # Get raw calls
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//       Mutect2 \
//       -R ${fasta}\
//       -I ${bamTumor}  -tumor ${idSampleTumor} \
//       -I ${bamNormal} -normal ${idSampleNormal} \
//       ${intervalsOptions} \
//       ${softClippedOption} \
//       --germline-resource ${germlineResource} \
//       ${PON} \
//       -O ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
//     """
// }

// mutect2Output = mutect2Output.groupTuple(by:[0,1,2])
// mutect2Stats = mutect2Stats.groupTuple(by:[0,1])

// // STEP GATK MUTECT2.2 - MERGING STATS

// process MergeMutect2Stats {
//     tag "${idSamplePair}"

//     publishDir "${params.outdir}/VariantCalling/${idSamplePair}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSamplePair, file(statsFiles), file(vcf) from mutect2Stats // Actual stats files and corresponding VCF chunks
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi
//         file(intervals) from intervals

//     output:
//         set idPatient, idSamplePair, file("${idSamplePair}.vcf.gz.stats") into mergedStatsFile

//     when: 'mutect2' in tools

//     script:
//     stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         MergeMutectStats \
//         ${stats} \
//         -O ${idSamplePair}.vcf.gz.stats
//     """
// }

// // we are merging the VCFs that are called separatelly for different intervals
// // so we can have a single sorted VCF containing all the calls for a given caller

// // STEP MERGING VCF - FREEBAYES & GATK HAPLOTYPECALLER

// vcfConcatenateVCFs = vcfFreeBayes.mix(vcfFreebayesSingle, vcfGenotypeGVCFs, gvcfHaplotypeCaller)
// vcfConcatenateVCFs = vcfConcatenateVCFs.dump(tag:'VCF to merge')

// process ConcatVCF {
//     label 'cpus_8'

//     tag "${variantCaller}-${idSample}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idPatient, idSample, file(vcf) from vcfConcatenateVCFs
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//     // we have this funny *_* pattern to avoid copying the raw calls to publishdir
//         set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated

//     when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

//     script:
//     if (variantCaller == 'HaplotypeCallerGVCF')
//         outputFile = "HaplotypeCaller_${idSample}.g.vcf"
//     else
//         outputFile = "${variantCaller}_${idSample}.vcf"
//     options = params.target_bed ? "-t ${targetBED}" : ""
//     intervalsOptions = params.no_intervals ? "-n" : ""
//     """
//     concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
//     """
// }

// vcfConcatenated = vcfConcatenated.dump(tag:'VCF')

// // STEP MERGING VCF - GATK MUTECT2 (UNFILTERED)

// mutect2Output = mutect2Output.dump(tag:'Mutect2 output VCF to merge')

// process ConcatVCF_Mutect2 {
//     label 'cpus_8'

//     tag "${idSample}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idPatient, idSample, file(vcf) from mutect2Output
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//     // we have this funny *_* pattern to avoid copying the raw calls to publishdir
//         set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenatedForFilter

//     when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

//     script:
//     outputFile = "Mutect2_unfiltered_${idSample}.vcf"
//     options = params.target_bed ? "-t ${targetBED}" : ""
//     intervalsOptions = params.no_intervals ? "-n" : ""
//     """
//     concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
//     """
// }

// vcfConcatenatedForFilter = vcfConcatenatedForFilter.dump(tag:'Mutect2 unfiltered VCF')

// // STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

// pairBamPileupSummaries = pairBamPileupSummaries.map{
//     idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor, intervalBed ->
//     [idPatient, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor, intervalBed]
// }.join(intervalStatsFiles, by:[0,1,2])

// process PileupSummariesForMutect2 {
//     tag "${idSampleTumor}_vs_${idSampleNormal}-${intervalBed.baseName}"

//     label 'cpus_1'

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(bamNormal), file(baiNormal), file(bamTumor), file(baiTumor), file(intervalBed), file(statsFile) from pairBamPileupSummaries
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table") into pileupSummaries

//     when: 'mutect2' in tools

//     script:
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GetPileupSummaries \
//         -I ${bamTumor} \
//         -V ${germlineResource} \
//         ${intervalsOptions} \
//         -O ${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table
//     """
// }

// pileupSummaries = pileupSummaries.groupTuple(by:[0,1,2])

// // STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

// process MergePileupSummaries {
//     label 'cpus_1'

//     tag "${idPatient}_${idSampleTumor}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(pileupSums) from pileupSummaries
//         file(dict) from dict

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}_pileupsummaries.table") into mergedPileupFile

//     when: 'mutect2' in tools

//     script:
//     allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
//     """
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         GatherPileupSummaries \
//         --sequence-dictionary ${dict} \
//         ${allPileups} \
//         -O ${idSampleTumor}_pileupsummaries.table
//     """
// }

// // STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

// pairBamCalculateContamination = pairBamCalculateContamination.map{
//     idPatient, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor ->
//     [idPatient, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor]
// }.join(mergedPileupFile, by:[0,1,2])

// process CalculateContamination {
//     label 'cpus_1'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(bamNormal), file(baiNormal), file(bamTumor), file(baiTumor), file(mergedPileup) from pairBamCalculateContamination

//      output:
//         set idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("${idSampleTumor}_contamination.table") into contaminationTable

//     when: 'mutect2' in tools

//     script:
//     """
//     # calculate contamination
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         CalculateContamination \
//         -I ${idSampleTumor}_pileupsummaries.table \
//         -O ${idSampleTumor}_contamination.table
//     """
// }

// // STEP GATK MUTECT2.6 - FILTERING CALLS

// mutect2CallsToFilter = vcfConcatenatedForFilter.map{
//     variantCaller, idPatient, idSamplePair, vcf, tbi ->
//     [idPatient, idSamplePair, vcf, tbi]
// }.join(mergedStatsFile, by:[0,1]).join(contaminationTable, by:[0,1])

// process FilterMutect2Calls {
//     label 'cpus_1'

//     tag "${idSamplePair}"

//     publishDir "${params.outdir}/VariantCalling/${idSamplePair}/Mutect2", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSamplePair, file(unfiltered), file(unfilteredIndex), file(stats), file(contaminationTable) from mutect2CallsToFilter
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(germlineResource) from germline_resource
//         file(germlineResourceIndex) from germline_resource_tbi
//         file(intervals) from intervals

//     output:
//         set val("Mutect2"), idPatient, idSamplePair, file("Mutect2_filtered_${idSamplePair}.vcf.gz"), file("Mutect2_filtered_${idSamplePair}.vcf.gz.tbi"), file("Mutect2_filtered_${idSamplePair}.vcf.gz.filteringStats.tsv") into filteredMutect2Output

//     when: 'mutect2' in tools

//     script:
//     """
//     # do the actual filtering
//     gatk --java-options "-Xmx${task.memory.toGiga()}g" \
//         FilterMutectCalls \
//         -V ${unfiltered} \
//         --contamination-table ${contaminationTable} \
//         --stats ${stats} \
//         -R ${fasta} \
//         -O Mutect2_filtered_${idSamplePair}.vcf.gz
//     """
// }

// // STEP SENTIEON TNSCOPE

// process Sentieon_TNscope {
//     label 'cpus_max'
//     label 'memory_max'
//     label 'sentieon'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), file(recalNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(recalTumor) from bam_pair_sentieon_TNscope
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(pon) from pon
//         file(ponIndex) from pon_tbi

//     output:
//         set val("SentieonTNscope"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf") into vcf_sentieon_TNscope

//     when: 'tnscope' in tools && params.sentieon

//     script:
//     PON = params.pon ? "--pon ${pon}" : ""
//     """
//     sentieon driver \
//         -t ${task.cpus} \
//         -r ${fasta} \
//         -i ${bamTumor} \
//         -q ${recalTumor} \
//         -i ${bamNormal} \
//         -q ${recalNormal} \
//         --algo TNscope \
//         --tumor_sample ${idSampleTumor} \
//         --normal_sample ${idSampleNormal} \
//         --dbsnp ${dbsnp} \
//         ${PON} \
//         TNscope_${idSampleTumor}_vs_${idSampleNormal}.vcf
//     """
// }

// vcf_sentieon_TNscope = vcf_sentieon_TNscope.dump(tag:'Sentieon TNscope')

// vcf_sentieon = vcf_sentieon_DNAseq.mix(vcf_sentieon_DNAscope, vcf_sentieon_DNAscope_SV, vcf_sentieon_TNscope)

// process CompressSentieonVCF {
//     tag "${idSample} - ${vcf}"

//     publishDir "${params.outdir}/VariantCalling/${idSample}/${variantCaller}", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idPatient, idSample, file(vcf) from vcf_sentieon

//     output:
//         set variantCaller, idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcf_sentieon_compressed

//     script:
//     """
//     bgzip < ${vcf} > ${vcf}.gz
//     tabix ${vcf}.gz
//     """
// }

// vcf_sentieon_compressed = vcf_sentieon_compressed.dump(tag:'Sentieon VCF indexed')

// // STEP STRELKA.2 - SOMATIC PAIR

// process Strelka {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Strelka", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamStrelka
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Strelka"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelka

//     when: 'strelka' in tools

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     """
//     ${beforeScript}
//     configureStrelkaSomaticWorkflow.py \
//         --tumor ${bamTumor} \
//         --normal ${bamNormal} \
//         --referenceFasta ${fasta} \
//         ${options} \
//         --runDir Strelka

//     python Strelka/runWorkflow.py -m local -j ${task.cpus}

//     mv Strelka/results/variants/somatic.indels.vcf.gz \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
//     mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
//     mv Strelka/results/variants/somatic.snvs.vcf.gz \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
//     mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
//         Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
//     """
// }

// vcfStrelka = vcfStrelka.dump(tag:'Strelka')

// // STEP MANTA.2 - SOMATIC PAIR

// process Manta {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Manta", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamManta
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Manta"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfManta
//         set idPatient, idSampleNormal, idSampleTumor, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

//     when: 'manta' in tools

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     """
//     ${beforeScript}
//     configManta.py \
//         --normalBam ${bamNormal} \
//         --tumorBam ${bamTumor} \
//         --reference ${fasta} \
//         ${options} \
//         --runDir Manta

//     python Manta/runWorkflow.py -m local -j ${task.cpus}

//     mv Manta/results/variants/candidateSmallIndels.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
//     mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
//     mv Manta/results/variants/candidateSV.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz
//     mv Manta/results/variants/candidateSV.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
//     mv Manta/results/variants/diploidSV.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz
//     mv Manta/results/variants/diploidSV.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
//     mv Manta/results/variants/somaticSV.vcf.gz \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz
//     mv Manta/results/variants/somaticSV.vcf.gz.tbi \
//         Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
//     """
// }

// vcfManta = vcfManta.dump(tag:'Manta')

// // Remmaping channels to match input for StrelkaBP
// pairBamStrelkaBP = pairBamStrelkaBP.map {
//     idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor ->
//     [idPatientNormal, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor]
// }.join(mantaToStrelka, by:[0,1,2]).map {
//     idPatientNormal, idSampleNormal, idSampleTumor, bamNormal, baiNormal, bamTumor, baiTumor, mantaCSI, mantaCSIi ->
//     [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor, mantaCSI, mantaCSIi]
// }

// // STEP STRELKA.3 - SOMATIC PAIR - BEST PRACTICES

// process StrelkaBP {
//     label 'cpus_max'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Strelka", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(mantaCSI), file(mantaCSIi) from pairBamStrelkaBP
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai
//         file(targetBED) from ch_target_bed

//     output:
//         set val("Strelka"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelkaBP

//     when: 'strelka' in tools && 'manta' in tools && !params.no_strelka_bp

//     script:
//     beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
//     options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
//     """
//     ${beforeScript}
//     configureStrelkaSomaticWorkflow.py \
//         --tumor ${bamTumor} \
//         --normal ${bamNormal} \
//         --referenceFasta ${fasta} \
//         --indelCandidates ${mantaCSI} \
//         ${options} \
//         --runDir Strelka

//     python Strelka/runWorkflow.py -m local -j ${task.cpus}

//     mv Strelka/results/variants/somatic.indels.vcf.gz \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
//     mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
//     mv Strelka/results/variants/somatic.snvs.vcf.gz \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
//     mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
//         StrelkaBP_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
//     """
// }

// vcfStrelkaBP = vcfStrelkaBP.dump(tag:'Strelka BP')

// // STEP CNVkit

// process CNVkit {
//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/CNVkit", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamCNVkit
//         file(targetBED) from ch_target_bed
//         file(fasta) from fasta

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}*"), file("${idSampleNormal}*") into cnvkitOut

//     when: 'cnvkit' in tools && params.target_bed

//     script:
//     """
//     cnvkit.py \
//       batch \
//       ${bamTumor} \
//       --normal ${bamNormal} \
//       --targets ${targetBED} \
//       --fasta ${fasta} \
//       --output-reference output_reference.cnn \
//       --output-dir ./ \
//       --diagram \
//       --scatter
//     """
// }

// // STEP MSISENSOR.1 - SCAN

// // Scan reference genome for microsatellites
// process MSIsensor_scan {
//     label 'cpus_1'
//     label 'memory_max'

//     tag "${fasta}"

//     input:
//     file(fasta) from fasta
//     file(fastaFai) from fai

//     output:
//     file "microsatellites.list" into msi_scan_ch

//     when: 'msisensor' in tools

//     script:
//     """
//     msisensor scan -d ${fasta} -o microsatellites.list
//     """
// }

// // STEP MSISENSOR.2 - SCORE

// // Score the normal vs somatic pair of bams

// process MSIsensor_msi {
//     label 'cpus_4'
//     label 'memory_max'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/MSIsensor", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamMsisensor
//         file msiSites from msi_scan_ch

//     output:
//         set val("Msisensor"), idPatient, file("${idSampleTumor}_vs_${idSampleNormal}_msisensor"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_dis"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_germline"), file("${idSampleTumor}_vs_${idSampleNormal}_msisensor_somatic") into msisensor_out_ch

//     when: 'msisensor' in tools

//     script:
//     """
//     msisensor msi -d ${msiSites} \
//                   -b 4 \
//                   -n ${bamNormal} \
//                   -t ${bamTumor} \
//                   -o ${idSampleTumor}_vs_${idSampleNormal}_msisensor
//     """
// }

// // STEP ASCAT.1 - ALLELECOUNTER

// // Run commands and code from Malin Larsson
// // Based on Jesper Eisfeldt's code
// process AlleleCounter {
//     label 'memory_singleCPU_2_task'

//     tag "${idSample}"

//     input:
//         set idPatient, idSample, file(bam), file(bai) from bamAscat
//         file(acLoci) from loci
//         file(dict) from dict
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSample, file("${idSample}.alleleCount") into alleleCounterOut

//     when: 'ascat' in tools

//     script:
//     """
//     alleleCounter \
//         -l ${acLoci} \
//         -r ${fasta} \
//         -b ${bam} \
//         -o ${idSample}.alleleCount;
//     """
// }

// alleleCountOutNormal = Channel.create()
// alleleCountOutTumor = Channel.create()

// alleleCounterOut
//     .choice(alleleCountOutTumor, alleleCountOutNormal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// alleleCounterOut = alleleCountOutNormal.combine(alleleCountOutTumor, by:0)

// alleleCounterOut = alleleCounterOut.map {
//     idPatientNormal, idSampleNormal, alleleCountOutNormal,
//     idSampleTumor, alleleCountOutTumor ->
//     [idPatientNormal, idSampleNormal, idSampleTumor, alleleCountOutNormal, alleleCountOutTumor]
// }

// // STEP ASCAT.2 - CONVERTALLELECOUNTS

// // R script from Malin Larssons bitbucket repo:
// // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
// process ConvertAlleleCounts {
//     label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCounterOut

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOut

//     when: 'ascat' in tools

//     script:
//     gender = gender_map[idPatient]
//     """
//     convertAlleleCounts.r ${idSampleTumor} ${alleleCountTumor} ${idSampleNormal} ${alleleCountNormal} ${gender}
//     """
// }

// // STEP ASCAT.3 - ASCAT

// // R scripts from Malin Larssons bitbucket repo:
// // https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
// process Ascat {
//     label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOut
//         file(acLociGC) from loci_gc

//     output:
//         set val("ASCAT"), idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.*.{png,txt}") into ascatOut

//     when: 'ascat' in tools

//     script:
//     gender = gender_map[idPatient]
//     purity_ploidy = (params.ascat_purity && params.ascat_ploidy) ? "--purity ${params.ascat_purity} --ploidy ${params.ascat_ploidy}" : ""
//     """
//     for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
//     run_ascat.r \
//         --tumorbaf ${bafTumor} \
//         --tumorlogr ${logrTumor} \
//         --normalbaf ${bafNormal} \
//         --normallogr ${logrNormal} \
//         --tumorname ${idSampleTumor} \
//         --basedir ${baseDir} \
//         --gcfile ${acLociGC} \
//         --gender ${gender} \
//         ${purity_ploidy}
//     """
// }

// ascatOut.dump(tag:'ASCAT')

// // STEP MPILEUP.1

// process Mpileup {
//     label 'cpus_1'
//     label 'memory_singleCPU_2_task'

//     tag "${idSample}-${intervalBed.baseName}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup" ? "VariantCalling/${idSample}/Control-FREEC/${it}" : null }

//     input:
//         set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamMpileup
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSample, file("${prefix}${idSample}.pileup") into mpileupMerge
//         set idPatient, idSample into tsv_mpileup

//     when: 'controlfreec' in tools || 'mpileup' in tools

//     script:
//     prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
//     intervalsOptions = params.no_intervals ? "" : "-l ${intervalBed}"

//     """
//     # Control-FREEC reads uncompresses the zipped file TWICE in single-threaded mode.
//     # we are therefore not using compressed pileups here
//     samtools mpileup \
//         -f ${fasta} ${bam} \
//         ${intervalsOptions} > ${prefix}${idSample}.pileup
//     """
// }

// (tsv_mpileup, tsv_mpileup_sample) = tsv_mpileup.groupTuple(by:[0, 1]).into(2)

// // Creating a TSV file to restart from this step
// tsv_mpileup.map { idPatient, idSample ->
//     gender = gender_map[idPatient]
//     status = status_map[idPatient, idSample]
//     mpileup = "${params.outdir}/VariantCalling/${idSample}/Control-FREEC/${idSample}.pileup"
//     "${idPatient}\t${gender}\t${status}\t${idSample}\t${mpileup}\n"
// }.collectFile(
//     name: 'control-freec_mpileup.tsv', sort: true, storeDir: "${params.outdir}/VariantCalling/TSV"
// )

// tsv_mpileup_sample
//     .collectFile(storeDir: "${params.outdir}/VariantCalling/TSV") {
//         idPatient, idSample ->
//         status = status_map[idPatient, idSample]
//         gender = gender_map[idPatient]
//         mpileup = "${params.outdir}/VariantCalling/${idSample}/Control-FREEC/${idSample}.pileup"
//         ["control-freec_mpileup_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${mpileup}\n"]
// }

// if (!params.no_intervals) {
//     mpileupMerge = mpileupMerge.groupTuple(by:[0, 1])
//     mpileupNoInt = Channel.empty()
// } else {
//     (mpileupMerge, mpileupNoInt) = mpileupMerge.into(2)
//     mpileupMerge.close()
// }

// // STEP MPILEUP.2 - MERGE
// process MergeMpileup {
//     label 'cpus_1'

//     tag "${idSample}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: { it == "${idSample}.pileup" ? "VariantCalling/${idSample}/Control-FREEC/${it}" : null }

//     input:
//         set idPatient, idSample, file(mpileup) from mpileupMerge

//     output:
//         set idPatient, idSample, file("${idSample}.pileup") into mpileupOut

//     when: !(params.no_intervals) && 'controlfreec' in tools || 'mpileup' in tools

//     script:
//     """
//     for i in `ls -1v *.pileup`;
//         do cat \$i >> ${idSample}.pileup
//     done
//     """
// }

// mpileupOut = mpileupOut.mix(mpileupNoInt)
// mpileupOut = mpileupOut.dump(tag:'mpileup')

// mpileupOutNormal = Channel.create()
// mpileupOutTumor = Channel.create()

// if (step == 'controlfreec') mpileupOut = input_sample

// mpileupOut
//     .choice(mpileupOutTumor, mpileupOutNormal) {status_map[it[0], it[1]] == 0 ? 1 : 0}

// mpileupOut = mpileupOutNormal.combine(mpileupOutTumor, by:0)

// mpileupOut = mpileupOut.map {
//     idPatientNormal, idSampleNormal, mpileupOutNormal,
//     idSampleTumor, mpileupOutTumor ->
//     [idPatientNormal, idSampleNormal, idSampleTumor, mpileupOutNormal, mpileupOutTumor]
// }

// // STEP CONTROLFREEC.1 - CONTROLFREEC

// process ControlFREEC {
//     label 'cpus_max'
//     //label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Control-FREEC", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(mpileupNormal), file(mpileupTumor) from mpileupOut
//         file(chrDir) from chr_dir
//         file(mappability) from mappability
//         file(chrLength) from chr_length
//         file(dbsnp) from dbsnp
//         file(dbsnpIndex) from dbsnp_tbi
//         file(fasta) from fasta
//         file(fastaFai) from fai

//     output:
//         set idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.pileup_CNVs"), file("${idSampleTumor}.pileup_ratio.txt"), file("${idSampleTumor}.pileup_normal_CNVs"), file("${idSampleTumor}.pileup_normal_ratio.txt"), file("${idSampleTumor}.pileup_BAF.txt"), file("${idSampleNormal}.pileup_BAF.txt") into controlFreecViz
//         set file("*.pileup*"), file("${idSampleTumor}_vs_${idSampleNormal}.config.txt") into controlFreecOut

//     when: 'controlfreec' in tools

//     script:
//     config = "${idSampleTumor}_vs_${idSampleNormal}.config.txt"
//     gender = gender_map[idPatient]
//     // if we are using coefficientOfVariation, we must delete the window parameter
//     // it is "window = 20000" in the default settings, without coefficientOfVariation set,
//     // but we do not like it. Note, it is not written in stone
//     coeff_or_window = params.cf_window ? "window = ${params.cf_window}" : "coefficientOfVariation = ${params.cf_coeff}"

//     """
//     touch ${config}
//     echo "[general]" >> ${config}
//     echo "BedGraphOutput = TRUE" >> ${config}
//     echo "chrFiles = \${PWD}/${chrDir.fileName}" >> ${config}
//     echo "chrLenFile = \${PWD}/${chrLength.fileName}" >> ${config}
//     echo "gemMappabilityFile = \${PWD}/${mappability}" >> ${config}
//     echo "${coeff_or_window}" >> ${config}
//     echo "contaminationAdjustment = TRUE" >> ${config}
//     echo "forceGCcontentNormalization = 1" >> ${config}
//     echo "maxThreads = ${task.cpus}" >> ${config}
//     echo "minimalSubclonePresence = 20" >> ${config}
//     echo "ploidy = ${params.cf_ploidy}" >> ${config}
//     echo "sex = ${gender}" >> ${config}
//     echo "" >> ${config}

//     echo "[control]" >> ${config}
//     echo "inputFormat = pileup" >> ${config}
//     echo "mateFile = \${PWD}/${mpileupNormal}" >> ${config}
//     echo "mateOrientation = FR" >> ${config}
//     echo "" >> ${config}

//     echo "[sample]" >> ${config}
//     echo "inputFormat = pileup" >> ${config}
//     echo "mateFile = \${PWD}/${mpileupTumor}" >> ${config}
//     echo "mateOrientation = FR" >> ${config}
//     echo "" >> ${config}

//     echo "[BAF]" >> ${config}
//     echo "SNPfile = ${dbsnp.fileName}" >> ${config}

//     freec -conf ${config}
//     """
// }

// controlFreecOut.dump(tag:'ControlFREEC')

// // STEP CONTROLFREEC.3 - VISUALIZATION

// process ControlFreecViz {
//     label 'memory_singleCPU_2_task'

//     tag "${idSampleTumor}_vs_${idSampleNormal}"

//     publishDir "${params.outdir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Control-FREEC", mode: params.publish_dir_mode

//     input:
//         set idPatient, idSampleNormal, idSampleTumor, file(cnvTumor), file(ratioTumor), file(cnvNormal), file(ratioNormal), file(bafTumor), file(bafNormal) from controlFreecViz

//     output:
//         set file("*.txt"), file("*.png"), file("*.bed") into controlFreecVizOut

//     when: 'controlfreec' in tools

//     """
//     echo "Shaping CNV files to make sure we can assess significance"
//     awk 'NF==9{print}' ${cnvTumor} > TUMOR.CNVs
//     awk 'NF==7{print}' ${cnvNormal} > NORMAL.CNVs

//     echo "############### Calculating significance values for TUMOR CNVs #############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/assess_significance.R | R --slave --args TUMOR.CNVs ${ratioTumor}

//     echo "############### Calculating significance values for NORMAL CNVs ############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/assess_significance.R | R --slave --args NORMAL.CNVs ${ratioNormal}

//     echo "############### Creating graph for TUMOR ratios ###############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/makeGraph.R | R --slave --args 2 ${ratioTumor} ${bafTumor}

//     echo "############### Creating graph for NORMAL ratios ##############"
//     cat /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/makeGraph.R | R --slave --args 2 ${ratioNormal} ${bafNormal}

//     echo "############### Creating BED files for TUMOR ##############"
//     perl /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/freec2bed.pl -f ${ratioTumor} > ${idSampleTumor}.bed

//     echo "############### Creating BED files for NORMAL #############"
//     perl /opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/freec2bed.pl -f ${ratioNormal} > ${idSampleNormal}.bed
//     """
// }

// controlFreecVizOut.dump(tag:'ControlFreecViz')

// // Remapping channels for QC and annotation

// (vcfStrelkaIndels, vcfStrelkaSNVS) = vcfStrelka.into(2)
// (vcfStrelkaBPIndels, vcfStrelkaBPSNVS) = vcfStrelkaBP.into(2)
// (vcfMantaSomaticSV, vcfMantaDiploidSV) = vcfManta.into(2)

// vcfKeep = Channel.empty().mix(
//     filteredMutect2Output.map{
//         variantCaller, idPatient, idSample, vcf, tbi, tsv ->
//         [variantcaller, idSample, vcf]
//     },
//     vcfConcatenated.map{
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf]
//     },
//     vcf_sentieon_compressed.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf]
//     },
//     vcfStrelkaSingle.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[1]]
//     },
//     vcfMantaSingle.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[2]]
//     },
//     vcfMantaDiploidSV.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[2]]
//     },
//     vcfMantaSomaticSV.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[3]]
//     },
//     vcfStrelkaIndels.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[0]]
//     },
//     vcfStrelkaSNVS.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[1]]
//     },
//     vcfStrelkaBPIndels.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[0]]
//     },
//     vcfStrelkaBPSNVS.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf[1]]
//     },
//     vcfTIDDIT.map {
//         variantcaller, idPatient, idSample, vcf, tbi ->
//         [variantcaller, idSample, vcf]
//     })

// (vcfBCFtools, vcfVCFtools, vcfAnnotation) = vcfKeep.into(3)

// // STEP VCF.QC

// process BcftoolsStats {
//     label 'cpus_1'

//     tag "${variantCaller} - ${vcf}"

//     publishDir "${params.outdir}/Reports/${idSample}/BCFToolsStats", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from vcfBCFtools

//     output:
//         file ("*.bcf.tools.stats.out") into bcftoolsReport

//     when: !('bcftools' in skip_qc)

//     script:
//     """
//     bcftools stats ${vcf} > ${reduceVCF(vcf.fileName)}.bcf.tools.stats.out
//     """
// }

// bcftoolsReport = bcftoolsReport.dump(tag:'BCFTools')

// process Vcftools {
//     label 'cpus_1'

//     tag "${variantCaller} - ${vcf}"

//     publishDir "${params.outdir}/Reports/${idSample}/VCFTools", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from vcfVCFtools

//     output:
//         file ("${reduceVCF(vcf.fileName)}.*") into vcftoolsReport

//     when: !('vcftools' in skip_qc)

//     script:
//     """
//     vcftools \
//     --gzvcf ${vcf} \
//     --TsTv-by-count \
//     --out ${reduceVCF(vcf.fileName)}

//     vcftools \
//     --gzvcf ${vcf} \
//     --TsTv-by-qual \
//     --out ${reduceVCF(vcf.fileName)}

//     vcftools \
//     --gzvcf ${vcf} \
//     --FILTER-summary \
//     --out ${reduceVCF(vcf.fileName)}
//     """
// }

// vcftoolsReport = vcftoolsReport.dump(tag:'VCFTools')

// /*
// ================================================================================
//                                    ANNOTATION
// ================================================================================
// */

// if (step == 'annotate') {
//     vcfToAnnotate = Channel.create()
//     vcfNoAnnotate = Channel.create()

//     if (tsv_path == []) {
//     // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
//     // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
//     // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
//     // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
//     // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
//     // This field is used to output final annotated VCFs in the correct directory
//       Channel.empty().mix(
//         Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
//           .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
//           .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/Mutect2/*.vcf.gz")
//           .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAseq/*.vcf.gz")
//           .flatten().map{vcf -> ['SentieonDNAseq', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAscope/*.vcf.gz")
//           .flatten().map{vcf -> ['SentieonDNAscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonTNscope/*.vcf.gz")
//           .flatten().map{vcf -> ['SentieonTNscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
//           .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
//         Channel.fromPath("${params.outdir}/VariantCalling/*/TIDDIT/*.vcf.gz")
//           .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
//       ).choice(vcfToAnnotate, vcfNoAnnotate) {
//         annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools) ? 0 : 1
//       }
//     } else if (annotate_tools == []) {
//     // Annotate user-submitted VCFs
//     // If user-submitted, Sarek assume that the idSample should be assumed automatically
//       vcfToAnnotate = Channel.fromPath(tsv_path)
//         .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
//     } else exit 1, "specify only tools or files to annotate, not both"

//     vcfNoAnnotate.close()
//     vcfAnnotation = vcfAnnotation.mix(vcfToAnnotate)
// }

// // as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

// (vcfSnpeff, vcfVep) = vcfAnnotation.into(2)

// vcfVep = vcfVep.map {
//   variantCaller, idSample, vcf ->
//   [variantCaller, idSample, vcf, null]
// }

// // STEP SNPEFF

// process Snpeff {
//     tag "${idSample} - ${variantCaller} - ${vcf}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
//         if (it == "${reducedVCF}_snpEff.ann.vcf") null
//         else "Reports/${idSample}/snpEff/${it}"
//     }

//     input:
//         set variantCaller, idSample, file(vcf) from vcfSnpeff
//         file(dataDir) from snpeff_cache
//         val snpeffDb from snpeff_db

//     output:
//         set file("${reducedVCF}_snpEff.genes.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv") into snpeffReport
//         set variantCaller, idSample, file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCF

//     when: 'snpeff' in tools || 'merge' in tools

//     script:
//     reducedVCF = reduceVCF(vcf.fileName)
//     cache = (params.snpeff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
//     """
//     snpEff -Xmx${task.memory.toGiga()}g \
//         ${snpeffDb} \
//         -csvStats ${reducedVCF}_snpEff.csv \
//         -nodownload \
//         ${cache} \
//         -canon \
//         -v \
//         ${vcf} \
//         > ${reducedVCF}_snpEff.ann.vcf

//     mv snpEff_summary.html ${reducedVCF}_snpEff.html
//     """
// }

// snpeffReport = snpeffReport.dump(tag:'snpEff report')

// // STEP COMPRESS AND INDEX VCF.1 - SNPEFF

// process CompressVCFsnpEff {
//     tag "${idSample} - ${vcf}"

//     publishDir "${params.outdir}/Annotation/${idSample}/snpEff", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from snpeffVCF

//     output:
//         set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into (compressVCFsnpEffOut)

//     script:
//     """
//     bgzip < ${vcf} > ${vcf}.gz
//     tabix ${vcf}.gz
//     """
// }

// compressVCFsnpEffOut = compressVCFsnpEffOut.dump(tag:'VCF')

// // STEP VEP.1

// process VEP {
//     label 'VEP'
//     label 'cpus_4'

//     tag "${idSample} - ${variantCaller} - ${vcf}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
//         if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
//         else null
//     }

//     input:
//         set variantCaller, idSample, file(vcf), file(idx) from vcfVep
//         file(dataDir) from vep_cache
//         val cache_version from vep_cache_version
//         file(cadd_InDels) from cadd_indels
//         file(cadd_InDels_tbi) from cadd_indels_tbi
//         file(cadd_WG_SNVs) from cadd_wg_snvs
//         file(cadd_WG_SNVs_tbi) from cadd_wg_snvs_tbi
//     output:
//         set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCF
//         file("${reducedVCF}_VEP.summary.html") into vepReport

//     when: 'vep' in tools

//     script:
//     reducedVCF = reduceVCF(vcf.fileName)
//     genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome

//     dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
//     cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
//     genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
//     """
//     mkdir ${reducedVCF}

//     vep \
//         -i ${vcf} \
//         -o ${reducedVCF}_VEP.ann.vcf \
//         --assembly ${genome} \
//         --species ${params.species} \
//         ${cadd} \
//         ${genesplicer} \
//         --cache \
//         --cache_version ${cache_version} \
//         --dir_cache ${dir_cache} \
//         --everything \
//         --filter_common \
//         --fork ${task.cpus} \
//         --format vcf \
//         --per_gene \
//         --stats_file ${reducedVCF}_VEP.summary.html \
//         --total_length \
//         --vcf

//     rm -rf ${reducedVCF}
//     """
// }

// vepReport = vepReport.dump(tag:'VEP')

// // STEP VEP.2 - VEP AFTER SNPEFF

// process VEPmerge {
//     label 'VEP'
//     label 'cpus_4'

//     tag "${idSample} - ${variantCaller} - ${vcf}"

//     publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
//         if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
//         else null
//     }

//     input:
//         set variantCaller, idSample, file(vcf), file(idx) from compressVCFsnpEffOut
//         file(dataDir) from vep_cache
//         val cache_version from vep_cache_version
//         file(cadd_InDels) from cadd_indels
//         file(cadd_InDels_tbi) from cadd_indels_tbi
//         file(cadd_WG_SNVs) from cadd_wg_snvs
//         file(cadd_WG_SNVs_tbi) from cadd_wg_snvs_tbi
//     output:
//         set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCFmerge
//         file("${reducedVCF}_VEP.summary.html") into vepReportMerge

//     when: 'merge' in tools

//     script:
//     reducedVCF = reduceVCF(vcf.fileName)
//     genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
//     dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
//     cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
//     genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
//     """
//     mkdir ${reducedVCF}

//     vep \
//         -i ${vcf} \
//         -o ${reducedVCF}_VEP.ann.vcf \
//         --assembly ${genome} \
//         --species ${params.species} \
//         ${cadd} \
//         ${genesplicer} \
//         --cache \
//         --cache_version ${cache_version} \
//         --dir_cache ${dir_cache} \
//         --everything \
//         --filter_common \
//         --fork ${task.cpus} \
//         --format vcf \
//         --per_gene \
//         --stats_file ${reducedVCF}_VEP.summary.html \
//         --total_length \
//         --vcf

//     rm -rf ${reducedVCF}
//     """
// }

// vepReportMerge = vepReportMerge.dump(tag:'VEP')

// vcfCompressVCFvep = vepVCF.mix(vepVCFmerge)

// // STEP COMPRESS AND INDEX VCF.2 - VEP

// process CompressVCFvep {
//     tag "${idSample} - ${vcf}"

//     publishDir "${params.outdir}/Annotation/${idSample}/VEP", mode: params.publish_dir_mode

//     input:
//         set variantCaller, idSample, file(vcf) from vcfCompressVCFvep

//     output:
//         set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into compressVCFOutVEP

//     script:
//     """
//     bgzip < ${vcf} > ${vcf}.gz
//     tabix ${vcf}.gz
//     """
// }

// compressVCFOutVEP = compressVCFOutVEP.dump(tag:'VCF')