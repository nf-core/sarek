/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSarek.initialise(params, log)

// Check input path parameters to see if they exist
checkPathParamList = [
    params.ac_loci,
    params.ac_loci_gc,
    params.cadd_indels,
    params.cadd_indels_tbi,
    params.cadd_wg_snvs,
    params.cadd_wg_snvs_tbi,
    params.chr_dir,
    params.chr_length,
    params.dbsnp,
    params.fasta,
    params.fasta_fai,
    params.germline_resource,
    params.input,
    params.known_indels,
    params.mappability,
    params.multiqc_config,
    params.pon,
    params.snpeff_cache,
    params.target_bed,
    params.vep_cache
]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
input_sample = Channel.empty()

if (params.input) csv_file = file(params.input)
else {
    log.warn "No samplesheet specified, attempting to restart from csv files present in ${params.outdir}"
    switch (params.step.toLowerCase()) {
        case 'mapping': break
        case 'prepare_recalibration': csv_file = file("${params.outdir}/preprocessing/csv/markduplicates_no_table.csv", checkIfExists: true); break
        case 'recalibrate':           csv_file = file("${params.outdir}/preprocessing/csv/markduplicates.csv", checkIfExists: true); break
        case 'variant_calling':       csv_file = file("${params.outdir}/preprocessing/csv/recalibrated.csv", checkIfExists: true); break
        // case 'controlfreec':          csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
        case 'annotate': break
        default: exit 1, "Unknown step ${params.step}"
    }
}

input_sample = extract_csv(csv_file)

save_bam_mapped = params.skip_markduplicates ? true : params.save_bam_mapped ? true : false

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

////////////////////////////////////////////////////
/* --  UPDATE MODULES OPTIONS BASED ON PARAMS  -- */
////////////////////////////////////////////////////

def modules = params.modules.clone()

if (params.save_reference)       modules['build_intervals'].publish_files         = ['bed':'intervals']
if (params.save_reference)       modules['bwa_index'].publish_files               = ['amb':'bwa', 'ann':'bwa', 'bwt':'bwa', 'pac':'bwa', 'sa':'bwa']
if (params.save_reference)       modules['bwamem2_index'].publish_files           = ['0123':'bwamem2', 'amb':'bwamem2', 'ann':'bwamem2', 'bwt.2bit.64':'bwamem2', 'bwt.8bit.32':'bwamem2', 'pac':'bwamem2']
if (params.save_reference)       modules['create_intervals_bed'].publish_files    = ['bed':'intervals']
if (params.save_reference)       modules['dict'].publish_files                    = ['dict':'dict']
if (params.save_reference)       modules['index_target_bed'].publish_files        = ['bed.gz':'target', 'bed.gz.tbi':'target']
if (params.save_reference)       modules['msisensorpro_scan'].publish_files       = ['list':'msi']
if (params.save_reference)       modules['samtools_faidx'].publish_files          = ['fai':'fai']
if (params.save_reference)       modules['tabix_dbsnp'].publish_files             = ['vcf.gz.tbi':'dbsnp']
if (params.save_reference)       modules['tabix_germline_resource'].publish_files = ['vcf.gz.tbi':'germline_resource']
if (params.save_reference)       modules['tabix_known_indels'].publish_files      = ['vcf.gz.tbi':'known_indels']
if (params.save_reference)       modules['tabix_pon'].publish_files               = ['vcf.gz.tbi':'pon']
if (save_bam_mapped)             modules['samtools_index_mapping'].publish_files  = ['bam':'mapped', 'bai':'mapped']
if (params.skip_markduplicates)  modules['baserecalibrator'].publish_files        = ['recal.table':'mapped']
if (params.skip_markduplicates)  modules['gatherbqsrreports'].publish_files       = ['recal.table':'mapped']
if (!params.skip_markduplicates) modules['baserecalibrator'].publish_files        = false

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
chr_dir           = params.chr_dir           ? file(params.chr_dir)           : ch_dummy_file
chr_length        = params.chr_length        ? file(params.chr_length)        : ch_dummy_file
dbsnp             = params.dbsnp             ? file(params.dbsnp)             : ch_dummy_file
fasta             = params.fasta             ? file(params.fasta)             : ch_dummy_file
fasta_fai         = params.fasta_fai         ? file(params.fasta_fai)         : ch_dummy_file
germline_resource = params.germline_resource ? file(params.germline_resource) : ch_dummy_file
known_indels      = params.known_indels      ? file(params.known_indels)      : ch_dummy_file
loci              = params.ac_loci           ? file(params.ac_loci)           : ch_dummy_file
loci_gc           = params.ac_loci_gc        ? file(params.ac_loci_gc)        : ch_dummy_file
mappability       = params.mappability       ? file(params.mappability)       : ch_dummy_file

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
snpeff_db         = params.snpeff_db         ?: Channel.empty()
snpeff_species    = params.species           ?: Channel.empty()
vep_cache_version = params.vep_cache_version ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
cadd_indels       = params.cadd_indels       ? file(params.cadd_indels)      : ch_dummy_file
cadd_indels_tbi   = params.cadd_indels_tbi   ? file(params.cadd_indels_tbi)  : ch_dummy_file
cadd_wg_snvs      = params.cadd_wg_snvs      ? file(params.cadd_wg_snvs)     : ch_dummy_file
cadd_wg_snvs_tbi  = params.cadd_wg_snvs_tbi  ? file(params.cadd_wg_snvs_tbi) : ch_dummy_file
pon               = params.pon               ? file(params.pon)              : ch_dummy_file
snpeff_cache      = params.snpeff_cache      ? file(params.snpeff_cache)     : ch_dummy_file
target_bed        = params.target_bed        ? file(params.target_bed)       : ch_dummy_file
vep_cache         = params.vep_cache         ? file(params.vep_cache)        : ch_dummy_file

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
read_structure1   = params.read_structure1   ?: Channel.empty()
read_structure2   = params.read_structure2   ?: Channel.empty()

////////////////////////////////////////////////////
/* --        INCLUDE LOCAL SUBWORKFLOWS        -- */
////////////////////////////////////////////////////

include { BUILD_INDICES } from '../subworkflows/local/build_indices' addParams(
    build_intervals_options:         modules['build_intervals'],
    bwa_index_options:               modules['bwa_index'],
    bwamem2_index_options:           modules['bwamem2_index'],
    create_intervals_bed_options:    modules['create_intervals_bed'],
    gatk4_dict_options:              modules['dict'],
    index_target_bed_options:        modules['index_target_bed'],
    msisensorpro_scan_options:       modules['msisensorpro_scan'],
    samtools_faidx_options:          modules['samtools_faidx'],
    tabix_dbsnp_options:             modules['tabix_dbsnp'],
    tabix_germline_resource_options: modules['tabix_germline_resource'],
    tabix_known_indels_options:      modules['tabix_known_indels'],
    tabix_pon_options:               modules['tabix_pon']
)
include { MAPPING } from '../subworkflows/nf-core/mapping' addParams(
    seqkit_split2_options:           modules['seqkit_split2'],
    bwamem1_mem_options:             modules['bwa_mem1_mem'],
    bwamem1_mem_tumor_options:       modules['bwa_mem1_mem_tumor'],
    bwamem2_mem_options:             modules['bwa_mem2_mem'],
    bwamem2_mem_tumor_options:       modules['bwa_mem2_mem_tumor'],
)

include { MAPPING_CSV } from '../subworkflows/local/mapping_csv'
include { MARKDUPLICATES } from '../subworkflows/nf-core/markduplicates' addParams(
    markduplicates_options:             modules['markduplicates'],
    markduplicatesspark_options:        modules['markduplicatesspark'],
    estimatelibrarycomplexity_options:  modules['estimatelibrarycomplexity']
)
include { MARKDUPLICATES_CSV } from '../subworkflows/local/markduplicates_csv'

include { CONVERT_BAM_TO_CRAM } from '../subworkflows/nf-core/convert_bam_to_cram.nf' addParams(
    merge_bam_options:               modules['merge_bam_mapping'],
    samtools_index_options:          modules['samtools_index_mapping'],
    qualimap_bamqc_options:          modules['qualimap_bamqc_mapping'],
    samtools_stats_options:          modules['samtools_stats_mapping'],
    samtools_view_options:           modules['samtools_view']
)

include { PREPARE_RECALIBRATION } from '../subworkflows/nf-core/prepare_recalibration' addParams(
    baserecalibrator_options:        modules['baserecalibrator'],
    gatherbqsrreports_options:       modules['gatherbqsrreports']
)
include { PREPARE_RECALIBRATION_CSV } from '../subworkflows/local/prepare_recalibration_csv'
include { RECALIBRATE } from '../subworkflows/nf-core/recalibrate' addParams(
    applybqsr_options:               modules['applybqsr'],
    merge_cram_options:              modules['merge_cram_recalibrate'],
    qualimap_bamqc_options:          modules['qualimap_bamqc_recalibrate'],
    samtools_index_options:          modules['samtools_index_recalibrate'],
    samtools_stats_options:          modules['samtools_stats_recalibrate']
)
include { RECALIBRATE_CSV } from '../subworkflows/local/recalibrate_csv'
include { GERMLINE_VARIANT_CALLING } from '../subworkflows/local/germline_variant_calling' addParams(
    concat_gvcf_options:             modules['concat_gvcf'],
    concat_haplotypecaller_options:  modules['concat_haplotypecaller'],
    genotypegvcf_options:            modules['genotypegvcf'],
    haplotypecaller_options:         modules['haplotypecaller'],
    strelka_options:                 modules['strelka_germline']
)
// include { TUMOR_VARIANT_CALLING } from '../subworkflows/local/tumor_variant_calling' addParams(
// )
include { PAIR_VARIANT_CALLING } from '../subworkflows/local/pair_variant_calling' addParams(
    manta_options:                   modules['manta_somatic'],
    msisensorpro_msi_options:        modules['msisensorpro_msi'],
    strelka_bp_options:              modules['strelka_somatic_bp'],
    strelka_options:                 modules['strelka_somatic']
)

////////////////////////////////////////////////////
/* --          INCLUDE NF-CORE MODULES         -- */
////////////////////////////////////////////////////

include { MULTIQC }                       from '../modules/nf-core/software/multiqc/main'

////////////////////////////////////////////////////
/* --       INCLUDE NF-CORE SUBWORKFLOWS       -- */
////////////////////////////////////////////////////

include { FASTQC_TRIMGALORE }             from '../subworkflows/nf-core/fastqc_trimgalore' addParams(
    fastqc_options:                  modules['fastqc'],
    trimgalore_options:              modules['trimgalore']
)

workflow SAREK {

    ////////////////////////////////////////////////////
    /* --               BUILD INDICES              -- */
    ////////////////////////////////////////////////////

    BUILD_INDICES(
        dbsnp,
        fasta,
        fasta_fai,
        germline_resource,
        known_indels,
        pon,
        target_bed)

    intervals = BUILD_INDICES.out.intervals

    bwa  = params.bwa       ? file(params.bwa)       : BUILD_INDICES.out.bwa
    dict = params.dict      ? file(params.dict)      : BUILD_INDICES.out.dict
    fai  = params.fasta_fai ? file(params.fasta_fai) : BUILD_INDICES.out.fai

    dbsnp_tbi             = params.dbsnp             ? params.dbsnp_index             ? file(params.dbsnp_index)             : BUILD_INDICES.out.dbsnp_tbi                  : []
    germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? file(params.germline_resource_index) : BUILD_INDICES.out.germline_resource_tbi      : []
    known_indels_tbi      = params.known_indels      ? params.known_indels_index      ? file(params.known_indels_index)      : BUILD_INDICES.out.known_indels_tbi.collect() : []
    pon_tbi               = params.pon               ? params.pon_index               ? file(params.pon_index)               : BUILD_INDICES.out.pon_tbi                    : []

    known_sites     = [dbsnp, known_indels]
    known_sites_tbi = dbsnp_tbi.mix(known_indels_tbi).collect()

    msisensorpro_scan = BUILD_INDICES.out.msisensorpro_scan
    target_bed_gz_tbi = BUILD_INDICES.out.target_bed_gz_tbi

    ////////////////////////////////////////////////////
    /* --               PREPROCESSING              -- */
    ////////////////////////////////////////////////////

    bam_mapped          = Channel.empty()
    bam_mapped_qc       = Channel.empty()
    bam_recalibrated_qc = Channel.empty()
    qc_reports          = Channel.empty()

    // STEP 0: QC & TRIM
    // `--skip_qc fastqc` to skip fastqc
    // trim only with `--trim_fastq`
    // additional options to be set up

    if (params.step == 'mapping') {
        FASTQC_TRIMGALORE(
            input_sample,
            !(params.trim_fastq))

        reads_input = FASTQC_TRIMGALORE.out.reads

        qc_reports = qc_reports.mix(
            FASTQC_TRIMGALORE.out.fastqc_html,
            FASTQC_TRIMGALORE.out.fastqc_zip,
            FASTQC_TRIMGALORE.out.trim_html,
            FASTQC_TRIMGALORE.out.trim_log,
            FASTQC_TRIMGALORE.out.trim_zip)

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        MAPPING(
            params.aligner,
            bwa,
            fai,
            fasta,
            reads_input)

        bam_mapped    = MAPPING.out.bam
//         bam_mapped_qc = MAPPING.out.qc
//         qc_reports = qc_reports.mix(bam_mapped_qc)

        //TODO: is this fine to not have here, if no MD done, then the crams will be provided, but if MD is done then only crams after MD
        // Create CSV to restart from this step
        // MAPPING_CSV(bam_mapped, save_bam_mapped, params.skip_markduplicates)

        if (params.skip_markduplicates) {
            bam_markduplicates = bam_mapped
        } else {
            // STEP 2: MARKING DUPLICATES
            MARKDUPLICATES(bam_mapped, params.use_gatk_spark, !('markduplicates' in params.skip_qc), fasta, fai, dict)
            bam_markduplicates = MARKDUPLICATES.out.bam

            // Create CSV to restart from this step
            // MARKDUPLICATES_CSV(bam_markduplicates)
        }
        CONVERT_BAM_TO_CRAM(
            'bamqc' in params.skip_qc,
            'samtools' in params.skip_qc,
            fasta,
            bam_markduplicates,
            target_bed
        )

        cram_markduplicates = CONVERT_BAM_TO_CRAM.out.cram

        //TODO: Check with Maxime and add this
        //CRAM_CSV(cram_markduplicates)

        qc_reports = qc_reports.mix(CONVERT_BAM_TO_CRAM.out.qc)
    }

    if (params.step.toLowerCase() == 'prepare_recalibration') cram_markduplicates = input_sample

    if (params.step.toLowerCase() in ['mapping', 'prepare_recalibration']) {
        // STEP 3: CREATING RECALIBRATION TABLES
        PREPARE_RECALIBRATION(
            cram_markduplicates,
            dict,
            fai,
            fasta,
            intervals,
            known_sites,
            known_sites_tbi,
            params.no_intervals)

        table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr
        PREPARE_RECALIBRATION_CSV(table_bqsr)

        cram_applybqsr = cram_markduplicates.join(table_bqsr)
    }

    if (params.step.toLowerCase() == 'recalibrate') cram_applybqsr = input_sample

    if (params.step.toLowerCase() in ['mapping', 'prepare_recalibration', 'recalibrate']) {
        // STEP 4: RECALIBRATING
        RECALIBRATE(
            ('bamqc' in params.skip_qc),
            ('samtools' in params.skip_qc),
            cram_applybqsr,
            dict,
            fai,
            fasta,
            intervals,
            target_bed)

        cram_recalibrated    = RECALIBRATE.out.cram
        cram_recalibrated_qc = RECALIBRATE.out.qc

        RECALIBRATE_CSV(cram_recalibrated)

        qc_reports = qc_reports.mix(cram_recalibrated_qc)

        cram_variant_calling = cram_recalibrated

        cram_variant_calling.dump(tag:'input')
    }

    if (params.step.toLowerCase() == 'variant_calling') cram_variant_calling = input_sample

    if (params.tools != null) {

        ////////////////////////////////////////////////////
        /* --         GERMLINE VARIANT CALLING         -- */
        ////////////////////////////////////////////////////
        target_bed_gz_tbi.dump(tag:"tbi_gz")
        GERMLINE_VARIANT_CALLING(
            cram_variant_calling,
            dbsnp,
            dbsnp_tbi,
            dict,
            fai,
            fasta,
            intervals,
            target_bed,
            target_bed_gz_tbi)

        ////////////////////////////////////////////////////
        /* --          SOMATIC VARIANT CALLING         -- */
        ////////////////////////////////////////////////////

        // TUMOR_VARIANT_CALLING(
        //     cram_variant_calling,
        //     dbsnp,
        //     dbsnp_tbi,
        //     dict,
        //     fai,
        //     fasta,
        //     intervals,
        //     target_bed,
        //     target_bed_gz_tbi)

        PAIR_VARIANT_CALLING(
            cram_variant_calling,
            dbsnp,
            dbsnp_tbi,
            dict,
            fai,
            fasta,
            intervals,
            msisensorpro_scan,
            target_bed,
            target_bed_gz_tbi)

        ////////////////////////////////////////////////////
        /* --                ANNOTATION                -- */
        ////////////////////////////////////////////////////
    }
}

def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true).map{ row ->
        def meta = [:]

        meta.patient = row.patient.toString()
        meta.sample  = row.sample.toString()

        // If no gender specified, gender is not considered (only used for somatic CNV)
        if (row.gender == null) {
            meta.gender = "NA"
        } else meta.gender = row.gender.toString()

        // If no status specified, sample is considered normal
        if (row.status == null) {
            meta.status = 0
        } else meta.status = row.status.toInteger()

        if (row.lane != null) {
        // mapping with fastq
            meta.id         = "${row.sample}-${row.lane}".toString()
            def read1       = file(row.fastq1, checkIfExists: true)
            def read2       = file(row.fastq2, checkIfExists: true)
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:ILLUMINA\""
            meta.read_group = read_group.toString()
            return [meta, [read1, read2]]
        } else if (row.table != null) {
        // recalibration
            meta.id = meta.sample
            def bam     = file(row.bam, checkIfExists: true)
            def bai     = file(row.bai, checkIfExists: true)
            def table   = file(row.table, checkIfExists: true)
            return [meta, bam, bai, table]
        } else {
        // prepare_recalibration or variant_calling
            meta.id = meta.sample
            def bam     = file(row.bam, checkIfExists: true)
            def bai     = file(row.bai, checkIfExists: true)
            return [meta, bam, bai]
        }
    }
}