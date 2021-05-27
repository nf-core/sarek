/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def valid_params = [
    aligner : ['bwa-mem', 'bwa-mem2'],
    step : ['annotate', 'controlfreec', 'mapping', 'preparerecalibration', 'recalibrate', 'variantcalling'],
    tools : ['ascat', 'cnvkit', 'controlfreec', 'dnascope', 'dnaseq', 'freebayes', 'haplotypecaller', 'manta', 'merge', 'mpileup', 'msisensor', 'mutect2', 'snpeff', 'strelka', 'tiddit', 'tnscope', 'vep'],
    toolsAnnotate : ['haplotypecaller', 'manta', 'mutect2', 'strelka', 'tiddit'],
    toolsQcSkip : ['bamqc', 'baserecalibrator', 'bcftools', 'documentation', 'fastqc', 'markduplicates', 'multiqc', 'samtools', 'sentieon', 'vcftools', 'versions']
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSarek.initialise(params, log, valid_params)

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
if (params.input) { input_sample = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
chr_dir           = params.chr_dir           ? file(params.chr_dir)           : ch_dummy_file
chr_length        = params.chr_length        ? file(params.chr_length)        : ch_dummy_file
dbsnp             = params.dbsnp             ? file(params.dbsnp)             : ch_dummy_file
fasta             = params.fasta             ? file(params.fasta)             : ch_dummy_file
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
    gatk_dict_options:               modules['dict'],
    index_target_bed_options:        modules['index_target_bed'],
    msisensor_scan_options:          modules['msisensor_scan'],
    samtools_faidx_options:          modules['samtools_faidx'],
    tabix_dbsnp_options:             modules['tabix_dbsnp'],
    tabix_germline_resource_options: modules['tabix_germline_resource'],
    tabix_known_indels_options:      modules['tabix_known_indels'],
    tabix_pon_options:               modules['tabix_pon']
)
include { MAPPING } from '../subworkflows/local/mapping' addParams(
    bwamem1_mem_options:             modules['bwa_mem1_mem'],
    bwamem1_mem_tumor_options:       modules['bwa_mem1_mem_tumor'],
    bwamem2_mem_options:             modules['bwa_mem2_mem'],
    bwamem2_mem_tumor_options:       modules['bwa_mem2_mem_tumor'],
    merge_bam_options:               modules['merge_bam_mapping'],
    qualimap_bamqc_options:          modules['qualimap_bamqc_mapping'],
    samtools_index_options:          modules['samtools_index_mapping'],
    samtools_stats_options:          modules['samtools_stats_mapping']
)
include { MARKDUPLICATES } from '../subworkflows/local/markduplicates' addParams(
    markduplicates_options:          modules['markduplicates']
)
include { PREPARE_RECALIBRATION } from '../subworkflows/local/prepare_recalibration' addParams(
    baserecalibrator_options:        modules['baserecalibrator'],
    gatherbqsrreports_options:       modules['gatherbqsrreports']
)
include { RECALIBRATE } from '../subworkflows/local/recalibrate' addParams(
    applybqsr_options:               modules['applybqsr'],
    merge_bam_options:               modules['merge_bam_recalibrate'],
    qualimap_bamqc_options:          modules['qualimap_bamqc_recalibrate'],
    samtools_index_options:          modules['samtools_index_recalibrate'],
    samtools_stats_options:          modules['samtools_stats_recalibrate']
)
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
    msisensor_msi_options:           modules['msisensor_msi'],
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

    msisensor_scan    = BUILD_INDICES.out.msisensor_scan
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

    FASTQC_TRIMGALORE(
        input_sample,
        ('fastqc' in params.skip_qc || params.step != "mapping"),
        !(params.trim_fastq))

    reads_input = FASTQC_TRIMGALORE.out.reads

    qc_reports = qc_reports.mix(
        FASTQC_TRIMGALORE.out.fastqc_html,
        FASTQC_TRIMGALORE.out.fastqc_zip,
        FASTQC_TRIMGALORE.out.trim_html,
        FASTQC_TRIMGALORE.out.trim_log,
        FASTQC_TRIMGALORE.out.trim_zip)

    // STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA-MEM

    MAPPING(
        ('bamqc' in params.skip_qc),
        ('samtools' in params.skip_qc),
        bwa,
        fai,
        fasta,
        reads_input,
        save_bam_mapped,
        params.step,
        target_bed)

    bam_mapped    = MAPPING.out.bam
    bam_mapped_qc = MAPPING.out.qc

    qc_reports = qc_reports.mix(bam_mapped_qc)

    // STEP 2: MARKING DUPLICATES

    bam_markduplicates = channel.empty()

    if (params.step == 'preparerecalibration') {
        if (params.skip_markduplicates) bam_markduplicates = bam_mapped
        else {
            MARKDUPLICATES(bam_mapped, !('markduplicates' in params.skip_qc))
            bam_markduplicates = MARKDUPLICATES.out.bam
        }
    }

    if (params.step == 'preparerecalibration') bam_markduplicates = input_sample

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
        known_indels_tbi)

    table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr

    // STEP 4: RECALIBRATING
    bam_applybqsr = bam_markduplicates.join(table_bqsr)

    if (params.step == 'recalibrate') bam_applybqsr = input_sample

    RECALIBRATE(
        ('bamqc' in params.skip_qc),
        ('samtools' in params.skip_qc),
        bam_applybqsr,
        dict,
        fai,
        fasta,
        intervals,
        target_bed)

    bam_recalibrated    = RECALIBRATE.out.bam
    bam_recalibrated_qc = RECALIBRATE.out.qc

    qc_reports = qc_reports.mix(bam_recalibrated_qc)

    bam_variant_calling = bam_recalibrated

    if (params.step == 'variantcalling') bam_variant_calling = input_sample

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
        target_bed_gz_tbi)

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
    //     target_bed_gz_tbi)

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
        target_bed_gz_tbi)

    ////////////////////////////////////////////////////
    /* --                ANNOTATION                -- */
    ////////////////////////////////////////////////////

}