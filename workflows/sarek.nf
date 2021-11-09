/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSarek.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.ac_loci,
    params.ac_loci_gc,
    params.cadd_indels,
    params.cadd_indels_tbi,
    params.cadd_wg_snvs,
    params.cadd_wg_snvs_tbi,
    params.chr_dir,
    params.chr_length,
    params.dbsnp,
    params.dbsnp_tbi,
    params.fasta,
    params.fasta_fai,
    params.germline_resource,
    params.germline_resource_tbi,
    params.input,
    params.known_indels,
    params.known_indels_tbi,
    params.mappability,
    params.multiqc_config,
    params.pon,
    params.pon_tbi,
    params.snpeff_cache,
    params.target_bed,
    params.vep_cache
]

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Get step and tools
def step = params.step ? params.step.replaceAll('-', '').replaceAll('_', '') : ''
def tools = params.tools ? params.tools.split(',').collect{ it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '') } : []
def skip_qc = params.skip_qc ? params.skip_qc.split(',').collect{ it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '') } : []

// Check mandatory parameters
if (params.input) csv_file = file(params.input)
else {
    log.warn "No samplesheet specified, attempting to restart from csv files present in ${params.outdir}"
    switch (step) {
        case 'mapping': exit 1, "Can't start with step $step without samplesheet"
        case 'preparerecalibration': csv_file = file("${params.outdir}/preprocessing/csv/markduplicates_no_table.csv", checkIfExists: true); break
        case 'recalibrate':          csv_file = file("${params.outdir}/preprocessing/csv/markduplicates.csv",          checkIfExists: true); break
        case 'variantcalling':       csv_file = file("${params.outdir}/preprocessing/csv/recalibrated.csv",            checkIfExists: true); break
        // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
        case 'annotate':             csv_file = file("${params.outdir}/variant_calling/csv/recalibrated.csv",          checkIfExists: true); break
        default: exit 1, "Unknown step $step"
    }
}

input_sample = extract_csv(csv_file)

def save_bam_mapped = params.skip_markduplicates ? true : params.save_bam_mapped ? true : false

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[params.genome]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Stage dummy file to be used as an optional input where required
ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()

// Don't overwrite global params.modules, create a copy instead and use that within the main script
def modules = params.modules.clone()

// Update modules options based on params
if (params.save_reference) {
    modules['build_intervals'].publish_files         = ['bed':'intervals']
    modules['bwa_index'].publish_files               = ['amb':'bwa', 'ann':'bwa', 'bwt':'bwa', 'pac':'bwa', 'sa':'bwa']
    modules['bwamem2_index'].publish_files           = ['0123':'bwamem2', 'amb':'bwamem2', 'ann':'bwamem2', 'bwt.2bit.64':'bwamem2', 'bwt.8bit.32':'bwamem2', 'pac':'bwamem2']
    modules['create_intervals_bed'].publish_files    = ['bed':'intervals']
    modules['dict'].publish_files                    = ['dict':'dict']
    modules['bgziptabix_target_bed'].publish_files   = ['bed.gz':'target', 'bed.gz.tbi':'target']
    modules['msisensorpro_scan'].publish_files       = ['list':'msi']
    modules['samtools_faidx'].publish_files          = ['fai':'fai']
    modules['tabix_dbsnp'].publish_files             = ['vcf.gz.tbi':'dbsnp']
    modules['tabix_germline_resource'].publish_files = ['vcf.gz.tbi':'germline_resource']
    modules['tabix_known_indels'].publish_files      = ['vcf.gz.tbi':'known_indels']
    modules['tabix_pon'].publish_files               = ['vcf.gz.tbi':'pon']
}
if (params.skip_markduplicates) {
    modules['baserecalibrator'].publish_files        = ['recal.table':'mapped']
    modules['gatherbqsrreports'].publish_files       = ['recal.table':'mapped']
} else {
    modules['baserecalibrator'].publish_files        = false
}

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
chr_dir           = params.chr_dir           ? Channel.fromPath(params.chr_dir).collect()           : ch_dummy_file
chr_length        = params.chr_length        ? Channel.fromPath(params.chr_length).collect()        : ch_dummy_file
dbsnp             = params.dbsnp             ? Channel.fromPath(params.dbsnp).collect()             : ch_dummy_file
fasta             = params.fasta             ? Channel.fromPath(params.fasta).collect()             : ch_dummy_file
germline_resource = params.germline_resource ? Channel.fromPath(params.germline_resource).collect() : ch_dummy_file
known_indels      = params.known_indels      ? Channel.fromPath(params.known_indels).collect()      : ch_dummy_file
loci              = params.ac_loci           ? Channel.fromPath(params.ac_loci).collect()           : ch_dummy_file
loci_gc           = params.ac_loci_gc        ? Channel.fromPath(params.ac_loci_gc).collect()        : ch_dummy_file
mappability       = params.mappability       ? Channel.fromPath(params.mappability).collect()       : ch_dummy_file

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
snpeff_db         = params.snpeff_db         ?: Channel.empty()
vep_cache_version = params.vep_cache_version ?: Channel.empty()
vep_genome        = params.vep_genome        ?: Channel.empty()
vep_species       = params.vep_species       ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
cadd_indels       = params.cadd_indels       ? Channel.fromPath(params.cadd_indels).collect()       : ch_dummy_file
cadd_indels_tbi   = params.cadd_indels_tbi   ? Channel.fromPath(params.cadd_indels_tbi).collect()   : ch_dummy_file
cadd_wg_snvs      = params.cadd_wg_snvs      ? Channel.fromPath(params.cadd_wg_snvs).collect()      : ch_dummy_file
cadd_wg_snvs_tbi  = params.cadd_wg_snvs_tbi  ? Channel.fromPath(params.cadd_wg_snvs_tbi).collect()  : ch_dummy_file
pon               = params.pon               ? Channel.fromPath(params.pon).collect()               : ch_dummy_file
snpeff_cache      = params.snpeff_cache      ? Channel.fromPath(params.snpeff_cache).collect()      : ch_dummy_file
target_bed        = params.target_bed        ? Channel.fromPath(params.target_bed).collect()        : ch_dummy_file
vep_cache         = params.vep_cache         ? Channel.fromPath(params.vep_cache).collect()         : ch_dummy_file

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
read_structure1   = params.read_structure1   ?: Channel.empty()
read_structure2   = params.read_structure2   ?: Channel.empty()

// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules

// Create samplesheets to restart from different steps
include { MAPPING_CSV }               from '../subworkflows/local/mapping_csv'
include { MARKDUPLICATES_CSV }        from '../subworkflows/local/markduplicates_csv'
include { PREPARE_RECALIBRATION_CSV } from '../subworkflows/local/prepare_recalibration_csv'
include { RECALIBRATE_CSV }           from '../subworkflows/local/recalibrate_csv'

// Build indices if needed
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams(
    bgziptabix_target_bed_options:     modules['bgziptabix_target_bed'],
    build_intervals_options:           modules['build_intervals'],
    bwa_index_options:                 modules['bwa_index'],
    bwamem2_index_options:             modules['bwamem2_index'],
    create_intervals_bed_options:      modules['create_intervals_bed'],
    gatk4_dict_options:                modules['dict'],
    msisensorpro_scan_options:         modules['msisensorpro_scan'],
    samtools_faidx_options:            modules['samtools_faidx'],
    tabix_dbsnp_options:               modules['tabix_dbsnp'],
    tabix_germline_resource_options:   modules['tabix_germline_resource'],
    tabix_known_indels_options:        modules['tabix_known_indels'],
    tabix_pon_options:                 modules['tabix_pon']
)

// Map input reads to reference genome (+QC)
include { MAPPING } from '../subworkflows/nf-core/mapping' addParams(
    bwamem1_mem_options:               modules['bwa_mem1_mem'],
    bwamem1_mem_tumor_options:         modules['bwa_mem1_mem_tumor'],
    bwamem2_mem_options:               modules['bwa_mem2_mem'],
    bwamem2_mem_tumor_options:         modules['bwa_mem2_mem_tumor'],
    merge_bam_options:                 modules['merge_bam_mapping'],
    samtools_index_options:            modules['samtools_index_mapping'],
    seqkit_split2_options:             modules['seqkit_split2']
)

// Mark duplicates (+QC) + convert to CRAM
include { MARKDUPLICATES } from '../subworkflows/nf-core/markduplicates' addParams(
    estimatelibrarycomplexity_options: modules['estimatelibrarycomplexity'],
    markduplicates_options:            modules['markduplicates'],
    markduplicatesspark_options:       modules['markduplicatesspark'],
    merge_bam_options:                 modules['merge_bam_mapping'],
    qualimap_bamqc_options:            modules['qualimap_bamqc_mapping'],
    samtools_index_options:            modules['samtools_index_cram'],
    samtools_stats_options:            modules['samtools_stats_mapping'],
    samtools_view_options:             modules['samtools_view']
)

// Create recalibration tables
include { PREPARE_RECALIBRATION } from '../subworkflows/nf-core/prepare_recalibration' addParams(
    baserecalibrator_options:          modules['baserecalibrator'],
    baserecalibrator_spark_options:    modules['baserecalibrator_spark'],
    gatherbqsrreports_options:         modules['gatherbqsrreports']
)

// Create recalibrated cram files to use for variant calling
include { RECALIBRATE } from '../subworkflows/nf-core/recalibrate' addParams(
    applybqsr_options:                 modules['applybqsr'],
    applybqsr_spark_options:           modules['applybqsr_spark'],
    merge_cram_options:                modules['merge_cram_recalibrate'],
    qualimap_bamqc_options:            modules['qualimap_bamqc_recalibrate'],
    samtools_index_options:            modules['samtools_index_recalibrate'],
    samtools_stats_options:            modules['samtools_stats_recalibrate']
)

// Variant calling on a single normal sample
include { GERMLINE_VARIANT_CALLING } from '../subworkflows/local/germline_variant_calling' addParams(
    concat_gvcf_options:               modules['concat_gvcf'],
    concat_haplotypecaller_options:    modules['concat_haplotypecaller'],
    genotypegvcf_options:              modules['genotypegvcf'],
    haplotypecaller_options:           modules['haplotypecaller'],
    strelka_options:                   modules['strelka_germline']
)
// // Variant calling on a single tumor sample
// // include { TUMOR_VARIANT_CALLING } from '../subworkflows/local/tumor_variant_calling' addParams(
// // )
// // Variant calling on tumor/normal pair
// include { PAIR_VARIANT_CALLING } from '../subworkflows/local/pair_variant_calling' addParams(
//     manta_options:                   modules['manta_somatic'],
//     msisensorpro_msi_options:        modules['msisensorpro_msi'],
//     strelka_bp_options:              modules['strelka_somatic_bp'],
//     strelka_options:                 modules['strelka_somatic'],
//     mutect2_somatic_options:         modules['mutect2_somatic']
// )

// Annotation
include { ANNOTATE } from '../subworkflows/local/annotate' addParams(
    annotation_cache:                  params.annotation_cache,
    bgziptabix_merge_vep_options:      modules['bgziptabix_merge_vep'],
    bgziptabix_snpeff_options:         modules['bgziptabix_snpeff'],
    bgziptabix_vep_options:            modules['bgziptabix_vep'],
    merge_vep_options:                 modules['merge_vep'],
    snpeff_options:                    modules['snpeff'],
    snpeff_tag:                        "${modules['snpeff'].tag_base}.${params.genome}",
    vep_options:                       modules['ensemblvep'],
    vep_tag:                           "${modules['ensemblvep'].tag_base}.${params.genome}"
)

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

// Config files
ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
modules['multiqc'].args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// SUBWORKFLOWS
//

include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore' addParams(
    fastqc_options:     modules['fastqc'],
    trimgalore_options: modules['trimgalore']
)

//
// MODULES: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main' addParams( options: [publish_files : ['_versions.yml':'']] )
include { MULTIQC }                     from '../modules/nf-core/modules/multiqc/main'                     addParams(options: modules['multiqc'])

def multiqc_report = []

workflow SAREK {

    ch_versions = Channel.empty()
    qc_reports  = Channel.empty()

    // Build indices if needed
    PREPARE_GENOME(
        dbsnp,
        fasta,
        params.fasta_fai,
        germline_resource,
        known_indels,
        pon,
        target_bed,
        tools,
        step)

    // Gather built indices or get them from the params
    bwa                   = params.fasta             ? params.bwa                   ? Channel.fromPath(params.bwa).collect()                   : PREPARE_GENOME.out.bwa                   : ch_dummy_file
    dict                  = params.fasta             ? params.dict                  ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : ch_dummy_file
    fai                   = params.fasta             ? params.fasta_fai             ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fai                   : ch_dummy_file
    dbsnp_tbi             = params.dbsnp             ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : ch_dummy_file
    germline_resource_tbi = params.germline_resource ? params.germline_resource_tbi ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : ch_dummy_file
    known_indels_tbi      = params.known_indels      ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : ch_dummy_file
    pon_tbi               = params.pon               ? params.pon_tbi               ? Channel.fromPath(params.pon_tbi).collect()               : BUILD_INDICES.out.pon_tbi               : ch_dummy_file
    msisensorpro_scan     = PREPARE_GENOME.out.msisensorpro_scan
    target_bed_gz_tbi     = PREPARE_GENOME.out.target_bed_gz_tbi

    //TODO @Rike, is this working for you?
    // known_sites is made by grouping both the dbsnp and the known indels ressources
    // Which can either or both be optional
    known_sites     = dbsnp.concat(known_indels).collect()
    known_sites_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals     = PREPARE_GENOME.out.intervals
    num_intervals = 0
    intervals.count().map{ num_intervals = it }

    // Get versions from all software used
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    // PREPROCESSING

    bam_mapped          = Channel.empty()
    bam_mapped_qc       = Channel.empty()
    bam_recalibrated_qc = Channel.empty()
    bam_variant_calling = Channel.empty()

    // STEP 0: QC & TRIM
    // `--skip_qc fastqc` to skip fastqc
    // trim only with `--trim_fastq`
    // additional options to be set up

    if (step == 'mapping') {
        FASTQC_TRIMGALORE(
            input_sample,
            ('fastqc' in skip_qc),
            !(params.trim_fastq))

        // Get reads after optional trimming (+QC)
        reads_input = FASTQC_TRIMGALORE.out.reads

        // Get all qc reports for MultiQC
        qc_reports = qc_reports.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{ it[1] }.ifEmpty([]))
        qc_reports = qc_reports.mix(FASTQC_TRIMGALORE.out.trim_log.collect{ it[1] }.ifEmpty([]))
        qc_reports = qc_reports.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{ it[1] }.ifEmpty([]))

        // Get versions from all software used
        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        MAPPING(
            params.aligner,
            bwa,
            fai,
            fasta,
            reads_input,
            params.skip_markduplicates,
            save_bam_mapped)

        // Get mapped reads (BAM) with and without index
        bam_mapped  = MAPPING.out.bam
        bam_indexed = MAPPING.out.bam_indexed

        // Create CSV to restart from this step
        MAPPING_CSV(bam_indexed, save_bam_mapped, params.skip_markduplicates)
    }

    if (step == 'preparerecalibration') {
        if (params.skip_markduplicates) bam_indexed         = input_sample
        else                            cram_markduplicates = input_sample
    }

    if (step in ['mapping', 'preparerecalibration']) {
        // STEP 2: Mark duplicates (+QC) + convert to CRAM
        MARKDUPLICATES(
            bam_mapped,
            bam_indexed,
            ('markduplicates' in params.use_gatk_spark),
            !('markduplicates' in skip_qc),
            fasta,
            fai,
            dict,
            params.skip_markduplicates,
            ('bamqc' in skip_qc),
            ('samtools' in skip_qc),
            target_bed)

        cram_markduplicates = MARKDUPLICATES.out.cram

        // Create CSV to restart from this step
        MARKDUPLICATES_CSV(cram_markduplicates)

        qc_reports = qc_reports.mix(MARKDUPLICATES.out.qc.collect{ it[1] }.ifEmpty([]))

        // STEP 3: Create recalibration tables
        if(!params.skip_bqsr){
            PREPARE_RECALIBRATION(
                cram_markduplicates,
                ('bqsr' in params.use_gatk_spark),
                dict,
                fai,
                fasta,
                intervals,
                num_intervals,
                known_sites,
                known_sites_tbi,
                params.no_intervals)

            table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr
            PREPARE_RECALIBRATION_CSV(table_bqsr)

            cram_applybqsr = cram_markduplicates.join(table_bqsr)
        }
    }

    if (step == 'recalibrate') bam_applybqsr = input_sample

    if (step in ['mapping', 'preparerecalibration', 'recalibrate']) {

        if(!params.skip_bqsr){
            // STEP 4: RECALIBRATING
            RECALIBRATE(
                ('bqsr' in params.use_gatk_spark),
                ('bamqc' in skip_qc),
                ('samtools' in skip_qc),
                cram_applybqsr,
                dict,
                fai,
                fasta,
                intervals,
                num_intervals,
                target_bed)

            cram_recalibrated    = RECALIBRATE.out.cram
            cram_recalibrated_qc = RECALIBRATE.out.qc

            RECALIBRATE_CSV(cram_recalibrated)

            qc_reports = qc_reports.mix(cram_recalibrated_qc.collect{ it[1] }.ifEmpty([]))
            cram_variant_calling = cram_recalibrated

        }else{
            cram_variant_calling = cram_markduplicates
        }

    }

    if (step in 'variantcalling') cram_variant_calling = input_sample

    if (tools != []) {
        vcf_to_annotate = Channel.empty()
        if (step in 'annotate') cram_variant_calling = Channel.empty()

        // GERMLINE VARIANT CALLING
        GERMLINE_VARIANT_CALLING(
            tools,
            cram_variant_calling,
            dbsnp,
            dbsnp_tbi,
            dict,
            fai,
            fasta,
            intervals,
            num_intervals,
            target_bed,
            target_bed_gz_tbi)

        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.haplotypecaller_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.strelka_vcf)

        // SOMATIC VARIANT CALLING

        // TUMOR ONLY VARIANT CALLING
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

        // PAIR VARIANT CALLING
        // PAIR_VARIANT_CALLING(
        //     tools,
        //     cram_variant_calling,
        //     dbsnp,
        //     dbsnp_tbi,
        //     dict,
        //     fai,
        //     fasta,
        //     intervals,
        //     msisensorpro_scan,
        //     target_bed,
        //     target_bed_gz_tbi,
        //     germline_resource,
        //     germline_resource_tbi,
        //     pon,
        //     pon_tbi)

        // ANNOTATE
        if (step == 'annotate') vcf_to_annotate = input_sample

        if ('merge' in tools || 'snpeff' in tools || 'vep' in tools) {
            ANNOTATE(
                vcf_to_annotate,
                tools,
                snpeff_db,
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache)
        }
    }

    ch_version_yaml = Channel.empty()
    if (!('versions' in skip_qc)) {
        CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile())
        ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }

    workflow_summary    = WorkflowSarek.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_config)
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_version_yaml)
    ch_multiqc_files = ch_multiqc_files.mix(qc_reports)

    multiqc_report = Channel.empty()
    if (!('multiqc' in skip_qc)) {
        MULTIQC(ch_multiqc_files.collect())
        multiqc_report = MULTIQC.out.report.toList()
    }
}

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.patient && row.sample)) log.warn "Missing or unknown field in csv file header"
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing
        def meta = [:]

        //TODO since it is mandatory: error/warning if not present?
        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no gender specified, gender is not considered
        // gender is only mandatory for somatic CNV
        if (row.gender) meta.gender = row.gender.toString()
        else meta.gender = "NA"

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        // mapping with fastq
        if (row.lane && row.fastq_2) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:ILLUMINA\""
            meta.numLanes = numLanes.toInteger()
            meta.read_group = read_group.toString()
            return [meta, [fastq_1, fastq_2]]
        // recalibration
        } else if (row.table && row.cram) {
            meta.id   = meta.sample
            def cram  = file(row.cram,  checkIfExists: true)
            def crai  = file(row.crai,  checkIfExists: true)
            def table = file(row.table, checkIfExists: true)
            return [meta, cram, crai, table]
        // recalibration when skipping MarkDuplicates
        } else if (row.table && row.bam) {
            meta.id   = meta.sample
            def bam   = file(row.bam,   checkIfExists: true)
            def bai   = file(row.bai,   checkIfExists: true)
            def table = file(row.table, checkIfExists: true)
            return [meta, bam, bai, table]
        // prepare_recalibration or variant_calling
        } else if (row.cram) {
            meta.id = meta.sample
            def cram = file(row.cram, checkIfExists: true)
            def crai = file(row.crai, checkIfExists: true)
            return [meta, cram, crai]
        // prepare_recalibration when skipping MarkDuplicates
        } else if (row.bam) {
            meta.id = meta.sample
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)
            return [meta, bam, bai]
        // annotation
        } else if (row.vcf) {
            meta.id = meta.sample
            def vcf = file(row.vcf, checkIfExists: true)
            return [meta, vcf]
        } else {
            log.warn "Missing or unknown field in csv file header"
        }
    }
}
