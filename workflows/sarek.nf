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
    params.bwa,
    params.cadd_indels,
    params.cadd_indels_tbi,
    params.cadd_wg_snvs,
    params.cadd_wg_snvs_tbi,
    params.chr_dir,
    params.chr_length,
    params.dbsnp,
    params.dbsnp_tbi,
    params.dict,
    params.fasta,
    params.fasta_fai,
    params.germline_resource,
    params.germline_resource_tbi,
    params.input,
    params.intervals,
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

// Check mandatory parameters
if (params.input) csv_file = file(params.input)
else {
    log.warn "No samplesheet specified, attempting to restart from csv files present in ${params.outdir}"
    switch (params.step) {
        case 'mapping': exit 1, "Can't start with step $params.step without samplesheet"
        case 'preparerecalibration': csv_file = file("${params.outdir}/preprocessing/csv/markduplicates_no_table.csv", checkIfExists: true); break
        case 'recalibrate':          csv_file = file("${params.outdir}/preprocessing/csv/markduplicates.csv",          checkIfExists: true); break
        case 'variantcalling':       csv_file = file("${params.outdir}/preprocessing/csv/recalibrated.csv",            checkIfExists: true); break
        // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
        case 'annotate':             csv_file = file("${params.outdir}/variant_calling/csv/recalibrated.csv",          checkIfExists: true); break
        default: exit 1, "Unknown step $params.step"
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

// // Stage dummy file to be used as an optional input where required
// [] = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()

// if (params.skip_markduplicates) {
//     modules['baserecalibrator'].publish_files        = ['recal.table':'mapped']
//     modules['gatherbqsrreports'].publish_files       = ['recal.table':'mapped']
// } else {
//     modules['baserecalibrator'].publish_files        = false
// }

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
chr_dir           = params.chr_dir           ? Channel.fromPath(params.chr_dir).collect()           : []
chr_length        = params.chr_length        ? Channel.fromPath(params.chr_length).collect()        : []
dbsnp             = params.dbsnp             ? Channel.fromPath(params.dbsnp).collect()             : Channel.empty()
fasta             = params.fasta             ? Channel.fromPath(params.fasta).collect()             : []
germline_resource = params.germline_resource ? Channel.fromPath(params.germline_resource).collect() : []
known_indels      = params.known_indels      ? Channel.fromPath(params.known_indels).collect()      : Channel.empty()
loci              = params.ac_loci           ? Channel.fromPath(params.ac_loci).collect()           : []
loci_gc           = params.ac_loci_gc        ? Channel.fromPath(params.ac_loci_gc).collect()        : []
mappability       = params.mappability       ? Channel.fromPath(params.mappability).collect()       : []

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
snpeff_db         = params.snpeff_db         ?: Channel.empty()
vep_cache_version = params.vep_cache_version ?: Channel.empty()
vep_genome        = params.vep_genome        ?: Channel.empty()
vep_species       = params.vep_species       ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
cadd_indels       = params.cadd_indels       ? Channel.fromPath(params.cadd_indels).collect()       : []
cadd_indels_tbi   = params.cadd_indels_tbi   ? Channel.fromPath(params.cadd_indels_tbi).collect()   : []
cadd_wg_snvs      = params.cadd_wg_snvs      ? Channel.fromPath(params.cadd_wg_snvs).collect()      : []
cadd_wg_snvs_tbi  = params.cadd_wg_snvs_tbi  ? Channel.fromPath(params.cadd_wg_snvs_tbi).collect()  : []
pon               = params.pon               ? Channel.fromPath(params.pon).collect()               : []
snpeff_cache      = params.snpeff_cache      ? Channel.fromPath(params.snpeff_cache).collect()      : []
target_bed        = params.target_bed        ? Channel.fromPath(params.target_bed).collect()        : []
vep_cache         = params.vep_cache         ? Channel.fromPath(params.vep_cache).collect()         : []

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
umi_read_structure   = params.umi_read_structure   ? "${params.umi_read_structure} ${params.umi_read_structure}": Channel.empty()


// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules

// Create samplesheets to restart from different steps
include { MAPPING_CSV               } from '../subworkflows/local/mapping_csv'
include { MARKDUPLICATES_CSV        } from '../subworkflows/local/markduplicates_csv'
include { PREPARE_RECALIBRATION_CSV } from '../subworkflows/local/prepare_recalibration_csv'
include { RECALIBRATE_CSV           } from '../subworkflows/local/recalibrate_csv'

// Build indices if needed
include { PREPARE_GENOME            } from '../subworkflows/local/prepare_genome'

include { ALIGNMENT_TO_FASTQ } from '../subworkflows/local/bam2fastq'

// Map input reads to reference genome (+QC)
include { GATK4_MAPPING             } from '../subworkflows/nf-core/gatk4_mapping/main'

// Mark duplicates (+QC) + convert to CRAM
include { MARKDUPLICATES            } from '../subworkflows/nf-core/markduplicates'

// Create recalibration tables
include { PREPARE_RECALIBRATION     } from '../subworkflows/nf-core/prepare_recalibration'

// Create recalibrated cram files to use for variant calling
include { RECALIBRATE               } from '../subworkflows/nf-core/recalibrate'

// Variant calling on a single normal sample
include { GERMLINE_VARIANT_CALLING  } from '../subworkflows/local/germline_variant_calling'

// Variant calling on a single tumor sample
// include { TUMOR_VARIANT_CALLING     } from '../subworkflows/local/tumor_variant_calling'
// Variant calling on tumor/normal pair
// include { PAIR_VARIANT_CALLING      } from '../subworkflows/local/pair_variant_calling'

// Annotation
include { ANNOTATE                     } from '../subworkflows/local/annotate' addParams(
    annotation_cache:                  params.annotation_cache
)

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { CREATE_UMI_CONSENSUS } from '../subworkflows/nf-core/fgbio_create_umi_consensus/main'

// Config files
ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

//
// SUBWORKFLOWS
//

include { FASTQC_TRIMGALORE    } from '../subworkflows/nf-core/fastqc_trimgalore'

//
// MODULES: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main'

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
        params.tools,
        params.step)

    // Gather built indices or get them from the params
    bwa                   = params.fasta             ? params.bwa                   ? Channel.fromPath(params.bwa).collect()                   : PREPARE_GENOME.out.bwa                   : []
    dict                  = params.fasta             ? params.dict                  ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
    fasta_fai             = params.fasta             ? params.fasta_fai             ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []
    dbsnp_tbi             = params.dbsnp             ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.empty()
    germline_resource_tbi = params.germline_resource ? params.germline_resource_tbi ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
    known_indels_tbi      = params.known_indels      ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.empty()
    pon_tbi               = params.pon               ? params.pon_tbi               ? Channel.fromPath(params.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : []
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
    // `--d fastqc` to skip fastqc
    // trim only with `--trim_fastq`
    // additional options to be set up

    if (params.step == 'mapping') {

        FASTQC_TRIMGALORE(
            input_sample,
            ('fastqc' in params.skip_qc),
            !(params.trim_fastq))

        // Get reads after optional trimming (+QC)
        reads_input = FASTQC_TRIMGALORE.out.reads

        // Get all qc reports for MultiQC
        qc_reports = qc_reports.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
        qc_reports = qc_reports.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))
        qc_reports = qc_reports.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))

        // Get versions from all software used
        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

        //Since read need additional mapping afterwards, I would argue for haveing the process here
        if(params.umi_read_structure){
            CREATE_UMI_CONSENSUS(reads_input, fasta, bwa, umi_read_structure, params.group_by_umi_strategy, params.aligner)
            ALIGNMENT_TO_FASTQ( CREATE_UMI_CONSENSUS.out.consensusbam, [] )
            ALIGNMENT_TO_FASTQ.out.reads.set{reads_input}
        }

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        GATK4_MAPPING(
            params.aligner,
            bwa,
            fasta,
            fasta_fai,
            reads_input,
            params.skip_markduplicates,
            save_bam_mapped)

        // Get mapped reads (BAM) with and without index
        bam_mapped  = GATK4_MAPPING.out.bam
        bam_indexed = GATK4_MAPPING.out.bam_indexed

        // Create CSV to restart from this step
        MAPPING_CSV(bam_indexed, save_bam_mapped, params.skip_markduplicates)

        // Get versions from all software used
        ch_versions = ch_versions.mix(GATK4_MAPPING.out.versions)
    }

    if (params.step == 'preparerecalibration') {
        if (params.skip_markduplicates) bam_indexed         = input_sample
        else                            cram_markduplicates = input_sample
    }

    if (params.step in ['mapping', 'preparerecalibration']) {
        // STEP 2: Mark duplicates (+QC) + convert to CRAM
        MARKDUPLICATES(
            bam_mapped,
            bam_indexed,
            ('markduplicates' in params.use_gatk_spark),
            !('markduplicates' in params.skip_qc),
            dict,
            fasta,
            fasta_fai,
            params.skip_markduplicates,
            ('bamqc' in params.skip_qc),
            ('samtools' in params.skip_qc),
            target_bed)

        cram_markduplicates = MARKDUPLICATES.out.cram

        // Create CSV to restart from this step
        MARKDUPLICATES_CSV(cram_markduplicates)

        qc_reports = qc_reports.mix(MARKDUPLICATES.out.qc.collect{it[1]}.ifEmpty([]))

        // STEP 3: Create recalibration tables
        if(!params.skip_bqsr){
            PREPARE_RECALIBRATION(
                cram_markduplicates,
                ('bqsr' in params.use_gatk_spark),
                dict,
                fasta,
                fasta_fai,
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

    if (params.step == 'recalibrate') bam_applybqsr = input_sample

    if (params.step in ['mapping', 'preparerecalibration', 'recalibrate']) {
        if(!params.skip_bqsr){
            // STEP 4: RECALIBRATING
            RECALIBRATE(
                ('bqsr' in params.use_gatk_spark),
                ('bamqc' in params.skip_qc),
                ('samtools' in params.skip_qc),
                cram_applybqsr,
                dict,
                fasta,
                fasta_fai,
                intervals,
                num_intervals,
                target_bed)

            cram_recalibrated    = RECALIBRATE.out.cram
            cram_recalibrated_qc = RECALIBRATE.out.qc

            RECALIBRATE_CSV(cram_recalibrated)

            qc_reports = qc_reports.mix(cram_recalibrated_qc.collect{it[1]}.ifEmpty([]))
            cram_variant_calling = cram_recalibrated

        }else{
            cram_variant_calling = cram_markduplicates
        }

    }

    if (params.step in 'variantcalling') cram_variant_calling = input_sample

    if (params.tools) {
        vcf_to_annotate = Channel.empty()
        if (params.step in 'annotate') cram_variant_calling = Channel.empty()

        // GERMLINE VARIANT CALLING
        GERMLINE_VARIANT_CALLING(
            params.tools,
            cram_variant_calling,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
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
        //     fasta,
        //     fasta_fai,
        //     intervals,
        //     target_bed,
        //     target_bed_gz_tbi)

        // PAIR VARIANT CALLING
        // PAIR_VARIANT_CALLING(
        //     params.tools,
        //     cram_variant_calling,
        //     dbsnp,
        //     dbsnp_tbi,
        //     dict,
        //     fasta_fai,
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
        if (params.step == 'annotate') vcf_to_annotate = input_sample

        if (params.tools.contains('merge') || params.tools.contains('snpeff') || params.tools.contains('vep')) {

            ANNOTATE(
                vcf_to_annotate,
                params.tools,
                snpeff_db,
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache)
        }
    }

    ch_version_yaml = Channel.empty()
    if (!('versions' in params.skip_qc)) {
        CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
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
    if (!('multiqc' in params.skip_qc)) {
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
