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
    //params.target_bed,
    params.vep_cache
]

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Check mandatory parameters
if (params.input) csv_file = file(params.input)
else {
    log.warn "No samplesheet specified, attempting to restart from csv files present in ${params.outdir}"
    switch (params.step) {
        case 'mapping': exit 1, "Can't start with step $params.step without samplesheet"
        case 'prepare_recalibration': csv_file = file("${params.outdir}/preprocessing/csv/markduplicates_no_table.csv", checkIfExists: true); break
        case 'recalibrate':          csv_file = file("${params.outdir}/preprocessing/csv/markduplicates.csv",          checkIfExists: true); break
        case 'variant_calling':       csv_file = file("${params.outdir}/preprocessing/csv/recalibrated.csv",            checkIfExists: true); break
        // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
        case 'annotate':             csv_file = file("${params.outdir}/variant_calling/csv/recalibrated.csv",          checkIfExists: true); break
        default: exit 1, "Unknown step $params.step"
    }
}

input_sample = extract_csv(csv_file)

if (params.wes) {
    if (!params.intervals.endsWith("bed")) {
        exit 1, "Target file must be in BED format"
    }
} else {
    if (!params.intervals.endsWith("bed") && !params.intervals.endsWith("interval_list")) {
        exit 1, "Interval file must end with .bed or .interval_list"
    }
}

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

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
chr_dir           = params.chr_dir           ? Channel.fromPath(params.chr_dir).collect()           : []
chr_length        = params.chr_length        ? Channel.fromPath(params.chr_length).collect()        : []
dbsnp             = params.dbsnp             ? Channel.fromPath(params.dbsnp).collect()             : Channel.empty()
fasta             = params.fasta             ? Channel.fromPath(params.fasta).collect()             : Channel.empty()
fasta_fai         = params.fasta_fai         ? Channel.fromPath(params.fasta_fai).collect()         : Channel.empty()
germline_resource = params.germline_resource ? Channel.fromPath(params.germline_resource).collect() : Channel.empty()
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
pon               = params.pon               ? Channel.fromPath(params.pon).collect()               : Channel.empty()
snpeff_cache      = params.snpeff_cache      ? Channel.fromPath(params.snpeff_cache).collect()      : []
//target_bed        = params.target_bed        ? Channel.fromPath(params.target_bed).collect()        : []
vep_cache         = params.vep_cache         ? Channel.fromPath(params.vep_cache).collect()         : []

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
umi_read_structure   = params.umi_read_structure   ? "${params.umi_read_structure} ${params.umi_read_structure}": Channel.empty()


// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules

// Create samplesheets to restart from different steps
include { MAPPING_CSV                } from '../subworkflows/local/mapping_csv'
include { MARKDUPLICATES_CSV         } from '../subworkflows/local/markduplicates_csv'
include { PREPARE_RECALIBRATION_CSV  } from '../subworkflows/local/prepare_recalibration_csv'
include { RECALIBRATE_CSV            } from '../subworkflows/local/recalibrate_csv'

// Build indices if needed
include { PREPARE_GENOME             } from '../subworkflows/local/prepare_genome'

// Build intervals if needed
include { PREPARE_INTERVALS          } from '../subworkflows/local/prepare_intervals'

// Convert BAM files to FASTQ files
include { ALIGNMENT_TO_FASTQ         } from '../subworkflows/local/bam2fastq'

// Split FASTQ files
include { SPLIT_FASTQ                } from '../subworkflows/local/split_fastq'

// Map input reads to reference genome
include { GATK4_MAPPING              } from '../subworkflows/nf-core/gatk4/mapping/main'

// Merge and index BAM files (optional)
include { MERGE_INDEX_BAM            } from '../subworkflows/nf-core/merge_index_bam'

// Mark Duplicates (+QC)
include { MARKDUPLICATES             } from '../subworkflows/nf-core/gatk4/markduplicates/main'

// Mark Duplicates_SPARK (+QC)
include { MARKDUPLICATES_SPARK       } from '../subworkflows/nf-core/gatk4/markduplicates_spark/main'

// Convert to CRAM (+QC)
include { BAM_TO_CRAM                } from '../subworkflows/nf-core/bam_to_cram'

// CRAM QC
include { CRAM_QC                    } from '../subworkflows/nf-core/cram_qc'

// Create recalibration tables
include { PREPARE_RECALIBRATION      } from '../subworkflows/nf-core/gatk4/prepare_recalibration/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { RECALIBRATE                } from '../subworkflows/nf-core/gatk4/recalibrate/main'

// Variant calling on a single normal sample
include { GERMLINE_VARIANT_CALLING   } from '../subworkflows/local/germline_variant_calling'

// Variant calling on a single tumor sample
include { TUMOR_ONLY_VARIANT_CALLING } from '../subworkflows/local/tumor_variant_calling'

// Variant calling on tumor/normal pair
include { PAIR_VARIANT_CALLING       } from '../subworkflows/local/pair_variant_calling'

// Annotation
include { ANNOTATE                   } from '../subworkflows/local/annotate' addParams(
    annotation_cache:                  params.annotation_cache
)

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

// Config files
ch_multiqc_config        = Channel.fromPath("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

//
// SUBWORKFLOWS
//

include { FASTQC_TRIMGALORE    } from '../subworkflows/nf-core/fastqc_trimgalore'

// Create umi consensus bams from fastq
include { CREATE_UMI_CONSENSUS } from '../subworkflows/nf-core/fgbio_create_umi_consensus/main'

//
// MODULES: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main'

def multiqc_report = []

workflow SAREK {

    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()

    // Build indices if needed
    PREPARE_GENOME(
        dbsnp,
        fasta,
        fasta_fai,
        germline_resource,
        known_indels,
        pon)

    // Gather built indices or get them from the params
    bwa                    = params.fasta                   ? params.bwa                   ? Channel.fromPath(params.bwa).collect()                   : PREPARE_GENOME.out.bwa                   : []
    dict                   = params.fasta                   ? params.dict                  ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
    fasta_fai              = params.fasta                   ? params.fasta_fai             ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []
    dbsnp_tbi              = params.dbsnp                   ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.empty()
    germline_resource_tbi  = params.germline_resource       ? params.germline_resource_tbi ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
    known_indels_tbi       = params.known_indels            ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.empty()
    pon_tbi                = params.pon                     ? params.pon_tbi               ? Channel.fromPath(params.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : []
    msisensorpro_scan      = PREPARE_GENOME.out.msisensorpro_scan

    //TODO @Rike, is this working for you? Now it is, fixed a bug in prepare_genome.nf after chasing smoke for a while
    // known_sites is made by grouping both the dbsnp and the known indels ressources
    // Which can either or both be optional
    // Actually BQSR has been throwing erros if no sides were provided so it must be at least one
    known_sites     = dbsnp.concat(known_indels).collect()
    known_sites_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals_bed_combined        = (params.intervals && params.wes) ? Channel.fromPath(params.intervals).collect() : []
    intervals                     = PREPARE_INTERVALS.out.intervals_bed                             // multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi          = PREPARE_INTERVALS.out.intervals_bed_gz_tbi                      // multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather
    intervals_bed_combined_gz_tbi = PREPARE_INTERVALS.out.intervals_combined_bed_gz_tbi.collect()   // one file containing all intervals interval.bed.gz/.tbi file
    intervals_bed_combined_gz     = intervals_bed_combined_gz_tbi.map{ bed, tbi -> [bed]}.collect() // one file containing all intervals interval.bed.gz file
    intervals_for_preprocessing   = (!params.wes || params.no_intervals) ? [] : PREPARE_INTERVALS.out.intervals_bed //TODO: intervals also with WGS data? Probably need a parameter if WGS for deepvariant tool, that would allow to check here too

    num_intervals = 0
    intervals.count().map{ num_intervals = it }

    // Gather used softwares versions
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

        if (params.is_bam_input) {
            ALIGNMENT_TO_FASTQ(input_sample, [])
            input_sample_converted = ALIGNMENT_TO_FASTQ.out.reads
            ch_versions = ch_versions.mix(ALIGNMENT_TO_FASTQ.out.versions)
        } else input_sample_converted = input_sample

        FASTQC_TRIMGALORE(input_sample_converted)

        // Get reads after optional trimming (+QC)
        reads_input = FASTQC_TRIMGALORE.out.reads

        // Gather QC reports
        ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
        ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))

        //Since read need additional mapping afterwards, I would argue for having the process here
        if (params.umi_read_structure) {
            CREATE_UMI_CONSENSUS(reads_input, fasta, bwa, umi_read_structure, params.group_by_umi_strategy, params.aligner)
            ALIGNMENT_TO_FASTQ(CREATE_UMI_CONSENSUS.out.consensusbam, [])
            reads_input = ALIGNMENT_TO_FASTQ.out.reads

            // Gather used softwares versions
            ch_versions = ch_versions.mix(CREATE_UMI_CONSENSUS.out.versions)
            ch_versions = ch_versions.mix(ALIGNMENT_TO_FASTQ.out.versions)
        }

        // OPTIONNAL SPLIT OF FASTQ FILES WITH SEQKIT_SPLIT2
        SPLIT_FASTQ(reads_input)

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        GATK4_MAPPING(
            SPLIT_FASTQ.out.reads,
            bwa,
            fasta,
            fasta_fai)

        // Get mapped reads (BAM) with and without index
        // without index: always contains mapped_bams, only used if duplicate marking is done
        // with Index: Duplicate marking is skipped and/or bams are saved, else empty Channel
        bam_mapped  = GATK4_MAPPING.out.bam

        if (params.save_bam_mapped || (params.skip_tools && params.skip_tools.contains('markduplicates'))) {
            MERGE_INDEX_BAM(bam_mapped)
            bam_indexed = MERGE_INDEX_BAM.out.bam_bai

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MERGE_INDEX_BAM.out.versions)
        }

        // Create CSV to restart from this step
        // TODO: How is this handeled if not save_bam_mapped is set (no index should be present)
        //MAPPING_CSV(bam_indexed, save_bam_mapped, params.skip_tools.contains('markduplicates'))

        // Gather used softwares versions
        ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
        ch_versions = ch_versions.mix(SPLIT_FASTQ.out.versions)
        ch_versions = ch_versions.mix(GATK4_MAPPING.out.versions)
    }

    // Comment out till we get the tests to pass
    if (params.step == 'prepare_recalibration') bam_mapped = input_sample.map{ meta, bam, bai -> [meta, bam] }

    if (params.step in ['mapping', 'prepare_recalibration']) {
    cram_markduplicates_no_spark = Channel.empty()
    cram_markduplicates_spark    = Channel.empty()
    cram_no_markduplicates       = Channel.empty()

        // STEP 2: markduplicates (+QC) + convert to CRAM

        if (params.skip_tools && params.skip_tools.contains('markduplicates')) {
            if (params.step == 'prepare_recalibration') bam_indexed = input_sample

            BAM_TO_CRAM(
                bam_indexed,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_no_markduplicates = BAM_TO_CRAM.out.cram

            // Gather used softwares versions
            ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)

            // Gather QC reports
            ch_reports  = ch_reports.mix(BAM_TO_CRAM.out.qc)
        } else if (params.use_gatk_spark && params.use_gatk_spark.contains('markduplicates')) {
            MARKDUPLICATES_SPARK(
                bam_mapped,
                dict,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)
            cram_markduplicates_spark = MARKDUPLICATES_SPARK.out.cram

            // Gather QC reports
            ch_reports  = ch_reports.mix(MARKDUPLICATES_SPARK.out.qc.collect{it[1]}.ifEmpty([]))

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MARKDUPLICATES_SPARK.out.versions)
        } else {
            MARKDUPLICATES(
                bam_mapped,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)
            cram_markduplicates_no_spark = MARKDUPLICATES.out.cram

            // Gather QC reports
            ch_reports  = ch_reports.mix(MARKDUPLICATES.out.qc.collect{it[1]}.ifEmpty([]))

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)
        }

        cram_for_prepare_recalibration = Channel.empty().mix(
            cram_markduplicates_no_spark,
            cram_markduplicates_spark,
            cram_no_markduplicates)

        // Run Samtools stats on CRAM
        CRAM_QC(cram_for_prepare_recalibration, fasta)

        // Gather QC reports
        ch_reports  = ch_reports.mix(CRAM_QC.out.qc.collect{it[1]}.ifEmpty([]))

        // Gather used softwares versions
        ch_versions = ch_versions.mix(CRAM_QC.out.versions)

        // Create CSV to restart from this step
        MARKDUPLICATES_CSV(cram_for_prepare_recalibration)

        // STEP 3: Create recalibration tables
        if (!(params.skip_tools && params.skip_tools.contains('baserecalibrator'))) {
            PREPARE_RECALIBRATION(
                cram_for_prepare_recalibration,
                dict,
                fasta,
                fasta_fai,
                intervals,
                num_intervals,
                known_sites,
                known_sites_tbi,
                params.no_intervals)

            PREPARE_RECALIBRATION_CSV(PREPARE_RECALIBRATION.out.table_bqsr)

            cram_applybqsr = cram_for_prepare_recalibration.join(PREPARE_RECALIBRATION.out.table_bqsr)

            // Gather used softwares versions
            ch_versions = ch_versions.mix(PREPARE_RECALIBRATION.out.versions)
        }
    }

    if (params.step == 'recalibrate') cram_applybqsr = input_sample

    if (params.step in ['mapping', 'prepare_recalibration', 'recalibrate']) {

        if (!(params.skip_tools && params.skip_tools.contains('baserecalibrator'))) {
            // STEP 4: RECALIBRATING
            RECALIBRATE(
                cram_applybqsr,
                dict,
                fasta,
                fasta_fai,
                intervals,
                num_intervals,
                params.no_intervals,
                intervals_for_preprocessing)

            RECALIBRATE_CSV(RECALIBRATE.out.cram)

            cram_variant_calling = RECALIBRATE.out.cram

            // Gather QC reports
            ch_reports  = ch_reports.mix(RECALIBRATE.out.qc.collect{it[1]}.ifEmpty([]))

            // Gather used softwares versions
            ch_versions = ch_versions.mix(RECALIBRATE.out.versions)

        } else cram_variant_calling = cram_for_prepare_recalibration

    }

    if (params.step in 'variant_calling') cram_variant_calling = input_sample

    if (params.tools) {

        if (params.step in 'annotate') cram_variant_calling = Channel.empty()

        //
        // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
        //
        cram_variant_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }.set{cram_variant_calling_status}

        // All Germline samples
        cram_variant_calling_normal_to_cross = cram_variant_calling_status.normal.map{ meta, cram, crai -> [meta.patient, meta, cram, crai] }

        // All tumor samples
        cram_variant_calling_pair_to_cross = cram_variant_calling_status.tumor.map{ meta, cram, crai -> [meta.patient, meta, cram, crai] }

        // Tumor only samples
        // 1. Group together all tumor samples by patient ID [patient1, [meta1, meta2], [cram1,crai1, cram2, crai2]]

        // Downside: this only works by waiting for all tumor samples to finish preprocessing, since no group size is provided
        cram_variant_calling_tumor_grouped = cram_variant_calling_pair_to_cross.groupTuple()

        // 2. Join with normal samples, in each channel there is one key per patient now. Patients without matched normal end up with: [patient1, [meta1, meta2], [cram1,crai1, cram2, crai2], null]
        cram_variant_calling_tumor_joined = cram_variant_calling_tumor_grouped.join(cram_variant_calling_normal_to_cross, remainder: true)

        // 3. Filter out entries with last entry null
        cram_variant_calling_tumor_filtered = cram_variant_calling_tumor_joined.filter{ it ->  !(it.last()) }

        // 4. Transpose [patient1, [meta1, meta2], [cram1,crai1, cram2, crai2]] back to [patient1, meta1, [cram1,crai1], null] [patient1, meta2, [cram2,crai2], null]
        // and remove patient ID field & null value for further processing [meta1, [cram1,crai1]] [meta2, [cram2,crai2]]
        cram_variant_calling_tumor_only = cram_variant_calling_tumor_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

        // Tumor - normal pairs
        // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
        cram_variant_calling_pair = cram_variant_calling_normal_to_cross.cross(cram_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]
                meta.patient    = normal[0]
                meta.normal_id  = normal[1].sample
                meta.tumor_id   = tumor[1].sample
                meta.gender     = normal[1].gender
                meta.id         = "${meta.tumor_id}_vs_${meta.normal_id}".toString()

                [meta, normal[2], normal[3], tumor[2], tumor[3]]
            }

        // GERMLINE VARIANT CALLING
        GERMLINE_VARIANT_CALLING(
            cram_variant_calling_status.normal,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            intervals,
            intervals_bed_gz_tbi,
            intervals_bed_combined_gz_tbi,
            intervals_bed_combined_gz,
            num_intervals)
            // params.joint_germline)

        // TUMOR ONLY VARIANT CALLING
        TUMOR_ONLY_VARIANT_CALLING(
            params.tools,
            cram_variant_calling_tumor_only,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            intervals,
            intervals_bed_gz_tbi,
            intervals_bed_combined_gz_tbi,
            intervals_bed_combined_gz,
            num_intervals,
            params.no_intervals,
            germline_resource,
            germline_resource_tbi,
            pon,
            pon_tbi
        )

        // PAIR VARIANT CALLING
        PAIR_VARIANT_CALLING(
            params.tools,
            cram_variant_calling_pair,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            intervals,
            intervals_bed_gz_tbi,
            intervals_bed_combined_gz_tbi,
            intervals_bed_combined_gz,
            num_intervals,
            params.no_intervals,
            msisensorpro_scan,
            germline_resource,
            germline_resource_tbi,
            pon,
            pon_tbi)

        // Gather vcf files for annotation
        vcf_to_annotate = Channel.empty()
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.deepvariant_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.freebayes_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.haplotypecaller_gvcf)
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.manta_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.strelka_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(TUMOR_ONLY_VARIANT_CALLING.out.freebayes_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(TUMOR_ONLY_VARIANT_CALLING.out.mutect2_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(TUMOR_ONLY_VARIANT_CALLING.out.manta_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(TUMOR_ONLY_VARIANT_CALLING.out.strelka_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(PAIR_VARIANT_CALLING.out.mutect2_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(PAIR_VARIANT_CALLING.out.manta_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(PAIR_VARIANT_CALLING.out.strelka_vcf)

        // Gather used softwares versions
        ch_versions = ch_versions.mix(GERMLINE_VARIANT_CALLING.out.versions)
        ch_versions = ch_versions.mix(PAIR_VARIANT_CALLING.out.versions)
        ch_versions = ch_versions.mix(TUMOR_ONLY_VARIANT_CALLING.out.versions)

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

            // Gather used softwares versions
            ch_versions = ch_versions.mix(ANNOTATE.out.versions)
        }
    }

    ch_version_yaml = Channel.empty()
        if (!(params.skip_tools && params.skip_tools.contains('versions'))) {
        CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
        ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }

    // workflow_summary    = WorkflowSarek.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_config)
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_version_yaml)
    // ch_multiqc_files = ch_multiqc_files.mix(ch_reports)

    // multiqc_report = Channel.empty()
    // if (!(params.skip_tools && params.skip_tools.contains('multiqc'))) {
    //     MULTIQC(ch_multiqc_files.collect())
    //     multiqc_report = MULTIQC.out.report.toList()

    // }
}

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

// Function to extract information (meta data + file(s)) from csv file(s)

is_bam_input = false
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
        // start from BAM
        } else if (row.lane && row.bam) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:ILLUMINA\""
            meta.numLanes = numLanes.toInteger()
            meta.read_group = read_group.toString()
            is_bam_input = true
            return [meta, bam]
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
