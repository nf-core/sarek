/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSarek.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [
    params.ac_loci,
    params.ac_loci_gc,
    params.bwa,
    params.bwamem2,
    params.cadd_indels,
    params.cadd_indels_tbi,
    params.cadd_wg_snvs,
    params.cadd_wg_snvs_tbi,
    params.chr_dir,
    params.dbsnp,
    params.dbsnp_tbi,
    params.dict,
    params.dragmap,
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
        //case 'markduplicates': csv_file = file("${params.outdir}/preprocessing/csv/markduplicates_no_table.csv", checkIfExists: true); break
        case 'prepare_recalibration': csv_file = file("${params.outdir}/preprocessing/csv/markduplicates_no_table.csv", checkIfExists: true); break
        case 'recalibrate':           csv_file = file("${params.outdir}/preprocessing/csv/markduplicates.csv",          checkIfExists: true); break
        case 'variant_calling':       csv_file = file("${params.outdir}/preprocessing/csv/recalibrated.csv",            checkIfExists: true); break
        // case 'controlfreec':         csv_file = file("${params.outdir}/variant_calling/csv/control-freec_mpileup.csv", checkIfExists: true); break
        case 'annotate':              csv_file = file("${params.outdir}/variant_calling/csv/recalibrated.csv",          checkIfExists: true); break
        default: exit 1, "Unknown step $params.step"
    }
}

ch_input_sample = extract_csv(csv_file)

if (params.wes) {
    if (params.intervals && !params.intervals.endsWith("bed")) exit 1, "Target file must be in BED format"
} else {
    if (params.intervals && !params.intervals.endsWith("bed") && !params.intervals.endsWith("interval_list")) exit 1, "Interval file must end with .bed or .interval_list"
}

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[params.genome]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
chr_dir            = params.chr_dir            ? Channel.fromPath(params.chr_dir).collect()                  : []
dbsnp              = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()                    : Channel.empty()
fasta              = params.fasta              ? Channel.fromPath(params.fasta).collect()                    : Channel.empty()
fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()                : Channel.empty()
germline_resource  = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect()        : Channel.empty()
known_indels       = params.known_indels       ? Channel.fromPath(params.known_indels).collect()             : Channel.empty()
loci               = params.ac_loci            ? Channel.fromPath(params.ac_loci).collect()                  : []
loci_gc            = params.ac_loci_gc         ? Channel.fromPath(params.ac_loci_gc).collect()               : []
mappability        = params.mappability        ? Channel.fromPath(params.mappability).collect()              : []

// Initialize value channels based on params, defined in the params.genomes[params.genome] scope
snpeff_db          = params.snpeff_db          ?: Channel.empty()
vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
vep_genome         = params.vep_genome         ?: Channel.empty()
vep_species        = params.vep_species        ?: Channel.empty()

// Initialize files channels based on params, not defined within the params.genomes[params.genome] scope
cadd_indels        = params.cadd_indels        ? Channel.fromPath(params.cadd_indels).collect()              : []
cadd_indels_tbi    = params.cadd_indels_tbi    ? Channel.fromPath(params.cadd_indels_tbi).collect()          : []
cadd_wg_snvs       = params.cadd_wg_snvs       ? Channel.fromPath(params.cadd_wg_snvs).collect()             : []
cadd_wg_snvs_tbi   = params.cadd_wg_snvs_tbi   ? Channel.fromPath(params.cadd_wg_snvs_tbi).collect()         : []
pon                = params.pon                ? Channel.fromPath(params.pon).collect()                      : Channel.empty()
snpeff_cache       = params.snpeff_cache       ? Channel.fromPath(params.snpeff_cache).collect()             : []
//target_bed         = params.target_bed         ? Channel.fromPath(params.target_bed).collect()               : []
vep_cache          = params.vep_cache          ? Channel.fromPath(params.vep_cache).collect()                : []

// Initialize value channels based on params, not defined within the params.genomes[params.genome] scope
umi_read_structure = params.umi_read_structure ? "${params.umi_read_structure} ${params.umi_read_structure}" : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Create samplesheets to restart from different steps
include { MAPPING_CSV                                    } from '../subworkflows/local/mapping_csv'
include { MARKDUPLICATES_CSV                             } from '../subworkflows/local/markduplicates_csv'
include { PREPARE_RECALIBRATION_CSV                      } from '../subworkflows/local/prepare_recalibration_csv'
include { RECALIBRATE_CSV                                } from '../subworkflows/local/recalibrate_csv'

// Build indices if needed
include { PREPARE_GENOME                                 } from '../subworkflows/local/prepare_genome'

// Build intervals if needed
include { PREPARE_INTERVALS                              } from '../subworkflows/local/prepare_intervals'

// Convert BAM files to FASTQ files
include { ALIGNMENT_TO_FASTQ as ALIGNMENT_TO_FASTQ_INPUT } from '../subworkflows/nf-core/alignment_to_fastq'
include { ALIGNMENT_TO_FASTQ as ALIGNMENT_TO_FASTQ_UMI   } from '../subworkflows/nf-core/alignment_to_fastq'

// Split FASTQ files
include { SPLIT_FASTQ                                    } from '../subworkflows/local/split_fastq'

// Run FASTQC
include { RUN_FASTQC                                     } from '../subworkflows/nf-core/run_fastqc'

// Run TRIMGALORE
include { RUN_TRIMGALORE                                 } from '../subworkflows/nf-core/run_trimgalore'

// Create umi consensus bams from fastq
include { CREATE_UMI_CONSENSUS                           } from '../subworkflows/nf-core/fgbio_create_umi_consensus/main'

// Map input reads to reference genome
include { GATK4_MAPPING                                  } from '../subworkflows/nf-core/gatk4/mapping/main'

// Merge and index BAM files (optional)
include { MERGE_INDEX_BAM                                } from '../subworkflows/nf-core/merge_index_bam'

include {SAMTOOLS_VIEW as SAMTOOLS_CRAMTOBAM } from '../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_BAMTOCRAM    } from '../modules/nf-core/modules/samtools/bamtocram/main'

// Mark Duplicates (+QC)
include { MARKDUPLICATES                                 } from '../subworkflows/nf-core/gatk4/markduplicates/main'

// Mark Duplicates SPARK (+QC)
include { MARKDUPLICATES_SPARK                           } from '../subworkflows/nf-core/gatk4/markduplicates_spark/main'

// Convert to CRAM (+QC)
include { BAM_TO_CRAM                                    } from '../subworkflows/nf-core/bam_to_cram'
include { MAPPING_CRAM_QC                                } from '../subworkflows/nf-core/mapping_cram_qc'
// QC on CRAM
include { SAMTOOLS_STATS as SAMTOOLS_STATS_CRAM          } from '../modules/nf-core/modules/samtools/stats/main'
include { CRAM_QC                                        } from '../subworkflows/nf-core/cram_qc'

// Create recalibration tables
include { PREPARE_RECALIBRATION                          } from '../subworkflows/nf-core/gatk4/prepare_recalibration/main'

// Create recalibration tables SPARK
include { PREPARE_RECALIBRATION_SPARK                    } from '../subworkflows/nf-core/gatk4/prepare_recalibration_spark/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { RECALIBRATE                                    } from '../subworkflows/nf-core/gatk4/recalibrate/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { RECALIBRATE_SPARK                              } from '../subworkflows/nf-core/gatk4/recalibrate_spark/main'

// Variant calling on a single normal sample
include { GERMLINE_VARIANT_CALLING                       } from '../subworkflows/local/germline_variant_calling'

// Variant calling on a single tumor sample
include { TUMOR_ONLY_VARIANT_CALLING                     } from '../subworkflows/local/tumor_variant_calling'

// Variant calling on tumor/normal pair
include { PAIR_VARIANT_CALLING                           } from '../subworkflows/local/pair_variant_calling'

include { VCF_QC                                         } from '../subworkflows/nf-core/vcf_qc'

// Annotation
include { ANNOTATE                                       } from '../subworkflows/local/annotate'

// REPORTING VERSIONS OF SOFTWARE USED
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

// MULTIQC
include { MULTIQC                                        } from '../modules/nf-core/modules/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = Channel.fromPath(file("$projectDir/assets/multiqc_config.yml", checkIfExists: true))
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_sarek_logo            = Channel.fromPath(file("$projectDir/assets/nf-core-sarek_logo_light.png", checkIfExists: true))
def multiqc_report = []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAREK {

    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()

    // Build indices if needed
    PREPARE_GENOME(
        chr_dir,
        dbsnp,
        fasta,
        fasta_fai,
        germline_resource,
        known_indels,
        pon)

    // Gather built indices or get them from the params
    bwa                    = params.fasta                   ? params.bwa                   ? Channel.fromPath(params.bwa).collect()                   : PREPARE_GENOME.out.bwa                   : []
    chr_files              = PREPARE_GENOME.out.chr_files
    bwamem2                = params.fasta                   ? params.bwamem2               ? Channel.fromPath(params.bwamem2).collect()               : PREPARE_GENOME.out.bwamem2               : []
    dragmap                = params.fasta                   ? params.dragmap               ? Channel.fromPath(params.dragmap).collect()               : PREPARE_GENOME.out.hashtable             : []
    dict                   = params.fasta                   ? params.dict                  ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
    fasta_fai              = params.fasta                   ? params.fasta_fai             ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []
    dbsnp_tbi              = params.dbsnp                   ? params.dbsnp_tbi             ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.empty()
    germline_resource_tbi  = params.germline_resource       ? params.germline_resource_tbi ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
    known_indels_tbi       = params.known_indels            ? params.known_indels_tbi      ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.empty()
    pon_tbi                = params.pon                     ? params.pon_tbi               ? Channel.fromPath(params.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : []
    msisensorpro_scan      = PREPARE_GENOME.out.msisensorpro_scan

    // Gather index for mapping given the chosen aligner
    ch_map_index = params.aligner == "bwa-mem" ? bwa :
        params.aligner == "bwa-mem2" ? bwamem2 :
        dragmap

    // known_sites is made by grouping both the dbsnp and the known indels ressources
    // Which can either or both be optional
    // Actually BQSR has been throwing erros if no sides were provided so it must be at least one
    known_sites     = dbsnp.concat(known_indels).collect()
    known_sites_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals_bed_combined        = (params.intervals && params.wes) ? Channel.fromPath(params.intervals).collect() : []
    intervals_bed_combined_gz_tbi = PREPARE_INTERVALS.out.intervals_combined_bed_gz_tbi.collect()   // one file containing all intervals interval.bed.gz/.tbi file
    intervals_bed_combined_gz     = intervals_bed_combined_gz_tbi.map{ bed, tbi -> [bed]}.collect() // one file containing all intervals interval.bed.gz file

    intervals                     = PREPARE_INTERVALS.out.intervals_bed        // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi          = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [interval_bed, tbi, num_intervals] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    //TODO: intervals also with WGS data? Probably need a parameter if WGS for deepvariant tool, that would allow to check here too
    intervals_for_preprocessing   = (params.wes && !params.no_intervals) ? intervals_bed_combined : []

    // Gather used softwares versions
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    // PREPROCESSING

    if (params.step == 'mapping') {

        // Figure out if input is bam or fastq
        ch_input_sample.branch{
            bam:   it[0].data_type == "bam"
            fastq: it[0].data_type == "fastq"
        }.set{ch_input_sample_type}

        // convert any bam input to fastq
        ALIGNMENT_TO_FASTQ_INPUT(ch_input_sample_type.bam, [])

        // gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        ch_input_fastq = ch_input_sample_type.fastq.mix(ALIGNMENT_TO_FASTQ_INPUT.out.reads)

        // STEP 0: QC & TRIM
        // `--skip_tools fastqc` to skip fastqc
        // trim only with `--trim_fastq`
        // additional options to be set up

        // QC
        if (!(params.skip_tools && params.skip_tools.contains('fastqc'))) {
            RUN_FASTQC(ch_input_fastq)

            ch_reports  = ch_reports.mix(RUN_FASTQC.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
            ch_versions = ch_versions.mix(RUN_FASTQC.out.versions)
        }

        // Trimming
        if (params.trim_fastq) {
            RUN_TRIMGALORE(ch_input_fastq)

            ch_reads = RUN_TRIMGALORE.out.reads

            ch_reports  = ch_reports.mix(RUN_TRIMGALORE.out.trim_html.collect{it[1]}.ifEmpty([]))
            ch_reports  = ch_reports.mix(RUN_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))
            ch_reports  = ch_reports.mix(RUN_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))

            ch_versions = ch_versions.mix(RUN_TRIMGALORE.out.versions)
        } else {
            ch_reads = ch_input_fastq
        }

        // UMI consensus calling
        if (params.umi_read_structure) {
            CREATE_UMI_CONSENSUS(ch_reads,
                fasta,
                ch_map_index,
                umi_read_structure,
                params.group_by_umi_strategy)

            // convert back to fastq for further preprocessing
            ALIGNMENT_TO_FASTQ_UMI(CREATE_UMI_CONSENSUS.out.consensusbam, [])

            ch_input_sample_to_split = ALIGNMENT_TO_FASTQ_UMI.out.reads

            // Gather used softwares versions
            ch_versions = ch_versions.mix(ALIGNMENT_TO_FASTQ_UMI.out.versions)
            ch_versions = ch_versions.mix(CREATE_UMI_CONSENSUS.out.versions)
        } else {
            ch_input_sample_to_split = ch_reads
        }

        // SPLIT OF FASTQ FILES WITH SEQKIT_SPLIT2
        if (params.split_fastq > 1) {
            SPLIT_FASTQ(ch_input_sample_to_split)

            ch_reads_to_map = SPLIT_FASTQ.out.reads

            // Gather used softwares versions
            ch_versions = ch_versions.mix(SPLIT_FASTQ.out.versions)
        } else {
            ch_reads_to_map = ch_input_sample_to_split
        }

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        // reads will be sorted

        ch_reads_to_map = ch_reads_to_map.map{ meta, reads ->
            new_meta = meta.clone()

            // update ID when no multiple lanes or splitted fastqs
            new_meta.id = meta.size * meta.numLanes == 1 ? meta.sample : meta.id

            [new_meta, reads]
        }

        GATK4_MAPPING(ch_reads_to_map, ch_map_index, true)

        // Grouping the bams from the same samples not to stall the workflow
        ch_bam_mapped = GATK4_MAPPING.out.bam.map{ meta, bam ->
            new_meta = meta.clone()

            numLanes = meta.numLanes ?: 1
            size     = meta.size     ?: 1

            // remove no longer necessary fields
            new_meta.remove('read_group') // Now in the BAM header
            new_meta.remove('numLanes')   // Was only needed for mapping
            new_meta.remove('size')       // Was only needed for mapping

            // update ID to be based on the sample name
            new_meta.id = meta.sample

            // update data_type
            //TODO: This is never used again as far as I see, could probably be removed
            new_meta.data_type = 'bam'

            // Use groupKey to make sure that the correct group can advance as soon as it is complete
            // and not stall the workflow until all reads from all channels are mapped
            def groupKey = groupKey(new_meta, numLanes * size)

            //Returns the values we need
            [groupKey, new_meta, bam]
        }.groupTuple(by:[0,1])
        .map{ groupKey, new_meta, bam -> [new_meta, bam] }

        // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
        // Except if and only if skipping markduplicates or saving mapped bams
        if (params.save_bam_mapped || (params.skip_tools && params.skip_tools.contains('markduplicates'))) {

            // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
            MERGE_INDEX_BAM(ch_bam_mapped)

            // Create CSV to restart from this step
            MAPPING_CSV(MERGE_INDEX_BAM.out.bam_bai)

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MERGE_INDEX_BAM.out.versions)
        }

        // Gather used softwares versions
        ch_versions = ch_versions.mix(ALIGNMENT_TO_FASTQ_INPUT.out.versions)
        ch_versions = ch_versions.mix(GATK4_MAPPING.out.versions)
    }

    if (params.step in ['mapping', 'markduplicates']) {

        //TODO Allow bam already supported, also allow cram:
        // 1. SAMTOOLS_CRAMTOBAM ( to speed up computation)
        // 2. Need fasta for cram compression (maybe just using --fasta, because this reference will be used elsewhere)
        ch_cram_markduplicates_no_spark = Channel.empty()
        ch_cram_markduplicates_spark    = Channel.empty()
        ch_cram_no_markduplicates       = Channel.empty()

        // STEP 2: markduplicates (+QC) + convert to CRAM

        // ch_bam_for_markduplicates will countain bam mapped with GATK4_MAPPING when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        // ch_bam_for_markduplicates = params.step == 'mapping'? ch_bam_mapped : ch_input_sample.map{ meta, input, index -> [meta, input] }

        ch_bam_for_markduplicates = Channel.empty()

        if(params.step == 'mapping'){

            ch_bam_for_markduplicates = ch_bam_mapped

        }else{

            ch_input_sample.map{ meta, input, index -> [meta, input, index] }.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{convert}

            if (!(params.skip_tools && params.skip_tools.contains('markduplicates'))){

                SAMTOOLS_CRAMTOBAM(convert.cram.map{ meta, input, index -> [meta, input]}, fasta)
                ch_versions = ch_versions.mix(SAMTOOLS_CRAMTOBAM.out.versions)

                ch_bam_for_markduplicates = ch_bam_for_markduplicates.mix(
                                                SAMTOOLS_CRAMTOBAM.out.bam,
                                                convert.bam)
            }
        }

        if (params.skip_tools && params.skip_tools.contains('markduplicates')) {

            // ch_bam_indexed will countain bam mapped with GATK4_MAPPING when step is mapping
            // which are then merged and indexed
            // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
            ch_bam_indexed = params.step == 'mapping' ? MERGE_INDEX_BAM.out.bam_bai : convert.bam
            ch_cram_indexed = convert.cram

            MAPPING_CRAM_QC (ch_cram_indexed,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            BAM_TO_CRAM(ch_bam_indexed,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_no_markduplicates = Channel.empty().mix(BAM_TO_CRAM.out.cram, ch_cram_indexed)

            // Gather QC reports
            ch_reports  = ch_reports.mix(BAM_TO_CRAM.out.qc.collect{it[1]}.ifEmpty([]))
            ch_reports  = ch_reports.mix(MAPPING_CRAM_QC.out.qc.collect{it[1]}.ifEmpty([]))

            // Gather used softwares versions
            ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)
            ch_versions = ch_versions.mix(MAPPING_CRAM_QC.out.versions)
        } else if (params.use_gatk_spark && params.use_gatk_spark.contains('markduplicates')) {
            MARKDUPLICATES_SPARK(ch_bam_for_markduplicates,
                dict,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)
            ch_cram_markduplicates_spark = MARKDUPLICATES_SPARK.out.cram

            // Gather QC reports
            ch_reports  = ch_reports.mix(MARKDUPLICATES_SPARK.out.qc.collect{it[1]}.ifEmpty([]))

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MARKDUPLICATES_SPARK.out.versions)
        } else {
            MARKDUPLICATES(ch_bam_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_markduplicates_no_spark = MARKDUPLICATES.out.cram

            ch_cram_markduplicates_no_spark.view()
            // Gather QC reports
            ch_reports  = ch_reports.mix(MARKDUPLICATES.out.qc.collect{it[1]}.ifEmpty([]))

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)
        }

        // ch_cram_for_prepare_recalibration contains either:
        // - crams from markduplicates
        // - crams from markduplicates_spark
        // - crams converted from bam mapped when skipping markduplicates
        ch_cram_for_prepare_recalibration = Channel.empty().mix(
            ch_cram_markduplicates_no_spark,
            ch_cram_markduplicates_spark,
            ch_cram_no_markduplicates)

        ch_cram_for_prepare_recalibration.view()
        // Run Samtools stats on CRAM
        SAMTOOLS_STATS_CRAM(ch_cram_for_prepare_recalibration, fasta)

        // Create CSV to restart from this step
        MARKDUPLICATES_CSV(ch_cram_for_prepare_recalibration)

        // Gather QC reports
        ch_reports  = ch_reports.mix(SAMTOOLS_STATS_CRAM.out.stats.collect{it[1]}.ifEmpty([]))

        // Gather used softwares versions
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_CRAM.out.versions)
    }

    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration']) {

        // Run if starting from step "prepare_recalibration"
        if(params.step == 'prepare_recalibration'){

            //Support if starting from BAM or CRAM files
            ch_input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{convert}

            //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            SAMTOOLS_BAMTOCRAM(convert.bam, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions)

            ch_cram_for_prepare_recalibration = Channel.empty().mix(SAMTOOLS_BAMTOCRAM.out.cram_crai, convert.cram)
        }

        // STEP 3: Create recalibration tables
        if (!(params.skip_tools && params.skip_tools.contains('baserecalibrator'))) {
            ch_table_bqsr_no_spark = Channel.empty()
            ch_table_bqsr_spark    = Channel.empty()

            if (params.use_gatk_spark && params.use_gatk_spark.contains('baserecalibrator')) {
            PREPARE_RECALIBRATION_SPARK(ch_cram_for_prepare_recalibration,
                dict,
                fasta,
                fasta_fai,
                intervals,
                known_sites,
                known_sites_tbi)

                ch_table_bqsr_spark = PREPARE_RECALIBRATION_SPARK.out.table_bqsr

                // Gather used softwares versions
                ch_versions = ch_versions.mix(PREPARE_RECALIBRATION_SPARK.out.versions)
            } else {

            PREPARE_RECALIBRATION(ch_cram_for_prepare_recalibration,
                dict,
                fasta,
                fasta_fai,
                intervals,
                known_sites,
                known_sites_tbi)

                ch_table_bqsr_no_spark = PREPARE_RECALIBRATION.out.table_bqsr

                // Gather used softwares versions
                ch_versions = ch_versions.mix(PREPARE_RECALIBRATION.out.versions)
            }

            // ch_table_bqsr contains either:
            // - bqsr table from baserecalibrator
            // - bqsr table from baserecalibrator_spark
            ch_table_bqsr = Channel.empty().mix(
                ch_table_bqsr_no_spark,
                ch_table_bqsr_spark)

            // Create CSV to restart from this step
            PREPARE_RECALIBRATION_CSV(ch_table_bqsr)

            ch_reports  = ch_reports.mix(ch_table_bqsr.map{ meta, table -> table})
        }
    }

    // STEP 4: RECALIBRATING
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        // Run if starting from step "prepare_recalibration"
        if(params.step == 'recalibrate'){

            //Support if starting from BAM or CRAM files
            ch_input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{convert}

            //If BAM file, split up table and mapped file to convert BAM to CRAM
            ch_bam_table = convert.bam.map{ meta, bam, bai, table -> [meta, table]}
            ch_bam_bam   = convert.bam.map{ meta, bam, bai, table -> [meta, bam, bai]}
            //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            SAMTOOLS_BAMTOCRAM(ch_bam_bam, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions)

            ch_cram_applybqsr = Channel.empty().mix(SAMTOOLS_BAMTOCRAM.out.cram_crai.join(ch_bam_table) // Join together converted cram with input tables
                                                    , convert.cram)

        }else ch_cram_applybqsr = ch_cram_for_prepare_recalibration.join(ch_table_bqsr)

        if (!(params.skip_tools && params.skip_tools.contains('baserecalibrator'))) {
            ch_cram_variant_calling_no_spark = Channel.empty()
            ch_cram_variant_calling_spark    = Channel.empty()

            if (params.use_gatk_spark && params.use_gatk_spark.contains('baserecalibrator')) {

                RECALIBRATE_SPARK(ch_cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals)

                ch_cram_variant_calling_spark = RECALIBRATE_SPARK.out.cram

                // Gather used softwares versions
                ch_versions = ch_versions.mix(RECALIBRATE_SPARK.out.versions)

            } else {

                RECALIBRATE(ch_cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals)

                ch_cram_variant_calling_no_spark = RECALIBRATE.out.cram

                // Gather used softwares versions
                ch_versions = ch_versions.mix(RECALIBRATE.out.versions)
            }
            cram_variant_calling = Channel.empty().mix(
                ch_cram_variant_calling_no_spark,
                ch_cram_variant_calling_spark)

            CRAM_QC(cram_variant_calling,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            // Create CSV to restart from this step
            RECALIBRATE_CSV(cram_variant_calling)

            // Gather QC reports
            ch_reports  = ch_reports.mix(CRAM_QC.out.qc.collect{it[1]}.ifEmpty([]))

            // Gather used softwares versions
            ch_versions = ch_versions.mix(CRAM_QC.out.versions)
        } else cram_variant_calling = ch_cram_for_prepare_recalibration

    }

    if (params.step == 'variant_calling') cram_variant_calling = ch_input_sample

    if (params.tools) {

        if (params.step == 'annotate') cram_variant_calling = Channel.empty()

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

        if(params.only_paired_variant_calling){
            // Normal only samples

            // 1. Join with tumor samples, in each channel there is one key per patient now. Patients without matched tumor end up with: [patient1, [meta1], [cram1,crai1], null] as there is only one matched normal possible
            cram_variant_calling_normal_joined = cram_variant_calling_normal_to_cross.join(cram_variant_calling_tumor_grouped, remainder: true)

            // 2. Filter out entries with last entry null
            cram_variant_calling_normal_filtered = cram_variant_calling_normal_joined.filter{ it ->  !(it.last()) }

            // 3. Remove patient ID field & null value for further processing [meta1, [cram1,crai1]] [meta2, [cram2,crai2]] (no transposing needed since only one normal per patient ID)
            cram_variant_calling_status_normal = cram_variant_calling_normal_filtered.map{ it -> [it[1], it[2], it[3]] }

        }else{
            cram_variant_calling_status_normal = cram_variant_calling_status.normal
        }

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
            cram_variant_calling_status_normal,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            intervals,
            intervals_bed_gz_tbi,
            intervals_bed_combined_gz_tbi,
            intervals_bed_combined_gz)
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
            intervals_bed_combined,
            germline_resource,
            germline_resource_tbi,
            pon,
            pon_tbi,
            chr_files,
            mappability
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
            intervals_bed_combined,
            msisensorpro_scan,
            germline_resource,
            germline_resource_tbi,
            pon,
            pon_tbi,
            chr_files,
            mappability)

        // Gather vcf files for annotation and QC
        vcf_to_annotate = Channel.empty()
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.deepvariant_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.freebayes_vcf)
        vcf_to_annotate = vcf_to_annotate.mix(GERMLINE_VARIANT_CALLING.out.haplotypecaller_vcf)
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

        //QC
        VCF_QC(vcf_to_annotate, intervals_bed_combined)

        ch_versions = ch_versions.mix(VCF_QC.out.versions)
        ch_reports  = ch_reports.mix(VCF_QC.out.bcftools_stats.collect{it[1]}.ifEmpty([]))
        ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_counts.collect{it[1]}.ifEmpty([]))
        ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_qual.collect{it[1]}.ifEmpty([]))
        ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_filter_summary.collect{it[1]}.ifEmpty([]))

        // ANNOTATE
        if (params.step == 'annotate') vcf_to_annotate = ch_input_sample

        if (params.tools.contains('merge') || params.tools.contains('snpeff') || params.tools.contains('vep')) {

            ANNOTATE(vcf_to_annotate,
                params.tools,
                snpeff_db,
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache)

            // Gather used softwares versions
            ch_versions = ch_versions.mix(ANNOTATE.out.versions)
            ch_reports  = ch_reports.mix(ANNOTATE.out.reports)

        }
    }

    ch_version_yaml = Channel.empty()
    if (!(params.skip_tools && params.skip_tools.contains('versions'))) {
        CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
        ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }

    if (!(params.skip_tools && params.skip_tools.contains('multiqc'))) {
        workflow_summary    = WorkflowSarek.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files =  Channel.empty().mix(ch_version_yaml,
                                            ch_multiqc_custom_config.collect().ifEmpty([]),
                                            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                                            ch_reports.collect(),
                                            ch_multiqc_config,
                                            ch_sarek_logo)

        MULTIQC(ch_multiqc_files.collect())
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:${params.seq_platform}\""
            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()

            println read_group

            meta.data_type  = "fastq"
            meta.test = "test"
            meta.size       = 1 // default number of splitted fastq
            return [meta, [fastq_1, fastq_2]]
        // start from BAM
        } else if (row.lane && row.bam) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:${params.seq_platform}\""
            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = "bam"
            meta.size       = 1 // default number of splitted fastq
            return [meta, bam]
        // recalibration
        } else if (row.table && row.cram) {
            meta.id   = meta.sample
            def cram  = file(row.cram,  checkIfExists: true)
            def crai  = file(row.crai,  checkIfExists: true)
            def table = file(row.table, checkIfExists: true)
            meta.data_type  = "cram"
            return [meta, cram, crai, table]
        // recalibration when skipping MarkDuplicates
        } else if (row.table && row.bam) {
            meta.id   = meta.sample
            def bam   = file(row.bam,   checkIfExists: true)
            def bai   = file(row.bai,   checkIfExists: true)
            def table = file(row.table, checkIfExists: true)
            meta.data_type  = "bam"
            return [meta, bam, bai, table]
        // prepare_recalibration or variant_calling
        } else if (row.cram) {
            meta.id = meta.sample
            def cram = file(row.cram, checkIfExists: true)
            def crai = file(row.crai, checkIfExists: true)
            meta.data_type  = "cram"
            return [meta, cram, crai]
        // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
        } else if (row.bam) {
            meta.id = meta.sample
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)
            meta.data_type  = "bam"
            return [meta, bam, bai]
        // annotation
        } else if (row.vcf) {
            meta.id = meta.sample
            def vcf = file(row.vcf, checkIfExists: true)
            meta.data_type     = "vcf"
            meta.variantcaller = ""
            return [meta, vcf]
        } else {
            log.warn "Missing or unknown field in csv file header"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
