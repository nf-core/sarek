/*
========================================================================================
    BUILDING INDICES
========================================================================================
*/

params.build_intervals_options         = [:]
params.bwa_index_options               = [:]
params.bwamem2_index_options           = [:]
params.create_intervals_bed_options    = [:]
params.gatk4_dict_options              = [:]
params.index_target_bed_options        = [:]
params.msisensorpro_scan_options       = [:]
params.samtools_faidx_options          = [:]
params.tabix_dbsnp_options             = [:]
params.tabix_germline_resource_options = [:]
params.tabix_known_indels_options      = [:]
params.tabix_pon_options               = [:]

// Initialize channels based on params or indices that were just built

include { BUILD_INTERVALS }                              from '../../modules/local/build_intervals/main'                           addParams(options: params.build_intervals_options)
include { BWA_INDEX as BWAMEM1_INDEX }                   from '../../modules/nf-core/software/bwa/index/main'                      addParams(options: params.bwa_index_options)
include { BWAMEM2_INDEX }                                from '../../modules/nf-core/software/bwamem2/index/main'                  addParams(options: params.bwamem2_index_options)
include { CREATE_INTERVALS_BED }                         from '../../modules/local/create_intervals_bed/main'                      addParams(options: params.create_intervals_bed_options)
include { GATK4_CREATESEQUENCEDICTIONARY as GATK4_DICT } from '../../modules/nf-core/software/gatk4/createsequencedictionary/main' addParams(options: params.gatk4_dict_options)
include { TABIX_TABIX as TABIX_DBSNP }                   from '../../modules/nf-core/software/tabix/tabix/main'                    addParams(options: params.tabix_dbsnp_options)
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE }       from '../../modules/nf-core/software/tabix/tabix/main'                    addParams(options: params.tabix_germline_resource_options)
include { TABIX_TABIX as TABIX_KNOWN_INDELS }            from '../../modules/nf-core/software/tabix/tabix/main'                    addParams(options: params.tabix_known_indels_options)
include { TABIX_TABIX as TABIX_PON }                     from '../../modules/nf-core/software/tabix/tabix/main'                    addParams(options: params.tabix_pon_options)
include { INDEX_TARGET_BED }                             from '../../modules/local/index_target_bed/main'                          addParams(options: params.index_target_bed_options)
include { MSISENSORPRO_SCAN }                            from '../../modules/nf-core/software/msisensorpro/scan/main'              addParams(options: params.msisensorpro_scan_options)
include { SAMTOOLS_FAIDX }                               from '../../modules/nf-core/software/samtools/faidx/main'                 addParams(options: params.samtools_faidx_options)

workflow BUILD_INDICES {
    take:
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon
        target_bed        // channel: [optionnal] target_bed

    main:

    result_bwa  = Channel.empty()
    version_bwa = Channel.empty()
    if (!(params.bwa) && 'mapping' in params.step.toLowerCase())
        if (params.aligner == "bwa-mem") (result_bwa, version_bwa) = BWAMEM1_INDEX(fasta)
        else                             (result_bwa, version_bwa) = BWAMEM2_INDEX(fasta)

    result_dict = Channel.empty()
    version_dict = Channel.empty()
    if (!(params.dict) && !('annotate' in params.step.toLowerCase()) && !('controlfreec' in params.step.toLowerCase()))
        (result_dict, version_dict) = GATK4_DICT(fasta)

    result_fai = Channel.empty()
    version_fai = Channel.empty()
    if (!(params.fasta_fai) && !('annotate' in params.step.toLowerCase()))
        (result_fai, version_fai) = SAMTOOLS_FAIDX(fasta)

    result_dbsnp_tbi = Channel.empty()
    version_dbsnp_tbi = Channel.empty()
    if (!(params.dbsnp_index) && params.dbsnp && ('mapping' in params.step.toLowerCase() || 'prepare_recalibration' in params.step.toLowerCase() || 'controlfreec' in params.tools.toString().toLowerCase() || 'haplotypecaller' in params.tools.toString().toLowerCase() || 'mutect2' in params.tools.toString().toLowerCase() || 'tnscope' in params.tools.toString().toLowerCase()))
        (result_dbsnp_tbi, version_dbsnp_tbi) = TABIX_DBSNP([[id:"${dbsnp.fileName}"], dbsnp])
    result_dbsnp_tbi = result_dbsnp_tbi.map {meta, tbi -> [tbi]}

    result_target_bed = Channel.empty()
    version_target_bed = Channel.empty()
    if ((params.target_bed) && ('manta' in params.tools.toString().toLowerCase() || 'strelka' in params.tools.toString().toLowerCase()))
        (result_target_bed, version_target_bed) = INDEX_TARGET_BED(target_bed)

    result_germline_resource_tbi = Channel.empty()
    version_germline_resource_tbi = Channel.empty()
    if (!(params.germline_resource_index) && params.germline_resource && 'mutect2' in params.tools.toString().toLowerCase())
        (result_germline_resource_tbi, version_germline_resource_tbi) = TABIX_GERMLINE_RESOURCE([[id:"${germline_resource.fileName}"], germline_resource])

    result_known_indels_tbi = Channel.empty()
    version_known_indels_tbi = Channel.empty()
    if (!(params.known_indels_index) && params.known_indels && ('mapping' in params.step.toLowerCase() || 'prepare_recalibration' in params.step.toLowerCase()))
        (result_known_indels_tbi, version_known_indels_tbi) = TABIX_KNOWN_INDELS([[id:"${known_indels.fileName}"], known_indels])
    result_known_indels_tbi = result_known_indels_tbi.map {meta, tbi -> [tbi]}

    result_msisensorpro_scan = Channel.empty()
    version_msisensorpro_scan = Channel.empty()
    if ('msisensorpro' in params.tools.toString().toLowerCase())
        (result_msisensorpro_scan, version_msisensorpro_scan) = MSISENSORPRO_SCAN(fasta)

    result_pon_tbi = Channel.empty()
    version_pon_tbi = Channel.empty()
    if (!(params.pon_index) && params.pon && ('tnscope' in params.tools.toString().toLowerCase() || 'mutect2' in params.tools.toString().toLowerCase()))
        (result_pon_tbi, version_pon_tbi) = TABIX_PON([[id:"${pon.fileName}"], pon])

    result_intervals = Channel.empty()
    if (params.no_intervals) {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        result_intervals = Channel.from(file("${params.outdir}/no_intervals.bed"))
    } else if (!('annotate' in params.step.toLowerCase()) && !('controlfreec' in params.step.toLowerCase()))
        if (!params.intervals)
            result_intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS(result_fai))
        else
            result_intervals = CREATE_INTERVALS_BED(file(params.intervals))

    if (!params.no_intervals) {
        result_intervals = result_intervals.flatten()
            .map { intervalFile ->
                def duration = 0.0
                for (line in intervalFile.readLines()) {
                    final fields = line.split('\t')
                    if (fields.size() >= 5) duration += fields[4].toFloat()
                    else {
                        start = fields[1].toInteger()
                        end = fields[2].toInteger()
                        duration += (end - start) / params.nucleotides_per_second
                    }
                }
                [duration, intervalFile]
            }.toSortedList({ a, b -> b[0] <=> a[0] })
            .flatten().collate(2)
            .map{duration, intervalFile -> intervalFile}
    }

    emit:
        bwa                   = result_bwa
        bwa_version           = version_bwa
        dbsnp_tbi             = result_dbsnp_tbi
        dict                  = result_dict
        fai                   = result_fai
        germline_resource_tbi = result_germline_resource_tbi
        intervals             = result_intervals
        known_indels_tbi      = result_known_indels_tbi
        msisensorpro_scan     = result_msisensorpro_scan
        pon_tbi               = result_pon_tbi
        target_bed_gz_tbi     = result_target_bed
}