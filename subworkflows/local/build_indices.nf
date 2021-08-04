/*
========================================================================================
    BUILDING INDICES
========================================================================================
*/

params.bgziptabix_options              = [:]
params.build_intervals_options         = [:]
params.bwa_index_options               = [:]
params.bwamem2_index_options           = [:]
params.create_intervals_bed_options    = [:]
params.gatk4_dict_options              = [:]
params.msisensorpro_scan_options       = [:]
params.samtools_faidx_options          = [:]
params.tabix_dbsnp_options             = [:]
params.tabix_germline_resource_options = [:]
params.tabix_known_indels_options      = [:]
params.tabix_pon_options               = [:]

// Initialize channels based on params or indices that were just built

include { BUILD_INTERVALS }                              from '../../modules/local/build_intervals/main'                           addParams(options: params.build_intervals_options)
include { BWA_INDEX as BWAMEM1_INDEX }                   from '../../modules/nf-core/modules/bwa/index/main'                       addParams(options: params.bwa_index_options)
include { BWAMEM2_INDEX }                                from '../../modules/nf-core/modules/bwamem2/index/main'                   addParams(options: params.bwamem2_index_options)
include { CREATE_INTERVALS_BED }                         from '../../modules/local/create_intervals_bed/main'                      addParams(options: params.create_intervals_bed_options)
include { GATK4_CREATESEQUENCEDICTIONARY as GATK4_DICT } from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'  addParams(options: params.gatk4_dict_options)
include { MSISENSORPRO_SCAN }                            from '../../modules/local/msisensorpro/scan/main'                         addParams(options: params.msisensorpro_scan_options)
include { SAMTOOLS_FAIDX }                               from '../../modules/nf-core/modules/samtools/faidx/main'                  addParams(options: params.samtools_faidx_options)
include { TABIX_BGZIPTABIX }                             from '../../modules/nf-core/modules/tabix/bgziptabix/main'                addParams(options: params.bgziptabix_options)
include { TABIX_TABIX as TABIX_DBSNP }                   from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_dbsnp_options)
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE }       from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_germline_resource_options)
include { TABIX_TABIX as TABIX_KNOWN_INDELS }            from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_known_indels_options)
include { TABIX_TABIX as TABIX_PON }                     from '../../modules/nf-core/modules/tabix/tabix/main'                     addParams(options: params.tabix_pon_options)

workflow BUILD_INDICES {
    take:
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        fasta_fai         // channel: [optional]  fasta_fai
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon
        target_bed        // channel: [optionnal] target_bed
        tools
        step

    main:

    result_bwa  = Channel.empty()
    version_bwa = Channel.empty()
    if (!(params.bwa) && 'mapping' in step)
        if (params.aligner == "bwa-mem") (result_bwa, version_bwa) = BWAMEM1_INDEX(fasta)
        else                             (result_bwa, version_bwa) = BWAMEM2_INDEX(fasta)

    result_dict = Channel.empty()
    version_dict = Channel.empty()
    if (!(params.dict) && !('annotate' in step) && !('controlfreec' in step))
        (result_dict, version_dict) = GATK4_DICT(fasta)

    result_fai = Channel.empty()
    if (fasta_fai) result_fai = fasta_fai
    version_fai = Channel.empty()
    if (!(params.fasta_fai) && !('annotate' in step))
        (result_fai, version_fai) = SAMTOOLS_FAIDX(fasta)

    result_dbsnp_tbi = Channel.empty()
    version_dbsnp_tbi = Channel.empty()
    if (!(params.dbsnp_tbi) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools|| 'mutect2' in tools || 'tnscope' in tools)) {
        dbsnp_id = dbsnp.map {it -> [[id:"$it.baseName"], it]}
        (result_dbsnp_tbi, version_dbsnp_tbi) = TABIX_DBSNP(dbsnp_id)
        result_dbsnp_tbi = result_dbsnp_tbi.map {meta, tbi -> [tbi]}
    }

    result_target_bed = Channel.empty()
    version_target_bed = Channel.empty()
    if ((params.target_bed) && ('manta' in tools || 'strelka' in tools)) {
        target_bed_id = target_bed.map {it -> [[id:"$it.baseName"], it]}
        (result_target_bed, version_target_bed) = TABIX_BGZIPTABIX(target_bed_id)
        result_target_bed = result_target_bed.map {meta, bed, tbi -> [bed, tbi]}
    }

    result_germline_resource_tbi = Channel.empty()
    version_germline_resource_tbi = Channel.empty()
    if (!(params.germline_resource_tbi) && params.germline_resource && 'mutect2' in tools) {
        germline_resource_id = germline_resource.map {it -> [[id:"$it.baseName"], it]}
        (result_germline_resource_tbi, version_germline_resource_tbi) = TABIX_GERMLINE_RESOURCE(germline_resource_id)
    }

    result_known_indels_tbi = Channel.empty()
    version_known_indels_tbi = Channel.empty()
    if (!(params.known_indels_tbi) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step)) {
        known_indels_id = known_indels.map {it -> [[id:"$it.baseName"], it]}
        (result_known_indels_tbi, version_known_indels_tbi) = TABIX_KNOWN_INDELS(known_indels_id)
        result_known_indels_tbi = result_known_indels_tbi.map {meta, tbi -> [tbi]}
    }

    result_msisensorpro_scan = Channel.empty()
    version_msisensorpro_scan = Channel.empty()
    if ('msisensorpro' in tools)
        (result_msisensorpro_scan, version_msisensorpro_scan) = MSISENSORPRO_SCAN(fasta)

    result_pon_tbi = Channel.empty()
    version_pon_tbi = Channel.empty()
    if (!(params.pon_tbi) && params.pon && ('tnscope' in tools || 'mutect2' in tools)) {
        pon_id = pon.map {it -> [[id:"$it.baseName"], it]}
        (result_pon_tbi, version_pon_tbi) = TABIX_PON(pon_id)
    }

    result_intervals = Channel.empty()
    if (params.no_intervals) {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        result_intervals = Channel.from(file("${params.outdir}/no_intervals.bed"))
    } else if (!('annotate' in step) && !('controlfreec' in step))
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
        bwa                           = result_bwa
        bwa_version                   = version_bwa
        dbsnp_tbi                     = result_dbsnp_tbi
        dbsnp_tbi_version             = version_dbsnp_tbi
        dict                          = result_dict
        dict_version                  = version_dict
        fai                           = result_fai
        fai_version                   = version_fai
        germline_resource_tbi         = result_germline_resource_tbi
        germline_resource_tbi_version = version_germline_resource_tbi
        intervals                     = result_intervals
        known_indels_tbi              = result_known_indels_tbi.collect()
        known_indels_tbi_version      = version_known_indels_tbi
        msisensorpro_scan             = result_msisensorpro_scan
        msisensorpro_scan_version     = version_msisensorpro_scan
        pon_tbi                       = result_pon_tbi
        pon_tbi_version               = version_pon_tbi
        target_bed_gz_tbi             = result_target_bed
        target_bed_version            = version_target_bed
}
