/*
================================================================================
                                BUILDING INDICES
================================================================================
*/

params.build_intervals_options         = [:]
params.bwa_index_options               = [:]
params.bwamem2_index_options           = [:]
params.create_intervals_bed_options    = [:]
params.gatk_dict_options               = [:]
params.samtools_faidx_options          = [:]
params.tabix_dbsnp_options             = [:]
params.tabix_germline_resource_options = [:]
params.tabix_known_indels_options      = [:]
params.tabix_pon_options               = [:]

// Initialize channels based on params or indices that were just built

include { BUILD_INTERVALS }                            from '../process/build_intervals.nf'                           addParams(option: params.build_intervals_options)
include { BWA_INDEX }                                  from '../../nf-core/software/bwa/index/main.nf'                addParams(option: params.bwa_index_options)
include { BWAMEM2_INDEX }                              from '../../nf-core/software/bwamem2_index.nf'                 addParams(option: params.bwamem2_index_options)
include { CREATE_INTERVALS_BED }                       from '../process/create_intervals_bed.nf'                      addParams(option: params.create_intervals_bed_options)
include { GATK_CREATESEQUENCEDICTIONARY as GATK_DICT } from '../../nf-core/software/gatk/createsequencedictionary.nf' addParams(option: params.gatk_dict_options)
include { HTSLIB_TABIX as TABIX_DBSNP }                from '../../nf-core/software/htslib_tabix'                     addParams(option: params.tabix_dbsnp_options)
include { HTSLIB_TABIX as TABIX_GERMLINE_RESOURCE }    from '../../nf-core/software/htslib_tabix'                     addParams(option: params.tabix_germline_resource_options)
include { HTSLIB_TABIX as TABIX_KNOWN_INDELS }         from '../../nf-core/software/htslib_tabix'                     addParams(option: params.tabix_known_indels_options)
include { HTSLIB_TABIX as TABIX_PON }                  from '../../nf-core/software/htslib_tabix'                     addParams(option: params.tabix_pon_options)
include { SAMTOOLS_FAIDX }                             from '../../nf-core/software/samtools/faidx.nf'                addParams(option: params.samtools_faidx_options)

workflow BUILD_INDICES{
    take:
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon
        step              //   value: [mandatory] starting step
        tools             //    list: [optional]  tools to run

    main:

    result_bwa = Channel.empty()
    version_bwa = Channel.empty()
    if (!(params.bwa) && 'mapping' in step)
        if (params.aligner == "bwa-mem") (result_bwa, version_bwa) = BWA_INDEX(fasta)
        else                             result_bwa = BWAMEM2_INDEX(fasta)

    result_dict = Channel.empty()
    if (!(params.dict) && !('annotate' in step) && !('controlfreec' in step))
        result_dict = GATK_DICT(fasta)

    result_fai = Channel.empty()
    if (!(params.fasta_fai) && !('annotate' in step))
        result_fai = SAMTOOLS_FAIDX(fasta)

    result_dbsnp_tbi = Channel.empty()
    if (!(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || 'tnscope' in tools))
        result_dbsnp_tbi = TABIX_DBSNP(dbsnp)

    result_germline_resource_tbi = Channel.empty()
    if (!(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools)
        result_germline_resource_tbi = TABIX_GERMLINE_RESOURCE(germline_resource)

    result_known_indels_tbi = Channel.empty()
    if (!(params.known_indels_index) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step))
        result_known_indels_tbi = TABIX_KNOWN_INDELS(known_indels)

    result_pon_tbi = Channel.empty()
    if (!(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools))
        result_pon_tbi = TABIX_PON(pon)

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
        bwa                   = result_bwa
        bwa_version           = version_bwa
        dbsnp_tbi             = result_dbsnp_tbi
        dict                  = result_dict
        fai                   = result_fai
        germline_resource_tbi = result_germline_resource_tbi
        intervals             = result_intervals
        known_indels_tbi      = result_known_indels_tbi
        pon_tbi               = result_pon_tbi
}