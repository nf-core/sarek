/*
================================================================================
                                BUILDING INDICES
================================================================================
*/

// And then initialize channels based on params or indices that were just built

include { BUILD_INTERVALS }                            from '../process/build_intervals.nf'
include { BWA_INDEX }                                  from '../../nf-core/software/bwa/index/main.nf'
include { BWAMEM2_INDEX }                              from '../../nf-core/software/bwamem2_index.nf'
include { CREATE_INTERVALS_BED }                       from '../process/create_intervals_bed.nf'
include { GATK_CREATESEQUENCEDICTIONARY as GATK_DICT } from '../../nf-core/software/gatk/createsequencedictionary.nf'
include { HTSLIB_TABIX as TABIX_DBSNP;
          HTSLIB_TABIX as TABIX_GERMLINE_RESOURCE;
          HTSLIB_TABIX as TABIX_KNOWN_INDELS;
          HTSLIB_TABIX as TABIX_PON;}                  from '../../nf-core/software/htslib_tabix'
include { SAMTOOLS_FAIDX }                             from '../../nf-core/software/samtools/faidx.nf'

workflow BUILD_INDICES{
    take:
        dbsnp             // channel: [optional]  dbsnp
        fasta             // channel: [mandatory] fasta
        germline_resource // channel: [optional]  germline_resource
        known_indels      // channel: [optional]  known_indels
        pon               // channel: [optional]  pon
        step              //   value: [mandatory] starting step
        tools             //    list: [optional]  tools to run
        bwa_index_opts    //     map: options for BWA_INDEX module

    main:

    result_bwa = Channel.empty()
    version_bwa = Channel.empty()
    if (!(params.bwa) && 'mapping' in step)
        if (params.aligner == "bwa-mem") (result_bwa, version_bwa) = BWA_INDEX(fasta, bwa_index_opts)
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
            result_intervals = CREATE_INTERVALS_BED(result_fai)
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
