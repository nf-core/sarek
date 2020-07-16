/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built

include { BUILD_INTERVALS }                 from '../local/build_intervals.nf'
include { BWAMEM2_INDEX }                   from '../nf-core/bwamem2_index.nf'
include { CREATE_INTERVALS_BED }            from '../local/create_intervals_bed.nf'
include { GATK_CREATE_SEQUENCE_DICTIONARY } from '../local/gatk_dict.nf'
include {
    HTSLIB_TABIX as HTSLIB_TABIX_DBSNP;
    HTSLIB_TABIX as HTSLIB_TABIX_GERMLINE_RESOURCE;
    HTSLIB_TABIX as HTSLIB_TABIX_KNOWN_INDELS;
    HTSLIB_TABIX as HTSLIB_TABIX_PON;
} from '../nf-core/htslib_tabix'
include { SAMTOOLS_FAIDX }                  from '../nf-core/samtools_faidx.nf'

workflow BUILD_INDICES{
    take:
        dbsnp
        fasta
        germline_resource
        known_indels
        pon
        step
        tools

    main:

    if (!(params.bwa) && params.fasta && 'mapping' in step)
        result_bwa = BWAMEM2_INDEX(fasta)
    else
        result_bwa = Channel.empty()

    if (!(params.dict) && params.fasta && !('annotate' in step) && !('controlfreec' in step))
        result_dict = GATK_CREATE_SEQUENCE_DICTIONARY(fasta)
    else
        result_dict = Channel.empty()

    if (!(params.fasta_fai) && params.fasta && !('annotate' in step))
        result_fai = SAMTOOLS_FAIDX(fasta)
    else
        result_fai = Channel.empty()

    if (!(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || 'tnscope' in tools))
        result_dbsnp_tbi = HTSLIB_TABIX_DBSNP(dbsnp)
    else
        result_dbsnp_tbi = Channel.empty()

    if (!(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools)
        result_germline_resource_tbi = HTSLIB_TABIX_GERMLINE_RESOURCE(germline_resource)
    else
        result_germline_resource_tbi = Channel.empty()

    if (!(params.known_indels_index) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step))
        result_known_indels_tbi = HTSLIB_TABIX_KNOWN_INDELS(known_indels)
    else
        result_known_indels_tbi = Channel.empty()

    if (!(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools))
        result_pon_tbi = HTSLIB_TABIX_PON(pon)
    else
        result_pon_tbi = Channel.empty()

    if (params.no_intervals) {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        result_intervals = Channel.from(file("${params.outdir}/no_intervals.bed"))
    } else if (!('annotate' in step) && !('controlfreec' in step))
        if (!params.intervals)
            intervals = CREATE_INTERVALS_BED(BUILD_INTERVALS(SAMTOOLS_FAIDX.out))
        else
            intervals = CREATE_INTERVALS_BED(params.intervals)

    if (!params.no_intervals) {
        intervals.flatten()
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
            .multiMap{
                all: it
                empty: ""
            }.set{bed}
        result_intervals = bed.all
    }

    emit:
        bwa                   = result_bwa
        dbsnp_tbi             = result_dbsnp_tbi
        dict                  = result_dict
        fai                   = result_fai
        germline_resource_tbi = result_germline_resource_tbi
        intervals_bed         = result_intervals
        known_indels_tbi      = result_known_indels_tbi
        pon_tbi               = result_pon_tbi
}
