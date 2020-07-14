/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built


include HTSLIB_TABIX as HTSLIB_TABIX_DBSNP from '../../nf-core/htslib_tabix'
include HTSLIB_TABIX as HTSLIB_TABIX_GERMLINE_RESOURCE from '../../nf-core/htslib_tabix'
include HTSLIB_TABIX as HTSLIB_TABIX_KNOWN_INDELS from '../../nf-core/htslib_tabix'
include HTSLIB_TABIX as HTSLIB_TABIX_PON from '../../nf-core/htslib_tabix'

ch_known_indels_tbi = params.known_indels ? params.known_indels_index ? Channel.value(file(params.known_indels_index)) : known_indels_tbi.collect() : "null"


ch_pon_tbi = params.pon ? params.pon_index ? Channel.value(file(params.pon_index)) : pon_tbi : "null"



ch_intervals = params.no_intervals ? "null" : params.intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : intervalBuilt

workflow build_indices{

    BWAMEM2_INDEX(ch_fasta)

    GATK_CREATE_SEQUENCE_DICTIONARY(ch_fasta)

    SAMTOOLS_FAIDX(ch_fasta)

    HTSLIB_TABIX_DBSNP(ch_dbsnp) //ch_dbsnp

    HTSLIB_TABIX_GERMLINE_RESOURCE(ch_germline_resource) //ch_germline_resource

    HTSLIB_TABIX_KNOWN_INDELS(ch_known_indels) //ch_knwon_indels

}

ch_bwa = params.bwa ? Channel.value(file(params.bwa)) : BWAMEM2_INDEX.out

ch_dict = params.dict ? Channel.value(file(params.dict)) : GATK_CREATE_SEQUENCE_DICTIONARY.out

ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : SAMTOOLS_FAIDX.out

ch_dbsnp_tbi = params.dbsnp ? params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : dbsnp_tbi : "null"

ch_germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) : germline_resource_tbi : "null"

