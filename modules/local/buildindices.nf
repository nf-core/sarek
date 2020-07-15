/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built


include { HTSLIB_TABIX as HTSLIB_TABIX_DBSNP } from '../nf-core/htslib_tabix'
include { HTSLIB_TABIX as HTSLIB_TABIX_GERMLINE_RESOURCE } from '../nf-core/htslib_tabix'
include { HTSLIB_TABIX as HTSLIB_TABIX_KNOWN_INDELS } from '../nf-core/htslib_tabix'
include { HTSLIB_TABIX as HTSLIB_TABIX_PON } from '../nf-core/htslib_tabix'

workflow BUILD_INDICES{
    take:
    ch_fasta
    ch_dbsnp
    ch_germline_resource
    ch_known_indels

    main:
    BWAMEM2_INDEX(ch_fasta)
    GATK_CREATE_SEQUENCE_DICTIONARY(ch_fasta)
    SAMTOOLS_FAIDX(ch_fasta)
    HTSLIB_TABIX_DBSNP(ch_dbsnp) //ch_dbsnp
    HTSLIB_TABIX_GERMLINE_RESOURCE(ch_germline_resource) //ch_germline_resource
    HTSLIB_TABIX_KNOWN_INDELS(ch_known_indels) //ch_knwon_indels

    emit:
    bwamem2_index = BWAMEM2_INDEX.out
    gatk_dict = GATK_CREATE_SEQUENCE_DICTIONARY.out
    samtools_faidx = SAMTOOLS_FAIDX.out
    tabix_dbsnp = HTSLIB_TABIX_DBSNP.out
    tabix_germline = HTSLIB_TABIX_GERMLINE_RESOURCE.out
    tabix_indels = HTSLIB_TABIX_KNOWN_INDELS.out

}


