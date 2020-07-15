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
include { BWAMEM2_INDEX as BWAMEM2_INDEX } from '../nf-core/bwamem2_index.nf'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX } from '../nf-core/samtools_faidx.nf'
include { GATK_CREATE_SEQUENCE_DICTIONARY as GATK_CREATE_SEQUENCE_DICTIONARY } from './gatk_dict.nf'

workflow BUILD_INDICES{
    take:
    step
    ch_fasta
    ch_dbsnp
    ch_germline_resource
    ch_known_indels

    main:

    if(!(params.bwa) && params.fasta && 'mapping' in step) 
        BWAMEM2_INDEX(ch_fasta)

    if(!(params.dict) && params.fasta && !('annotate' in step) && !('controlfreec' in step))
        ATK_CREATE_SEQUENCE_DICTIONARY(ch_fasta)

    if(!(params.fasta_fai) && params.fasta && !('annotate' in step))
        SAMTOOLS_FAIDX(ch_fasta)

    if(!(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || 'tnscope' in tools))
        HTSLIB_TABIX_DBSNP(ch_dbsnp) 
    
    if(!(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools)
        HTSLIB_TABIX_GERMLINE_RESOURCE(ch_germline_resource) 
    
    if(!(params.known_indels_index) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step))
        HTSLIB_TABIX_KNOWN_INDELS(ch_known_indels) 

    emit:
    bwamem2_index = BWAMEM2_INDEX.out
    gatk_dict = GATK_CREATE_SEQUENCE_DICTIONARY.out
    samtools_faidx = SAMTOOLS_FAIDX.out
    tabix_dbsnp = HTSLIB_TABIX_DBSNP.out
    tabix_germline = HTSLIB_TABIX_GERMLINE_RESOURCE.out
    tabix_indels = HTSLIB_TABIX_KNOWN_INDELS.out

}


