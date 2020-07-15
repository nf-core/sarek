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
include { GATK_CREATE_SEQUENCE_DICTIONARY as GATK_CREATE_SEQUENCE_DICTIONARY } from '../local/gatk_dict.nf'
include { BUILD_INTERVALS } from '../local/build_intervals.nf'


workflow BUILD_INDICES{
    take:
        ch_dbsnp
        ch_fasta
        ch_germline_resource
        ch_known_indels
        ch_pon
        step

    main:

    if (!(params.bwa) && params.fasta && 'mapping' in step) 
        BWAMEM2_INDEX(ch_fasta)

    if (!(params.dict) && params.fasta && !('annotate' in step) && !('controlfreec' in step))
        GATK_CREATE_SEQUENCE_DICTIONARY(ch_fasta)

    if (!(params.fasta_fai) && params.fasta && !('annotate' in step))
        SAMTOOLS_FAIDX(ch_fasta)

    if (!(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || 'tnscope' in tools))
        HTSLIB_TABIX_DBSNP(ch_dbsnp) 
    
    if (!(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools)
        HTSLIB_TABIX_GERMLINE_RESOURCE(ch_germline_resource) 
    
    if (!(params.known_indels_index) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step))
        HTSLIB_TABIX_KNOWN_INDELS(ch_known_indels) 

    if (!(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools))
        HTSLIB_TABIX_PON(ch_pon)

    if (!(params.intervals) && !('annotate' in step) && !('controlfreec' in step)){
        ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : SAMTOOLS_FAIDX.out
        ch_fai.dump(tag: 'ch_fai')
        BUILD_INTERVALS(ch_fai)
    }

    emit:
        bwa_built = BWAMEM2_INDEX.out 
        dictBuilt = GATK_CREATE_SEQUENCE_DICTIONARY.out
        fai_built = SAMTOOLS_FAIDX.out
        dbsnp_tbi = HTSLIB_TABIX_DBSNP.out
        germline_resource_tbi = HTSLIB_TABIX_GERMLINE_RESOURCE.out
        known_indels_tbi = HTSLIB_TABIX_KNOWN_INDELS.out
        pon_tbi = HTSLIB_TABIX_PON.out
        intervalBuilt = BUILD_INTERVALS.out
}
