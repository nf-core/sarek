// CONCATENATE Germline VCFs

// Concatenation of germline vcf-files
include { ADD_INFO_TO_VCF  } from '../../../modules/local/add_info_to_vcf/main'
include { TABIX_BGZIPTABIX as TABIX_EXT_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_NORM as GERMLINE_VCFS_NORM } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_CONCAT as GERMLINE_VCFS_CONCAT } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT as GERMLINE_VCFS_CONCAT_SORT } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX_GERMLINE_VCFS_CONCAT_SORT } from '../../../modules/nf-core/tabix/tabix/main'

workflow CONCATENATE_GERMLINE_VCFS {

    take:
    vcfs
    fasta

    main:
    versions = Channel.empty()

    // Add additional information to VCF files
    ADD_INFO_TO_VCF(vcfs)
    
    // Compress the VCF files with bgzip
    TABIX_EXT_VCF(ADD_INFO_TO_VCF.out.vcf)

    // Normalize the VCF files with BCFTOOLS_NORM
    GERMLINE_VCFS_NORM(vcf: ADD_INFO_TO_VCF.out.vcf, fasta: fasta)

    // Compress the normalized VCF files with bgzip
    TABIX_EXT_VCF(GERMLINE_VCFS_NORM.out.vcf)

    // Index the compressed normalized VCF files
    TABIX_GERMLINE_VCFS_CONCAT_SORT(TABIX_EXT_VCF.out.gz)

    // Gather vcfs and vcf-tbis for concatenating germline-vcfs
    germline_vcfs_with_tbis = TABIX_GERMLINE_VCFS_CONCAT_SORT.out.map { meta, vcf, tbi -> [meta.subMap('id'), vcf, tbi] }.groupTuple()

    // Concatenate the VCF files
    GERMLINE_VCFS_CONCAT(germline_vcfs_with_tbis)
    
    // Sort the concatenated VCF files
    GERMLINE_VCFS_CONCAT_SORT(GERMLINE_VCFS_CONCAT.out.vcf)
    
    // Index the sorted concatenated VCF files
    TABIX_GERMLINE_VCFS_CONCAT_SORT(GERMLINE_VCFS_CONCAT_SORT.out.vcf)

    // Gather versions of all tools used
    versions = versions.mix(ADD_INFO_TO_VCF.out.versions)
    versions = versions.mix(TABIX_EXT_VCF.out.versions)
    versions = versions.mix(GERMLINE_VCFS_NORM.out.versions)
    versions = versions.mix(GERMLINE_VCFS_CONCAT.out.versions)
    versions = versions.mix(GERMLINE_VCFS_CONCAT_SORT.out.versions)
    versions = versions.mix(TABIX_GERMLINE_VCFS_CONCAT_SORT.out.versions)

    emit:
    vcfs = TABIX_GERMLINE_VCFS_CONCAT_SORT.out.gz_tbi // post-processed VCFs
    versions // channel: [ versions.yml ]
}