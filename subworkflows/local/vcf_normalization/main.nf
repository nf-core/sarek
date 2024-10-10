// Normalize all unannotated VCFs

// Import modules
include { ADD_INFO_TO_VCF  } from '../../../modules/local/add_info_to_vcf/main'
include { TABIX_BGZIPTABIX as TABIX_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_NORM as VCFS_NORM } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SORT as VCFS_SORT } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX_VCFS_INDEX } from '../../../modules/nf-core/tabix/tabix/main'

// Workflow to normalize, compress, and index VCF files
workflow NORMALIZE_VCFS {

    take:
    vcfs
    fasta

    main:
    versions = Channel.empty()

    // Add additional information to VCF files
    ADD_INFO_TO_VCF(vcfs)
    
    // Normalize the VCF files with BCFTOOLS_NORM
    normalized_vcf = VCFS_NORM(vcf: ADD_INFO_TO_VCF.out.vcf)

    // Compress the normalized VCF files with bgzip
    compressed_vcf = TABIX_VCF(normalized_vcf)

    // Sort the compressed normalized VCF files
    sorted_vcf = VCFS_SORT(compressed_vcf)

    // Index the sorted VCF files
    sorted_indexed_vcf = TABIX_VCFS_INDEX(sorted_vcf)

    // Gather versions of all tools used
    versions = versions.mix(ADD_INFO_TO_VCF.out.versions)
    versions = versions.mix(VCFS_NORM.out.versions)
    versions = versions.mix(TABIX_VCF.out.versions)
    versions = versions.mix(VCFS_SORT.out.versions)
    versions = versions.mix(TABIX_VCFS_INDEX.out.versions)

    emit:
    normalized_vcfs = sorted_indexed_vcf // Post-processed sorted VCFs
    versions // Channel: [versions.yml]
}

