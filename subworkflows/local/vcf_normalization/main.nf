// Normalize all unannotated VCFs

// Import modules
include { ADD_INFO_TO_VCF  } from '../../../modules/local/add_info_to_vcf/main'
include { TABIX_BGZIPTABIX as TABIX_EXT_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_NORM as VCFS_NORM } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SORT as VCFS_NORM_SORT } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX as TABIX_VCFS_NORM_SORT } from '../../../modules/nf-core/tabix/tabix/main'

// Workflow to normalize, compress, and index VCF files
workflow NORMALIZE_VCFS {

    take:
    germline_vcfs
    tumor_only_vcfs
    somatic_vcfs
    fasta

    main:
    versions = Channel.empty()

    vcfs = germline_vcfs.mix(tumor_only_vcfs, somatic_vcfs)

    // Add additional information to VCF files
    ADD_INFO_TO_VCF(vcfs)

    // Compress the VCF files with bgzip
    TABIX_EXT_VCF(ADD_INFO_TO_VCF.out.vcf)

    // Normalize the VCF files with BCFTOOLS_NORM
    VCFS_NORM(TABIX_EXT_VCF.out.gz_tbi, fasta)

    // Sort the normalized VCF files
    VCFS_NORM_SORT(VCFS_NORM.out.vcf)

    // Index the sorted normalized VCF files
    TABIX_VCFS_NORM_SORT(VCFS_NORM_SORT.out.vcf)

    // Gather versions of all tools used
    versions = versions.mix(ADD_INFO_TO_VCF.out.versions)
    versions = versions.mix(VCFS_NORM.out.versions)
    versions = versions.mix(TABIX_EXT_VCF.out.versions)
    versions = versions.mix(VCFS_NORM_SORT.out.versions)
    versions = versions.mix(TABIX_VCFS_NORM_SORT.out.versions)

    emit:
    vcfs = VCFS_NORM_SORT.out.vcf // normalized vcfs
    versions // Channel: [versions.yml]
}

