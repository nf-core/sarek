// Normalize all unannotated VCFs

// Import modules
include { ADD_INFO_TO_VCF                   } from '../../../modules/local/add_info_to_vcf'
include { BCFTOOLS_NORM as VCFS_NORM        } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_SORT as VCFS_NORM_SORT   } from '../../../modules/nf-core/bcftools/sort'
include { TABIX_BGZIPTABIX as TABIX_EXT_VCF } from '../../../modules/nf-core/tabix/bgziptabix'

// Workflow to normalize, compress, and index VCF files
workflow NORMALIZE_VCFS {
    take:
    vcfs
    fasta

    main:
    versions = channel.empty()

    // Add additional information to VCF files
    ADD_INFO_TO_VCF(vcfs)

    // Compress the VCF files with bgzip
    TABIX_EXT_VCF(ADD_INFO_TO_VCF.out.vcf)

    // Normalize the VCF files with BCFTOOLS_NORM
    VCFS_NORM(TABIX_EXT_VCF.out.gz_index, fasta)

    // Sort the normalized VCF files
    VCFS_NORM_SORT(VCFS_NORM.out.vcf)

    // Gather versions of all tools used
    versions = versions.mix(ADD_INFO_TO_VCF.out.versions)
    versions = versions.mix(VCFS_NORM.out.versions)
    versions = versions.mix(VCFS_NORM_SORT.out.versions)

    emit:
    vcfs     = VCFS_NORM_SORT.out.vcf // normalized vcfs
    tbis     = VCFS_NORM_SORT.out.tbi // matching tbis
    versions // Channel: [versions.yml]
}
