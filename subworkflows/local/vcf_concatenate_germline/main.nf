//
// CONCATENATE Germline VCFs
//

// Concatenation of germline vcf-files
include { ADD_INFO_TO_VCF                                     } from '../../../modules/local/add_info_to_vcf/main'
include { TABIX_BGZIPTABIX as TABIX_EXT_VCF                   } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_CONCAT  as GERMLINE_VCFS_CONCAT            } from '../../../modules/nf-core/bcftools/concat/main'
include { BCFTOOLS_SORT    as GERMLINE_VCFS_CONCAT_SORT       } from '../../../modules/nf-core/bcftools/sort/main'
include { TABIX_TABIX      as TABIX_GERMLINE_VCFS_CONCAT_SORT } from '../../../modules/nf-core/tabix/tabix/main'

workflow CONCATENATE_GERMLINE_VCFS {

    take:
    vcfs

    main:
    versions = Channel.empty()

    // Concatenate vcf-files
    ADD_INFO_TO_VCF(vcfs)
    TABIX_EXT_VCF(ADD_INFO_TO_VCF.out.vcf)

    // Gather vcfs and vcf-tbis for concatenating germline-vcfs
    germline_vcfs_with_tbis = TABIX_EXT_VCF.out.gz_tbi.groupTuple()

    GERMLINE_VCFS_CONCAT(germline_vcfs_with_tbis)
    GERMLINE_VCFS_CONCAT_SORT(GERMLINE_VCFS_CONCAT.out.vcf)
    TABIX_GERMLINE_VCFS_CONCAT_SORT(GERMLINE_VCFS_CONCAT_SORT.out.vcf)

    // Gather versions of all tools used
    versions = versions.mix(ADD_INFO_TO_VCF.out.versions)
    versions = versions.mix(TABIX_EXT_VCF.out.versions)
    versions = versions.mix(GERMLINE_VCFS_CONCAT.out.versions)
    versions = versions.mix(GERMLINE_VCFS_CONCAT_SORT.out.versions)
    versions = versions.mix(TABIX_GERMLINE_VCFS_CONCAT_SORT.out.versions)

    emit:
    vcfs = GERMLINE_VCFS_CONCAT_SORT.out.vcf // concatenated vcfs

    versions // channel: [ versions.yml ]
}

