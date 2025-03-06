//
// CONCATENATE Germline VCFs
//

// Concatenation of germline vcf-files
include { ADD_INFO_TO_VCF                                } from '../../../modules/local/add_info_to_vcf'
include { BCFTOOLS_CONCAT as GERMLINE_VCFS_CONCAT        } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT as GERMLINE_VCFS_CONCAT_SORT     } from '../../../modules/nf-core/bcftools/sort'
include { TABIX_BGZIPTABIX as TABIX_EXT_VCF              } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX as TABIX_GERMLINE_VCFS_CONCAT_SORT } from '../../../modules/nf-core/tabix/tabix'

workflow CONCATENATE_GERMLINE_VCFS {
    take:
    vcfs

    main:
    versions = Channel.empty()

    // Concatenate vcf-files
    ADD_INFO_TO_VCF(vcfs)
    TABIX_EXT_VCF(ADD_INFO_TO_VCF.out.vcf)

    // Gather vcfs and vcf-tbis for concatenating germline-vcfs
    germline_vcfs_with_tbis = TABIX_EXT_VCF.out.gz_tbi.map { meta, vcf, tbi -> [meta.subMap('id'), vcf, tbi] }.groupTuple()

    GERMLINE_VCFS_CONCAT(germline_vcfs_with_tbis)
    GERMLINE_VCFS_CONCAT_SORT(GERMLINE_VCFS_CONCAT.out.vcf)
    TABIX_GERMLINE_VCFS_CONCAT_SORT(GERMLINE_VCFS_CONCAT_SORT.out.vcf)

    // Gather versions of all tools used
    versions = versions.mix(ADD_INFO_TO_VCF.out.versions)
    versions = versions.mix(TABIX_EXT_VCF.out.versions)
    versions = versions.mix(GERMLINE_VCFS_CONCAT.out.versions)
    versions = versions.mix(GERMLINE_VCFS_CONCAT.out.versions)
    versions = versions.mix(TABIX_GERMLINE_VCFS_CONCAT_SORT.out.versions)

    emit:
    vcfs     = germline_vcfs_with_tbis // post processed vcfs
    versions // channel: [ versions.yml ]
}
