//
// CONCATENATE Germline VCFs
//

// Concatenation of germline vcf-files
include { ADD_INFO_TO_VCF                                } from '../../../modules/local/add_info_to_vcf'
include { BCFTOOLS_CONCAT as GERMLINE_VCFS_CONCAT        } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT as GERMLINE_VCFS_CONCAT_SORT     } from '../../../modules/nf-core/bcftools/sort'
include { HTSLIB_BGZIPTABIX as TABIX_EXT_VCF             } from '../../../modules/nf-core/htslib/bgziptabix'

workflow CONCATENATE_GERMLINE_VCFS {
    take:
    vcfs

    main:
    // Concatenate vcf-files
    ADD_INFO_TO_VCF(vcfs)
    TABIX_EXT_VCF(ADD_INFO_TO_VCF.out.vcf.map{ meta, vcf -> [ meta, vcf, [], [] ] }, 'compress', true, 'vcf')

    // Gather vcfs and vcf-tbis for concatenating germline-vcfs
    germline_vcfs_with_tbis = TABIX_EXT_VCF.out.output.join(TABIX_EXT_VCF.out.index).groupTuple()

    GERMLINE_VCFS_CONCAT(germline_vcfs_with_tbis)
    GERMLINE_VCFS_CONCAT_SORT(GERMLINE_VCFS_CONCAT.out.vcf)

    emit:
    vcfs     = GERMLINE_VCFS_CONCAT_SORT.out.vcf // concatenated vcfs
    tbis     = GERMLINE_VCFS_CONCAT_SORT.out.index // matching tbis
}
