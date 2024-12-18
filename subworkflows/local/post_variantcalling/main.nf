//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_GERMLINE_VCFS } from '../vcf_concatenate_germline/main'
include { NORMALIZE_VCFS } from '../vcf_normalization/main'

workflow POST_VARIANTCALLING {

    take:
    germline_vcfs
    tumor_only_vcfs
    somatic_vcfs
    fasta
    concatenate_vcfs
    normalize_vcfs

    main:
    versions = Channel.empty()
    vcfs = Channel.empty()

    if (concatenate_vcfs){
        CONCATENATE_GERMLINE_VCFS(germline_vcfs)

        vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
        versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
    }

    if (normalize_vcfs){
        NORMALIZE_VCFS(germline_vcfs, tumor_only_vcfs, somatic_vcfs, fasta)

        vcfs = vcfs.mix(NORMALIZE_VCFS.out.vcfs)
        versions = versions.mix(NORMALIZE_VCFS.out.versions)
    }

    emit:
    vcfs // post processed vcfs

    versions // channel: [ versions.yml ]
}
