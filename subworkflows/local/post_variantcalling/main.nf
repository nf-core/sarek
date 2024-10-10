//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_GERMLINE_VCFS } from '../vcf_concatenate_germline/main'
include { NORMALIZE_VCFS } from '../vcf_normalization/main'

workflow POST_VARIANTCALLING {

    take:
    vcfs
    fasta
    concatenate_vcfs
    normalized_vcfs
   
    main:
    versions = Channel.empty()

    if (concatenate_vcfs){
        CONCATENATE_GERMLINE_VCFS(vcfs, fasta)

        vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
        versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
    }

    if (normalized_vcfs){
        NORMALIZE_VCFS(vcfs, fasta)

        vcfs = vcfs.mix(NORMALIZE_VCFS.out.vcfs)
        versions = versions.mix(NORMALIZE_VCFS.out.versions)
    }

    emit:
    vcfs // post processed vcfs

    versions // channel: [ versions.yml ]
}
