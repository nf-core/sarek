//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_GERMLINE_VCFS } from '../vcf_concatenate_germline/main'

include { checkInParam } from "${projectDir}/checkInParam"

workflow POST_VARIANTCALLING {

    take:
    vcfs
    concatenate_vcfs

    main:
    versions = Channel.empty()

    if(concatenate_vcfs){
        CONCATENATE_GERMLINE_VCFS(vcfs)

        vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
        versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
    }

    emit:
    vcfs // post processed vcfs

    versions // channel: [ versions.yml ]
}
