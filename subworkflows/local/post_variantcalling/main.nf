//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_VCFS                                     } from '../../../subworkflows/vcf_concatenate_germline/main'

workflow CONCATENATE_VCFS {
    take:
    vcfs
    concatenate_vcfs

    main:
    versions = Channel.empty()

    if(concatenate_vcfs){
        CONCATENATE_VCFS(vcfs)
        versions = versions.mix(CONCATENATE_VCFS.out.versions)
    }


    emit:
    vcfs = CONCATENATE_VCFS.out.vcfs // post processed vcfs

    versions // channel: [ versions.yml ]
}
