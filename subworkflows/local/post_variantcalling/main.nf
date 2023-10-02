//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_GERMLINE_VCFS } from '../vcf_concatenate_germline/main'
include { VARLOCIRAPTOR_CALLS       } from '../vcf_varlociraptor_calls/main'

workflow POST_VARIANTCALLING {

    take:
    tools
    cram
    fasta
    fasta_fai
    vcfs
    concatenate_vcfs
    varlociraptor_scenario

    main:
    versions = Channel.empty()

    if(concatenate_vcfs){
        CONCATENATE_GERMLINE_VCFS(vcfs)

        vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
        versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
    }

    if(tools.split(',').contains('varlociraptor')) {
        VARLOCIRAPTOR_CALLS(vcfs, cram, fasta, fasta_fai, varlociraptor_scenario)
        vcfs = VARLOCIRAPTOR_CALLS.out.vcfs
        versions = versions.mix(VARLOCIRAPTOR_CALLS.out.versions)
    }
    emit:
    vcfs // post processed vcfs

    versions // channel: [ versions.yml ]
}
