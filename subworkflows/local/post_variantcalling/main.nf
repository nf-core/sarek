//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_GERMLINE_VCFS } from '../vcf_concatenate_germline'
include { NORMALIZE_VCFS            } from '../vcf_normalization'
include { VCF_VARLOCIRAPTOR         } from '../vcf_varlociraptor'

workflow POST_VARIANTCALLING {
    take:
    tools
    ch_cram
    germline_vcfs
    tumor_only_vcfs
    somatic_vcfs
    fasta
    fai
    concatenate_vcfs
    normalize_vcfs

    main:
    versions = Channel.empty()
    vcfs = Channel.empty()

    if (concatenate_vcfs) {
        CONCATENATE_GERMLINE_VCFS(germline_vcfs)

        vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
        versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
    }

    if (normalize_vcfs) {
        NORMALIZE_VCFS(germline_vcfs, tumor_only_vcfs, somatic_vcfs, fasta)

        vcfs = vcfs.mix(NORMALIZE_VCFS.out.vcfs)

        versions = versions.mix(NORMALIZE_VCFS.out.versions)
    }

    if(tools.split(',').contains('varlociraptor')) {
        VCF_VARLOCIRAPTOR(tools, ch_cram, fasta, fai, somatic_vcfs)
        vcfs = VCF_VARLOCIRAPTOR.out.vcf
        versions = versions.mix(VCF_VARLOCIRAPTOR.out.versions)
    }

    emit:
    vcfs     // post processed vcfs
    versions // channel: [ versions.yml ]
}
