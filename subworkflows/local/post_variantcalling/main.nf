//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_GERMLINE_VCFS } from '../vcf_concatenate_germline'
include { NORMALIZE_VCFS            } from '../vcf_normalization'
include { VCF_VARLOCIRAPTOR_SOMATIC } from '../vcf_varlociraptor_somatic'

workflow POST_VARIANTCALLING {
    take:
    tools
    cram_germline
    germline_vcfs
    cram_tumor_only
    tumor_only_vcfs
    cram_somatic
    somatic_vcfs
    fasta
    fai
    concatenate_vcfs
    normalize_vcfs
    varlociraptor_chunk_size      // integer: [mandatory] [default: 15] number of chunks to split BCF files when preprocessing and calling variants
    varlociraptor_scenario_file   // string: [mandatory] [default: "$projectDir/assets/varlociraptor_somatic_with_priors.yte.yaml"] path to the varlociraptor scenario file

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

    // implement for somatic case first
    if(tools.split(',').contains('varlociraptor')) {
        VCF_VARLOCIRAPTOR_SOMATIC(cram_somatic, fasta, fai, varlociraptor_scenario_file, somatic_vcfs, varlociraptor_chunk_size)
        vcfs = VCF_VARLOCIRAPTOR_SOMATIC.out.vcf
        versions = versions.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.versions)
    }

    // TODO: varlocirator for germline and tumor_only

    emit:
    vcfs     // post processed vcfs
    versions // channel: [ versions.yml ]
}
