//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//

include { CONCATENATE_GERMLINE_VCFS                                } from '../vcf_concatenate_germline'
include { NORMALIZE_VCFS                                           } from '../vcf_normalization'
include { VCF_VARLOCIRAPTOR_SINGLE as VCF_VARLOCIRAPTOR_GERMLINE   } from '../vcf_varlociraptor_single'
include { VCF_VARLOCIRAPTOR_SOMATIC                                } from '../vcf_varlociraptor_somatic'
include { VCF_VARLOCIRAPTOR_SINGLE as VCF_VARLOCIRAPTOR_TUMOR_ONLY } from '../vcf_varlociraptor_single'

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

    //
    // VARLOCIRAPTOR
    //
    if(tools.split(',').contains('varlociraptor') && cram_germline) {
        VCF_VARLOCIRAPTOR_GERMLINE(cram_germline, fasta, fai, varlociraptor_scenario_file, germline_vcfs, varlociraptor_chunk_size, 'normal')
        vcfs = VCF_VARLOCIRAPTOR_GERMLINE.out.vcf
        versions = versions.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.versions)
    }
    if(tools.split(',').contains('varlociraptor') && cram_somatic) {
        VCF_VARLOCIRAPTOR_SOMATIC(cram_somatic, fasta, fai, varlociraptor_scenario_file, somatic_vcfs, varlociraptor_chunk_size)
        vcfs = VCF_VARLOCIRAPTOR_SOMATIC.out.vcf
        versions = versions.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.versions)
    }
    if(tools.split(',').contains('varlociraptor') && cram_tumor_only) {
        VCF_VARLOCIRAPTOR_TUMOR_ONLY(cram_tumor_only, fasta, fai, varlociraptor_scenario_file, tumor_only_vcfs, varlociraptor_chunk_size, 'tumor')
        vcfs = VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.vcf
        versions = versions.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.versions)
    }

    emit:
    vcfs     // post processed vcfs
    versions // channel: [ versions.yml ]
}
