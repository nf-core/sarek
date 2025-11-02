//
// POST VARIANT CALLING: processes run on variantcalled but not annotated VCFs
//
include { BCFTOOLS_VIEW as FILTER_VCFS                             } from '../../../modules/nf-core/bcftools/view'
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
    germline_tbis
    cram_tumor_only
    tumor_only_vcfs
    tumor_only_tbis
    cram_somatic
    somatic_vcfs
    somatic_tbis
    fasta
    fai
    concatenate_vcfs
    filter_vcfs
    normalize_vcfs
    varlociraptor_chunk_size // integer: [mandatory] [default: 15] number of chunks to split BCF files when preprocessing and calling variants
    varlociraptor_scenario_germline
    varlociraptor_scenario_somatic
    varlociraptor_scenario_tumor_only

    main:
    versions = Channel.empty()
    vcfs = Channel.empty()
    tbis = Channel.empty()

    //
    // VARLOCIRAPTOR
    //
    if (tools && tools.split(',').contains('varlociraptor')) {
        // GERMLINE
        VCF_VARLOCIRAPTOR_GERMLINE(cram_germline, fasta, fai, varlociraptor_scenario_germline, germline_vcfs, varlociraptor_chunk_size, 'normal')

        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.vcf)
        tbis = tbis.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.tbi)
        versions = versions.mix(VCF_VARLOCIRAPTOR_GERMLINE.out.versions)

        // SOMATIC
        VCF_VARLOCIRAPTOR_SOMATIC(cram_somatic, fasta, fai, varlociraptor_scenario_somatic, somatic_vcfs, germline_vcfs, varlociraptor_chunk_size)

        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.vcf)
        tbis = tbis.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.tbi)
        versions = versions.mix(VCF_VARLOCIRAPTOR_SOMATIC.out.versions)

        // TUMOR ONLY
        VCF_VARLOCIRAPTOR_TUMOR_ONLY(cram_tumor_only, fasta, fai, varlociraptor_scenario_tumor_only, tumor_only_vcfs, varlociraptor_chunk_size, 'tumor')

        vcfs = vcfs.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.vcf)
        tbis = tbis.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.tbi)
        versions = versions.mix(VCF_VARLOCIRAPTOR_TUMOR_ONLY.out.versions)

    } else if (filter_vcfs || normalize_vcfs || concatenate_vcfs ) {

        all_vcfs = Channel.empty().mix(germline_vcfs, tumor_only_vcfs, somatic_vcfs)
        all_tbis = Channel.empty().mix(germline_tbis, tumor_only_tbis, somatic_tbis)

        // 1. Filter by PASS and custom fields
        // 2. Normalize
        // 3. Aggregate variants (Union, intersection, or n-1)
        if(filter_vcfs) {

            // Join VCFs with their corresponding TBIs before filtering
            FILTER_VCFS( all_vcfs.join(all_tbis, failOnDuplicate: true, failOnMismatch: true), [], [], [])

            all_vcfs = FILTER_VCFS.out.vcf
            all_tbis = FILTER_VCFS.out.tbi
            versions = versions.mix(FILTER_VCFS.out.versions)
        }

        if (normalize_vcfs) {

            NORMALIZE_VCFS(all_vcfs, fasta)

            all_vcfs = NORMALIZE_VCFS.out.vcfs // [meta, vcf]
            all_tbis = NORMALIZE_VCFS.out.tbis // [meta, tbi]
            versions = versions.mix(NORMALIZE_VCFS.out.versions)
        }

        if (concatenate_vcfs) {
            CONCATENATE_GERMLINE_VCFS(germline_vcfs)

            vcfs = vcfs.mix(CONCATENATE_GERMLINE_VCFS.out.vcfs)
            tbis = tbis.mix(CONCATENATE_GERMLINE_VCFS.out.tbis)
            versions = versions.mix(CONCATENATE_GERMLINE_VCFS.out.versions)
        }
    } else {
        // No post-processing requested, pass through original VCFs
        vcfs = vcfs.mix(germline_vcfs,tumor_only_vcfs, somatic_vcfs)
    }

    emit:
    vcfs     // post processed vcfs [meta, vcf]
    tbis     // post processed tbis [meta, tbi]
    versions // channel: [ versions.yml ]
}
