include { BGZIP as BGZIP_VC_DEEPVARIANT_GVCF        } from '../../../modules/local/bgzip'
include { BGZIP as BGZIP_VC_DEEPVARIANT_VCF         } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_DEEPVARIANT_GVCF     } from '../../../modules/local/concat_vcf/main'
include { CONCAT_VCF as CONCAT_DEEPVARIANT_VCF      } from '../../../modules/local/concat_vcf/main'
include { DEEPVARIANT                               } from '../../../modules/nf-core/modules/deepvariant/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_GVCF  } from '../../../modules/nf-core/modules/tabix/tabix/main'
include { TABIX_TABIX as TABIX_VC_DEEPVARIANT_VCF   } from '../../../modules/nf-core/modules/tabix/tabix/main'

//TODO: benchmark if it is better to provide multiple bed files & run on multiple machines + mergeing afterwards || one containing all intervals and run on one larger machine
// Deepvariant: https://github.com/google/deepvariant/issues/510
workflow RUN_DEEPVARIANT {
    take:
    cram_recalibrated_intervals            // channel: [mandatory]
    fasta
    fasta_fai
    intervals_bed_combine_gz
    num_intervals

    main:

    ch_versions = Channel.empty()

    DEEPVARIANT(
        cram_recalibrated_intervals,
        fasta,
        fasta_fai)

    // Only when no intervals
    TABIX_VC_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
    TABIX_VC_DEEPVARIANT_GVCF(DEEPVARIANT.out.gvcf)

    // Only when using intervals
    BGZIP_VC_DEEPVARIANT_VCF(DEEPVARIANT.out.vcf)
    BGZIP_VC_DEEPVARIANT_GVCF(DEEPVARIANT.out.gvcf)

    CONCAT_DEEPVARIANT_VCF(
        BGZIP_VC_DEEPVARIANT_VCF.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    CONCAT_DEEPVARIANT_GVCF(
        BGZIP_VC_DEEPVARIANT_GVCF.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    deepvariant_vcf = channel.empty().mix(
        CONCAT_DEEPVARIANT_GVCF.out.vcf,
        CONCAT_DEEPVARIANT_VCF.out.vcf,
        DEEPVARIANT.out.gvcf.join(TABIX_VC_DEEPVARIANT_GVCF.out.tbi),
        DEEPVARIANT.out.vcf.join(TABIX_VC_DEEPVARIANT_VCF.out.tbi))

    ch_versions = ch_versions.mix(BGZIP_VC_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(BGZIP_VC_DEEPVARIANT_VCF.out.versions)
    ch_versions = ch_versions.mix(CONCAT_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(CONCAT_DEEPVARIANT_VCF.out.versions)
    ch_versions = ch_versions.mix(DEEPVARIANT.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_GVCF.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_DEEPVARIANT_VCF.out.versions)

    emit:
    deepvariant_vcf
    versions = ch_versions
}
