include { BGZIP as BGZIP_VC_FREEBAYES               } from '../../../modules/local/bgzip'
include { CONCAT_VCF as CONCAT_FREEBAYES            } from '../../../modules/local/concat_vcf/main'
include { FREEBAYES                                 } from '../../../modules/nf-core/modules/freebayes/main'
include { TABIX_TABIX as TABIX_VC_FREEBAYES         } from '../../../modules/nf-core/modules/tabix/tabix/main'


workflow RUN_FREEBAYES {
    take:
    cram_recalibrated_intervals_freebayes
    fasta
    fasta_fai

    main:

    ch_versions = Channel.empty()

    FREEBAYES(
        cram_recalibrated_intervals_freebayes,
        fasta,
        fasta_fai,
        [], [], [])

    // Only when no intervals
    TABIX_VC_FREEBAYES(FREEBAYES.out.vcf)

    // Only when using intervals
    BGZIP_VC_FREEBAYES(FREEBAYES.out.vcf)

    CONCAT_FREEBAYES(
        BGZIP_VC_FREEBAYES.out.vcf
            .map{ meta, vcf ->
                new_meta = meta.clone()
                new_meta.id = new_meta.sample
                [new_meta, vcf]
            }.groupTuple(size: num_intervals),
        fasta_fai,
        intervals_bed_combine_gz)

    freebayes_vcf = Channel.empty().mix(
        CONCAT_FREEBAYES.out.vcf,
        FREEBAYES.out.vcf.join(TABIX_VC_FREEBAYES.out.tbi))

    ch_versions = ch_versions.mix(BGZIP_VC_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(CONCAT_FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    ch_versions = ch_versions.mix(TABIX_VC_FREEBAYES.out.versions)

    emit:
    versions = ch_versions
}
