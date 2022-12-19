include { TABIX_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TIDDIT_SV                                 } from '../../../modules/nf-core/tiddit/sv/main'

workflow BAM_VARIANT_CALLING_SINGLE_TIDDIT {
    take:
        cram
        fasta
        bwa

    main:
    versions = Channel.empty()

    TIDDIT_SV(cram, fasta, bwa)

    TABIX_BGZIP_TIDDIT_SV(TIDDIT_SV.out.vcf)

    ploidy = TIDDIT_SV.out.ploidy
    vcf = TABIX_BGZIP_TIDDIT_SV.out.gz_tbi.map{ meta, gz, tbi ->

        new_meta = meta.tumor_id ? [
                                        id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                        normal_id:      meta.normal_id,
                                        num_intervals:  meta.num_intervals,
                                        patient:        meta.patient,
                                        sex:            meta.sex,
                                        tumor_id:       meta.tumor_id,
                                        variantcaller:  'tiddit'
                                    ]
                                    : [
                                        id:             meta.sample,
                                        num_intervals:  meta.num_intervals,
                                        patient:        meta.patient,
                                        sample:         meta.sample,
                                        sex:            meta.sex,
                                        status:         meta.status,
                                        variantcaller:  'tiddit'
                                    ]
        [ new_meta, gz ] }

    versions = versions.mix(TABIX_BGZIP_TIDDIT_SV.out.versions)
    versions = versions.mix(TIDDIT_SV.out.versions)

    emit:
    ploidy
    vcf

    versions
}
