include { TABIX_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV } from '../../../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TIDDIT_SV                                 } from '../../../../modules/nf-core/modules/tiddit/sv/main'

workflow RUN_TIDDIT {
    take:
        cram_recalibrated
        fasta
        bwa

    main:

    ch_versions = Channel.empty()

    TIDDIT_SV(
        cram_recalibrated,
        fasta,
        bwa
    )

    TABIX_BGZIP_TIDDIT_SV(TIDDIT_SV.out.vcf)
    tiddit_ploidy = TIDDIT_SV.out.ploidy
    tiddit_vcf_gz = TABIX_BGZIP_TIDDIT_SV.out.gz_tbi.map{ meta, gz, tbi -> [meta, gz]}

    ch_versions = ch_versions.mix(TABIX_BGZIP_TIDDIT_SV.out.versions)
    ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)

    emit:
    versions = ch_versions

    tiddit_vcf = tiddit_vcf_gz
    tiddit_ploidy
}
