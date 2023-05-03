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
    vcf = TABIX_BGZIP_TIDDIT_SV.out.gz_tbi.map{ meta, gz, tbi -> [ meta + [ variantcaller: 'tiddit'], gz ] }

    versions = versions.mix(TABIX_BGZIP_TIDDIT_SV.out.versions)
    versions = versions.mix(TIDDIT_SV.out.versions)

    emit:
    ploidy
    vcf

    versions
}
