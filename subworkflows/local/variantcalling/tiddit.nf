include { TABIX_BGZIPTABIX as TABIX_BGZIP_TIDDIT_SV } from '../../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TIDDIT_SV                                 } from '../../../modules/nf-core/modules/tiddit/sv/main'

//TODO: UNDER CONSTRUCTIONS
workflow TIDDIT {
    take:


    main:

    ch_versions = Channel.empty()
    // if (tools.contains('tiddit')) {
    //     TODO: Update tiddit on bioconda, the current version does not support cram usage, needs newest version:
    //     https://github.com/SciLifeLab/TIDDIT/issues/82#issuecomment-1022103264
    //     Issue opened, either this week or end of february

    //     TIDDIT_SV(
    //         cram_recalibrated,
    //         fasta,
    //         fasta_fai
    //     )

    //     TABIX_BGZIP_TIDDIT_SV(TIDDIT_SV.out.vcf)
    //     tiddit_vcf_gz_tbi = TABIX_BGZIP_TIDDIT_SV.out.gz_tbi
    //     tiddit_ploidy     = TIDDIT_SV.out.ploidy
    //     tiddit_signals    = TIDDIT_SV.out.signals
    //     tiddit_wig        = TIDDIT_SV.out.wig
    //     tiddit_gc_wig     = TIDDIT_SV.out.gc_wig

    //     ch_versions = ch_versions.mix(TABIX_BGZIP_TIDDIT_SV.out.versions)
    //     ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)
    // }
    emit:
    versions = ch_versions
}
