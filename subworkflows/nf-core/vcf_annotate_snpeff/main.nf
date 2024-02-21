//
// Run SNPEFF to annotate VCF files
//

include { SNPEFF_SNPEFF    } from '../../../modules/nf-core/snpeff/snpeff/main.nf'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main.nf'

workflow VCF_ANNOTATE_SNPEFF {
    take:
    ch_vcf          // channel: [ val(meta), path(vcf) ]
    val_snpeff_db   // string:  db version to use
    ch_snpeff_cache // channel: [ path(cache) ] (optional)

    main:
    ch_versions = Channel.empty()

    SNPEFF_SNPEFF(ch_vcf, val_snpeff_db, ch_snpeff_cache)
    TABIX_BGZIPTABIX(SNPEFF_SNPEFF.out.vcf)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SNPEFF_SNPEFF.out.versions)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions)

    emit:
    vcf_tbi  = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), path(vcf), path(tbi) ]
    reports  = SNPEFF_SNPEFF.out.report    // channel: [ path(html) ]
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}
