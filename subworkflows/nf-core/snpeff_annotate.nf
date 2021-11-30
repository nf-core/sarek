//
// Run snpEff to annotate VCF files
//

include { SNPEFF           } from '../../modules/nf-core/modules/snpeff/main'
include { TABIX_BGZIPTABIX } from '../../modules/local/tabix/bgziptabix/main'

workflow SNPEFF_ANNOTATE {
    take:
    vcf            // channel: [ val(meta), vcf ]
    snpeff_db      //   value: version of db to use
    snpeff_cache   //    path: path_to_snpeff_cache (optionnal)

    main:
    SNPEFF(vcf, snpeff_db, snpeff_cache)
    TABIX_BGZIPTABIX(SNPEFF.out.vcf)

    emit:
    vcf_tbi  = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = SNPEFF.out.report           //    path: *.html
    versions = SNPEFF.out.versions         //    path: versions.yml
}
