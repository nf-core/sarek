//
// Run snpEff to annotate VCF files
//

params.bgziptabix_snpeff = [:]
params.snpeff_options    = [:]
params.snpeff_tag        = [:]
params.use_cache         = [:]

include { SNPEFF } from '../../modules/nf-core/modules/snpeff/main' addParams(
    options:    params.snpeff_options,
    snpeff_tag: params.snpeff_tag,
    use_cache:  params.use_cache
)

include { TABIX_BGZIPTABIX } from '../../modules/local/tabix/bgziptabix/main' addParams(options: params.bgziptabix_snpeff_options)

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
