/*
 * Run snpEff to annotate VCF files
 */

params.bgziptabix_snpeff = [:]
params.snpeff_options    = [:]
params.snpeff_tag        = [:]
params.use_cache         = [:]

include { SNPEFF } from '../../modules/nf-core/software/snpeff/main' addParams(
    options:    params.snpeff_options,
    snpeff_tag: params.snpeff_tag,
    use_cache:  params.use_cache
)

include { TABIX_BGZIPTABIX } from '../../modules/nf-core/software/tabix/bgziptabix/main' addParams(options: params.bgziptabix_snpeff_options)

workflow SNPEFF_ANNOTATE {
    take:
    vcf            // channel: [ val(meta), vcf, tbi ]
    snpeff_db      //   value: version of db to use
    snpeff_cache   //    path: path_to_snpeff_cache (optionnal)
    // skip   // boolean: true/false

    main:
    SNPEFF(vcf, snpeff_db, snpeff_cache)
    TABIX_BGZIPTABIX(SNPEFF.out.vcf)

    emit:
    vcf            = TABIX_BGZIPTABIX.out.tbi // channel: [ val(meta), vcf, tbi ]
    snpeff_reports = SNPEFF.out.reports       //    path: *.html
    snpeff_version = SNPEFF.out.version       //    path: *.version.txt
}
