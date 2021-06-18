/*
 * Run snpEff to annotate VCF files
 */

params.snpeff_options     = [:]
params.bgziptabix_snpeff  = [:]

include { SNPEFF }           from '../../modules/nf-core/software/snpeff/main'           addParams(options: params.snpeff_options)
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/software/tabix/bgziptabix/main' addParams(options: params.bgziptabix_snpeff)

workflow SNPEFF_ANNOTATE {
    take:
    vcf            // channel: [ val(meta), vcf, tbi ]
    snpeff_db      //   value: version of db to use
    use_cache      //    bool: use cache (default: false)
    snpeff_cache   //    path: path_to_snpeff_cache (optionnal)
    snpeff_tag     //   value: tag of docker image (optionnal)
    // skip   // boolean: true/false

    main:
    params.use_cache = false
    use_cache = params.use_cache
    params.snpeff_tag = false
    snpeff_tag = params.snpeff_tag

    SNPEFF(vcf, snpeff_db, use_cache, snpeff_cache, snpeff_tag)
    TABIX_BGZIPTABIX(SNPEFF.out.vcf)

    emit:
    vcf            = TABIX_BGZIPTABIX.out.tbi // channel: [ val(meta), vcf, tbi ]
    snpeff_reports = SNPEFF.out.reports       //    path: *.html
    snpeff_version = SNPEFF.out.version       //    path: *.version.txt
}
