//
// Run SNPEFF to annotate VCF files
//

include { SNPEFF_SNPEFF    } from '../../../modules/nf-core/snpeff/snpeff'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix'

workflow VCF_ANNOTATE_SNPEFF {
    take:
    ch_vcf // channel: [ val(meta), path(vcf) ]
    val_snpeff_db // string:  db version to use
    ch_snpeff_cache // channel: [ path(cache) ] (optional)

    main:
    SNPEFF_SNPEFF(ch_vcf, val_snpeff_db, ch_snpeff_cache)
    TABIX_BGZIPTABIX(SNPEFF_SNPEFF.out.vcf)


    emit:
    vcf_tbi   = TABIX_BGZIPTABIX.out.gz_index // channel: [ val(meta), path(vcf), path(tbi) ]
    reports   = SNPEFF_SNPEFF.out.report // channel: [ path(html) ]
    summary   = SNPEFF_SNPEFF.out.summary_html // channel: [ path(html) ]
    genes_txt = SNPEFF_SNPEFF.out.genes_txt // channel: [ path(genes.txt) ]
}
