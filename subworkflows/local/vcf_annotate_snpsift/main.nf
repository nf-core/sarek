//
// Run SnpSift annotateMem to annotate VCF files with multiple databases
//

include { SNPSIFT_ANNMEM } from '../../../modules/local/snpsift/annmem/main'

workflow VCF_ANNOTATE_SNPSIFT {
    take:
    ch_vcf          // channel: [val(meta), path(vcf)]
    ch_db_tuple     // channel: [[databases], [tbis], [vardbs], [fields], [prefixes]]

    main:
    // Add empty tbi and fields placeholder to input VCF channel
    // Input format: tuple(meta, vcf, tbi, fields)
    ch_vcf_input = ch_vcf.map { meta, vcf -> [meta, vcf, [], ''] }

    // Annotate with all databases in one pass
    SNPSIFT_ANNMEM(
        ch_vcf_input,
        ch_db_tuple
    )

    emit:
    vcf_tbi = SNPSIFT_ANNMEM.out.vcf  // channel: [val(meta), path(vcf), path(tbi)]
}
