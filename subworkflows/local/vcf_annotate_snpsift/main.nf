//
// Run SnpSift annotateMem to annotate VCF files with multiple databases
//

include { SNPSIFT_ANNMEM } from '../../../modules/local/snpsift/annmem/main'

workflow VCF_ANNOTATE_SNPSIFT {
    take:
    ch_vcf          // channel: [val(meta), path(vcf)]
    ch_db_tuple     // channel: [[db_vcf], [db_vcf_tbi], [db_vardb], [db_fields], [db_prefixes]]

    main:
    // Add empty tbi placeholder to input VCF channel
    ch_vcf_input = ch_vcf.map { meta, vcf -> [meta, vcf, []] }

    // Annotate with all databases in one pass
    SNPSIFT_ANNMEM(
        ch_vcf_input,   // Sample VCF to annotate
        ch_db_tuple,    // Database(s) info
        false           // create = false for annotation mode
    )

    emit:
    vcf_tbi = SNPSIFT_ANNMEM.out.vcf  // channel: [val(meta), path(vcf), path(tbi)]
}
