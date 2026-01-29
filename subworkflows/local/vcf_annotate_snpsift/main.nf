//
// Run SnpSift annotateMem to annotate VCF files with multiple databases
//

include { SNPSIFT_ANNMEM } from '../../../modules/local/snpsift/annmem/main'

workflow VCF_ANNOTATE_SNPSIFT {
    take:
    ch_vcf          // channel: [val(meta), path(vcf)]
    ch_db_configs   // channel: [[vcf: file, tbi: file, fields: '', prefix: '', vardb: file], ...]

    main:
    // Add empty tbi placeholder to input VCF channel
    ch_vcf_with_tbi = ch_vcf.map { meta, vcf -> [meta, vcf, []] }

    // Extract parallel lists from unified config
    ch_extracted = ch_db_configs.map { configs ->
        [
            configs.collect { it.vcf },
            configs.collect { it.tbi },
            configs.collect { it.vardb },
            configs.collect { [fields: it.fields ?: '', prefix: it.prefix ?: ''] }
        ]
    }

    // Annotate with all databases in one pass
    SNPSIFT_ANNMEM(
        ch_vcf_with_tbi,
        ch_extracted.map { it[0] },  // databases
        ch_extracted.map { it[1] },  // databases_tbi
        ch_extracted.map { it[2] },  // vardbs
        ch_extracted.map { it[3] }   // db_configs_simple
    )

    emit:
    vcf_tbi = SNPSIFT_ANNMEM.out.vcf  // channel: [val(meta), path(vcf), path(tbi)]
}
