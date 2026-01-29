//
// Run SnpSift annotateMem to annotate VCF files with multiple databases
//

include { SNPSIFT_ANNMEM_CREATE_DB } from '../../../modules/local/snpsift/annmem/create_db/main'
include { SNPSIFT_ANNMEM          } from '../../../modules/local/snpsift/annmem/main'

workflow VCF_ANNOTATE_SNPSIFT {
    take:
    ch_vcf                  // channel: [val(meta), path(vcf)]
    val_db_configs          // List of maps: [[vcf: file, tbi: file, fields: '', prefix: '', vardb: null], ...]
    val_create_dbs          // Boolean: whether to create databases for entries without pre-built vardb

    main:
    // Add empty tbi placeholder to input VCF channel
    ch_vcf_with_tbi = ch_vcf.map { meta, vcf -> [meta, vcf, []] }

    // Extract components from unified config
    def databases = val_db_configs.collect { it.vcf }
    def databases_tbi = val_db_configs.collect { it.tbi }
    def db_configs_simple = val_db_configs.collect { [fields: it.fields ?: '', prefix: it.prefix ?: ''] }
    def prebuilt_vardbs = val_db_configs.collect { it.vardb }.findAll { it != null }

    if (val_create_dbs) {
        // Only create databases for entries WITHOUT pre-built vardb
        def dbs_to_create = val_db_configs
            .findAll { it.vardb == null }
            .collect { config ->
                [
                    [id: config.vcf.baseName],
                    config.vcf,
                    config.tbi,
                    config.fields ?: ''
                ]
            }

        if (dbs_to_create) {
            ch_dbs_to_create = Channel.fromList(dbs_to_create)

            SNPSIFT_ANNMEM_CREATE_DB(ch_dbs_to_create)

            // Combine created vardbs with pre-built ones
            ch_db_vardbs = SNPSIFT_ANNMEM_CREATE_DB.out.database
                .map { meta, vardb -> vardb }
                .collect()
                .map { created -> prebuilt_vardbs + created }
        } else {
            // All databases have pre-built vardbs
            ch_db_vardbs = Channel.value(prebuilt_vardbs)
        }
    } else {
        // Use pre-existing databases only
        // If user specified vardb paths in config, use those
        ch_db_vardbs = Channel.value(prebuilt_vardbs)
    }

    // Annotate with all databases in one pass
    SNPSIFT_ANNMEM(
        ch_vcf_with_tbi,
        databases,
        databases_tbi,
        ch_db_vardbs,
        db_configs_simple
    )

    emit:
    vcf_tbi = SNPSIFT_ANNMEM.out.vcf  // channel: [val(meta), path(vcf), path(tbi)]
}
