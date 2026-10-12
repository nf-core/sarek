//
// Prepare SnpSift annotation databases
//

include { SNPSIFT_ANNMEMCREATE } from '../../../modules/nf-core/snpsift/annmemcreate'

workflow PREPARE_SNPSIFT_DATABASES {
    take:
    val_db_configs  // List of maps: [[vcf: file, tbi: file, fields: '', prefix: '', vardb: null], ...]

    main:
    ch_configs = channel.fromList(val_db_configs)

    // Branch: create vardb if not provided
    ch_branched = ch_configs.branch { config ->
        has_vardb: config.vardb != null
        needs_vardb: true
    }

    // Create vardbs for databases that need them
    // Convert semicolon-separated fields to comma-separated (SnpSift expects commas)
    SNPSIFT_ANNMEMCREATE(
        ch_branched.needs_vardb.map { config -> [[id: config.vcf.baseName], config.vcf, config.tbi, config.fields ? config.fields.replace(';', ',') : ''] }
    )

    // Join created vardbs back with their configs
    ch_created = SNPSIFT_ANNMEMCREATE.out.database
        .map { meta, vardb -> [meta.id, vardb] }
        .join(ch_branched.needs_vardb.map { config -> [config.vcf.baseName, config] })
        .map { _id, vardb, config -> [config.vcf, config.tbi, vardb, config.fields ? config.fields.replace(';', ',') : '', config.prefix ?: ''] }

    // Configs with pre-built vardb
    ch_prebuilt = ch_branched.has_vardb
        .map { config -> [config.vcf, config.tbi, config.vardb, config.fields ? config.fields.replace(';', ',') : '', config.prefix ?: ''] }

    // Collect all into output tuple: [ db_vcf, db_vcf_tbi, db_vardb, db_fields, db_prefixes ]
    // Each row is [ vcf, tbi, vardb, fields, prefix ], so transposing regroups them by column.
    // The caller only runs this subworkflow when there is at least one config, so the list is
    // never empty here (transpose collapses an empty list to [] rather than five empty lists).
    ch_db_tuple = ch_prebuilt
        .mix(ch_created)
        .toList()
        .map { rows -> rows.transpose() }

    emit:
    db_tuple = ch_db_tuple
}
