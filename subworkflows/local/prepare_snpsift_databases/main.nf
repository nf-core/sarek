//
// Prepare SnpSift annotation databases
//

include { SNPSIFT_ANNMEM } from '../../../modules/nf-core/snpsift/annmem'

workflow PREPARE_SNPSIFT_DATABASES {
    take:
    val_db_configs  // List of maps: [[vcf: file, tbi: file, fields: '', prefix: '', vardb: null], ...]

    main:
    ch_configs = channel.fromList(val_db_configs)

    // Branch: create vardb if not provided
    ch_configs.branch {
        has_vardb: it.vardb != null
        needs_vardb: true
    }.set { ch_branched }

    // Create vardbs for databases that need them
    SNPSIFT_ANNMEM(
        ch_branched.needs_vardb.map { [[id: it.vcf.baseName], [], []] },
        ch_branched.needs_vardb.map { [it.vcf, it.tbi, [], it.fields ? [it.fields] : [[]], []] },
        true
    )

    // Join created vardbs back with their configs
    ch_created = SNPSIFT_ANNMEM.out.database
        .map { meta, vardb -> [meta.id, vardb] }
        .join(ch_branched.needs_vardb.map { [it.vcf.baseName, it] })
        .map { _id, vardb, config -> [config.vcf, config.tbi, vardb, config.fields ?: '', config.prefix ?: ''] }

    // Configs with pre-built vardb
    ch_prebuilt = ch_branched.has_vardb
        .map { [it.vcf, it.tbi, it.vardb, it.fields ?: '', it.prefix ?: ''] }

    // Collect all into output tuple
    ch_db_tuple = ch_prebuilt
        .mix(ch_created)
        .toList()
        .map { list ->
            [
                list.collect { it[0] },  // db_vcf
                list.collect { it[1] },  // db_vcf_tbi
                list.collect { it[2] },  // db_vardb
                list.collect { it[3] },  // db_fields
                list.collect { it[4] }   // db_prefixes
            ]
        }

    emit:
    db_tuple = ch_db_tuple
}
