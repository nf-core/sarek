//
// Prepare SnpSift annotation databases
//

include { SNPSIFT_ANNMEM } from '../../../modules/local/snpsift/annmem/main'

workflow PREPARE_SNPSIFT_DATABASES {
    take:
    val_db_configs  // List of maps: [[vcf: file, tbi: file, fields: '', prefix: '', vardb: null], ...]

    main:
    // Create indexed channel preserving original order
    ch_indexed_configs = Channel.fromList(
        val_db_configs.withIndex().collect { config, idx -> [idx, config] }
    )

    // Branch by whether vardb exists
    ch_indexed_configs.branch {
        has_vardb: it[1].vardb != null
        needs_vardb: it[1].vardb == null
    }.set { ch_branched }

    // Create databases for those that need it
    // First tuple: empty sample VCF (not used in create mode)
    // Second tuple: database VCF and fields to index
    ch_to_create = ch_branched.needs_vardb.map { idx, config ->
        [
            [[id: config.vcf.baseName, idx: idx], [], []],  // Empty sample VCF tuple
            [config.vcf, config.tbi, [], config.fields ? [config.fields] : [[]], []]  // Database info
        ]
    }

    // Call unified module in create mode
    SNPSIFT_ANNMEM(
        ch_to_create.map { it[0] },  // Empty sample tuple
        ch_to_create.map { it[1] },  // Database VCF to index
        true                         // create = true for database creation
    )

    // Join created vardb back with original config to build complete tuple
    // Output: [idx, vcf, tbi, vardb, fields, prefix]
    ch_created_complete = SNPSIFT_ANNMEM.out.database
        .map { meta, vardb -> [meta.idx, vardb] }
        .join(ch_branched.needs_vardb)
        .map { idx, vardb, config ->
            [idx, config.vcf, config.tbi, vardb, config.fields ?: '', config.prefix ?: '']
        }

    // Pre-built configs: [idx, vcf, tbi, vardb, fields, prefix]
    ch_prebuilt_complete = ch_branched.has_vardb
        .map { idx, config ->
            [idx, config.vcf, config.tbi, config.vardb, config.fields ?: '', config.prefix ?: '']
        }

    // Merge, sort by index, and build final tuple for SNPSIFT_ANNMEM
    ch_db_tuple = ch_prebuilt_complete
        .mix(ch_created_complete)
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { list ->
            [
                list.collect { it[1] },  // db_vcf
                list.collect { it[2] },  // db_vcf_tbi
                list.collect { it[3] },  // db_vardb
                list.collect { it[4] },  // db_fields
                list.collect { it[5] }   // db_prefixes
            ]
        }

    emit:
    db_tuple = ch_db_tuple  // channel: [[db_vcf], [db_vcf_tbi], [db_vardb], [db_fields], [db_prefixes]]
}
