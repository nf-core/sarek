//
// Run SnpSift annotateMem to annotate VCF files with multiple databases
//

include { SNPSIFT_ANNMEM_CREATE_DB } from '../../../modules/local/snpsift/annmem/create_db/main'
include { SNPSIFT_ANNMEM          } from '../../../modules/nf-core/snpsift/annmem/main'

workflow VCF_ANNOTATE_SNPSIFT {
    take:
    ch_vcf                  // channel: [val(meta), path(vcf), path(tbi)]
    val_databases           // List: [path(db1.vcf.gz), path(db2.vcf.gz), ...]
    val_databases_tbi       // List: [path(db1.vcf.gz.tbi), path(db2.vcf.gz.tbi), ...]
    val_db_configs          // List: [[fields: 'ID', prefix: ''], [fields: 'AF', prefix: 'ExAC_'], ...]
    val_create_dbs          // Boolean: whether to create databases (false = use existing)

    main:
    ch_versions = Channel.empty()
    ch_db_vardbs = Channel.empty()

    // Step 1: Create databases if requested
    if (val_create_dbs) {
        // Create channel with database info
        ch_dbs_to_create = Channel.fromList(
            val_databases.withIndex().collect { db, idx ->
                [
                    [id: db.baseName],
                    db,
                    val_databases_tbi[idx],
                    val_db_configs[idx].fields ?: ''
                ]
            }
        )

        SNPSIFT_ANNMEM_CREATE_DB(ch_dbs_to_create)

        ch_db_vardbs = SNPSIFT_ANNMEM_CREATE_DB.out.database
            .map { meta, vardb -> vardb }
            .collect()

        ch_versions = ch_versions.mix(SNPSIFT_ANNMEM_CREATE_DB.out.versions.first())
    } else {
        // Use pre-existing databases
        // SnpSift will look for .snpsift.vardb directories next to VCF files
        ch_db_vardbs = Channel.value([])
    }

    // Step 2: Annotate with all databases in one pass
    SNPSIFT_ANNMEM(
        ch_vcf,
        val_databases,
        val_databases_tbi,
        ch_db_vardbs,
        val_db_configs
    )

    ch_versions = ch_versions.mix(SNPSIFT_ANNMEM.out.versions)

    emit:
    vcf_tbi  = SNPSIFT_ANNMEM.out.vcf  // channel: [val(meta), path(vcf), path(tbi)]
    versions = ch_versions             // channel: path(versions.yml)
}
