include { TILEDBVCF_CREATE } from '../../../modules/local/tiledbvcf/create/main'
include { TILEDBVCF_INGEST } from '../../../modules/local/tiledbvcf/ingest/main'

workflow TILEDBVCF_CREATE {
    take:
    gvcf_files  // channel: [ val(meta), path(gvcf) ]

    main:
    ch_versions = Channel.empty()

    if (params.tiledb_create_store) {
        TILEDBVCF_CREATE (
            [id: params.tiledb_store_name], params.tiledb_store_name
        )
        ch_versions = ch_versions.mix(TILEDBVCF_CREATE.out.versions)
        tiledb_db = TILEDBVCF_CREATE.out.tiledb_db
    } else {
        tiledb_db = Channel.of([id: params.tiledb_store_name], file(params.tiledb_store_name))
    }

    if (params.tiledb_ingest_gvcfs) {
        TILEDBVCF_INGEST (
            tiledb_db,
            gvcf_files.collect()
        )
        ch_versions = ch_versions.mix(TILEDBVCF_INGEST.out.versions)
        final_tiledb_db = TILEDBVCF_INGEST.out.tiledb_db
    } else {
        final_tiledb_db = tiledb_db
    }

    emit:
    tiledb_db   = final_tiledb_db   // channel: [ val(meta), path(tiledb_db) ]
    versions    = ch_versions       // channel: [ versions.yml ]
}

