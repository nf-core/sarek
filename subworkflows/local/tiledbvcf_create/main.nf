include { TILEDBVCF_CREATE } from '../../../modules/local/tiledbvcf/create/main'
include { TILEDBVCF_INGEST } from '../../../modules/local/tiledbvcf/ingest/main'

workflow TILEDBVCF_CREATE_DATASET {
    main:
    ch_versions = Channel.empty()

    if (params.tiledb_create_dataset) {
        TILEDBVCF_CREATE (
            params.tiledb_dataset_name
        )
        ch_versions = ch_versions.mix(TILEDBVCF_CREATE.out.versions)
        tiledb_db = TILEDBVCF_CREATE.out.tiledb_db
    } else {
        tiledb_db = Channel.of([id: params.tiledb_dataset_name], file(params.tiledb_dataset_name))
    }

    emit:
    tiledb_db   = tiledb_db   // channel: [ val(meta), path(tiledb_db) ]
    versions    = ch_versions       // channel: [ versions.yml ]
}

