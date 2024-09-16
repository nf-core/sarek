include { TILEDBVCF_CREATE } from '../../../modules/local/tiledbvcf/create/main'

workflow TILEDBVCF_CREATE_DATASET {
    main:
    ch_versions = Channel.empty()

    if (params.tiledb_create_dataset) {
        // Create a channel that emits the required tuple
        ch_tiledb_input = Channel.of([meta: [id: 'sample1', other_info: 'info'], db_name: params.tiledb_dataset_name])

        TILEDBVCF_CREATE(ch_tiledb_input)
        ch_versions = ch_versions.mix(TILEDBVCF_CREATE.out.versions)
        tiledb_db = TILEDBVCF_CREATE.out.tiledb_db
    } else {
        tiledb_db = Channel.of([id: params.tiledb_dataset_name], val(params.tiledb_dataset_name))
    }

    emit:
    tiledb_db   = tiledb_db   // channel: [ val(meta), path(tiledb_db) ]
    versions    = ch_versions       // channel: [ versions.yml ]
}


