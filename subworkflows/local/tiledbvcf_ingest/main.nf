include { TILEDBVCF_INGEST } from '../../../modules/local/tiledbvcf/ingest/main'

workflow TILEDBVCF_INGEST_WORKFLOW {
    take:
    gvcf_files   // channel: [ val(meta), path(gvcf) ]
    tiledb_db    // channel: [ val(meta), path(tiledb_db) ]

    main:
    ch_versions = Channel.empty()

    // Group all gVCF files into a single list
    gvcf_files_grouped = gvcf_files
        .map { meta, gvcf -> gvcf }
        .collect()

    // Combine the grouped gVCF files with the TileDB database
    ch_input = tiledb_db.combine(gvcf_files_grouped)

    TILEDBVCF_INGEST(
        ch_input
    )

    ch_versions = ch_versions.mix(TILEDBVCF_INGEST.out.versions)

    emit:
    ingested_db = TILEDBVCF_INGEST.out.ingested_db  // channel: [ val(meta), path(ingested_db) ]
    versions    = ch_versions                       // channel: [ versions.yml ]
}
