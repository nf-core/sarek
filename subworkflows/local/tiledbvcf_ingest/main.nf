include { TILEDBVCF_INGEST } from '../../../modules/local/tiledbvcf/ingest/main'

workflow TILEDBVCF_INGEST_VCF {
    take:
    vcf_ch   // channel: [ val(meta), path(vcf) ]
    tiledb_db   

    main:
    ch_versions = Channel.empty()

    TILEDBVCF_INGEST(
        vcf_ch, tiledb_db
    )

    ch_versions = ch_versions.mix(TILEDBVCF_INGEST.out.versions)

    emit:
    ingested_db = TILEDBVCF_INGEST.out.ingested_db  // channel: [ val(meta), path(ingested_db) ]
    versions    = ch_versions                       // channel: [ versions.yml ]
}
