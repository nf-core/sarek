include { TILEDBVCF_INGEST } from '../../../modules/local/tiledbvcf/ingest/main'

workflow TILEDBVCF_INGEST_VCF {
    take:
    vcf_files   // channel: [ val(meta), path(vcf) ]
    tiledb_db    // channel: [ val(meta), path(tiledb_db) ]

    main:
    ch_versions = Channel.empty()

    // Group all vcf files into a single list
    vcf_files_grouped = vcf_files
        .map { meta, vcf -> vcf }
        .collect()

    TILEDBVCF_INGEST(
        vcf_files_grouped
    )

    ch_versions = ch_versions.mix(TILEDBVCF_INGEST.out.versions)

    emit:
    ingested_db = TILEDBVCF_INGEST.out.ingested_db  // channel: [ val(meta), path(ingested_db) ]
    versions    = ch_versions                       // channel: [ versions.yml ]
}
