include { TILEDBVCF_INGEST } from '../../../modules/local/tiledbvcf/ingest/main'
include { TABIX_TABIX as TABIX_TILEDBVCF_INGEST                   } from '../../../modules/nf-core/tabix/tabix/main'

workflow TILEDBVCF_INGEST_VCF {
    take:
    vcf_ch   // channel: [ val(meta), path(vcf) ]
    tiledb_db   

    main:
    ch_versions = Channel.empty()

    TABIX_TILEDBVCF_INGEST(vcf_ch)
    
    //create a channel with the tbi file
    tbi_ch = TABIX_TILEDBVCF_INGEST.out.tbi
    
    TILEDBVCF_INGEST(
        vcf_ch, tbi_ch, tiledb_db
    )

    ch_versions = ch_versions.mix(TILEDBVCF_INGEST.out.versions)

    emit:
    versions    = ch_versions                       // channel: [ versions.yml ]
}
