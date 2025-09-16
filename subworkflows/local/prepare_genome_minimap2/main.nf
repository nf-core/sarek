//
// Prepare genome files for MINIMAP2
//

include { MINIMAP2_INDEX } from '../../../modules/local/minimap2/index/main'

workflow PREPARE_GENOME_MINIMAP2 {
    take:
    ch_fasta  // channel: [ val(meta), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    // Build minimap2 index
    MINIMAP2_INDEX(ch_fasta)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    emit:
    index    = MINIMAP2_INDEX.out.index
    versions = ch_versions
}