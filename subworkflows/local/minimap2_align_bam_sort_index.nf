include { SAMTOOLS_SORT_INDEX } from '../../modules/local/samtools_sort_index'
include { SAMTOOLS_VIEW_BAM  } from '../../modules/local/samtools_view_bam'
include { MINIMAP2_ALIGN          } from '../../modules/local/minimap2/minimap2_align'


workflow MINIMAP2_ALIGN_BAM_SORT_INDEX {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index
        sort         // boolean: [mandatory] true -> sort, false -> don't sort

    main:
    MINIMAP2_ALIGN(ch_reads, ch_map_index) // If aligner is minimap2
    ch_align_sam = MINIMAP2_ALIGN.out.align_sam
    SAMTOOLS_VIEW_BAM  (ch_align_sam)
    SAMTOOLS_SORT_INDEX ( SAMTOOLS_VIEW_BAM.out.bam )
    // ch_align_sam
    //     .join( SAMTOOLS_SORT_INDEX.out.bam_bai )
    //     .map { it -> [ it[0], it[1], it[2], it[4], it[5] ] }
    //     .set { bam }
    
    ch_versions = MINIMAP2_ALIGN.out.versions.first()

    emit:
    bam = SAMTOOLS_SORT_INDEX.out.bam
    versions = ch_versions
}