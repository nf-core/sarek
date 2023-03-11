/*
 * Sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { SAMTOOLS_VIEW_BAM  } from '../../modules/local/samtools_view_bam'
include { SAMTOOLS_SORT      } from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_SORT_INDEX } from '../../modules/local/samtools_sort_index'
include { BAM_STATS_SAMTOOLS } from '../../subworkflows/nf-core/bam_stats_samtools'

workflow BAM_SORT_INDEX_SAMTOOLS {
    take:
    ch_sam // channel: [ val(meta), [ bam ] ]
    ch_fasta

    main:
    /*
     * Sam to bam conversion
     */
    SAMTOOLS_VIEW_BAM  ( ch_sam )
    SAMTOOLS_SORT_INDEX ( SAMTOOLS_VIEW_BAM.out.bam )
    ch_sam
        .join( SAMTOOLS_SORT_INDEX.out.bam_bai )
        .map { it -> [ it[0], it[1], it[2], it[4], it[5] ] }
        .set { sortbam }
    BAM_STATS_SAMTOOLS ( SAMTOOLS_SORT_INDEX.out.bam_bai, ch_fasta )

    /*
     * SUBWORKFLOW: Create stats using samtools
     */
    BAM_STATS_SAMTOOLS.out.stats
        .join ( BAM_STATS_SAMTOOLS.out.idxstats )
        .join ( BAM_STATS_SAMTOOLS.out.flagstat )
        .map  { it -> [ it[1], it[2], it[3] ] }
        .set  { sortbam_stats_multiqc }
    samtools_versions = BAM_STATS_SAMTOOLS.out.versions

    emit:
    sortbam
    sortbam_stats_multiqc
    samtools_versions
}
