/*
 * Run SAMtools stats, flagstat and idxstats
 */

include { SAMTOOLS_STATS    } from '../../modules/nf-core/modules/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/modules/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/modules/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [bai] ]
    ch_id_fasta

    main:
    /*
     * Stats with samtools
     */
    ch_id_fasta
       .map{ it -> it[1] }
       .set{ ch_fasta }
    SAMTOOLS_STATS    ( ch_bam_bai, ch_fasta )

    /*
     * Flagstat with samtools
     */
    SAMTOOLS_FLAGSTAT ( ch_bam_bai )

    /*
     * Idxstats with samtools
     */
    SAMTOOLS_IDXSTATS ( ch_bam_bai )

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions  = SAMTOOLS_STATS.out.versions     //    path: version.yml
}
