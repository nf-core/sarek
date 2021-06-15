//include { SAMTOOLS_MERGE }               from '../../modules/nf-core/software/samtools/merge/main'   addParams(options: params.merge_bam_options)
include { QUALIMAP_BAMQC }               from '../../modules/nf-core/software/qualimap/bamqc/main'   addParams(options: params.qualimap_bamqc_options)
//include { SAMTOOLS_INDEX }               from '../../modules/nf-core/software/samtools/index/main'   addParams(options: params.samtools_index_options)
include { SAMTOOLS_STATS }               from '../../modules/nf-core/software/samtools/stats/main'   addParams(options: params.samtools_stats_options)

workflow MAPPING {
    take:
        skip_bamqc      // boolean: true/false
        skip_samtools   // boolean: true/false


    qualimap_bamqc = Channel.empty()
    samtools_stats = Channel.empty()

    // if (!skip_bamqc) {
    //     QUALIMAP_BAMQC(cram_markduplicates, target_bed, params.target_bed)
    //     qualimap_bamqc = QUALIMAP_BAMQC.out
    // }

    if (!skip_samtools) {
        SAMTOOLS_STATS(cram_markduplicates)
        samtools_stats = SAMTOOLS_STATS.out.stats
    }

    bam_reports = samtools_stats.mix(qualimap_bamqc)
