params.merge_bam_options        = [:]
params.qualimap_bamqc_options   = [:]
params.samtools_stats_options   = [:]
params.samtools_view_options    = [:]
params.samtools_index_options   = [:]

include { SAMTOOLS_MERGE }                        from '../../modules/nf-core/software/samtools/merge/main'   addParams(options: params.merge_bam_options)
include { SAMTOOLS_INDEX }                        from '../../modules/nf-core/software/samtools/index/main'   addParams(options: params.samtools_index_options)
include { QUALIMAP_BAMQC }                        from '../../modules/nf-core/software/qualimap/bamqc/main'   addParams(options: params.qualimap_bamqc_options)
include { SAMTOOLS_STATS }                        from '../../modules/nf-core/software/samtools/stats/main'   addParams(options: params.samtools_stats_options)
include { SAMTOOLS_VIEW as SAMTOOLS_BAM_TO_CRAM } from '../../modules/nf-core/software/samtools/view/main.nf' addParams(options: params.samtools_view_options)


workflow CONVERT_BAM_TO_CRAM {
    take:
        skip_bamqc          // boolean: true/false
        skip_samtools       // boolean: true/false
        fasta               // channel: [mandatory] fasta
        bams                // channel: [mandatory] meta, bams
        target_bed          // channel: [optional]  target_bed

    main:

    //bams.dump(tag:'bams')
    if(params.skip_markduplicates){

        bams.branch{
            single:   it[1].size() == 1
            multiple: it[1].size() > 1
        }.set{ bam_to_merge }

        SAMTOOLS_MERGE(bam_to_merge.multiple)
        bam_merged = bam_to_merge.single.mix(SAMTOOLS_MERGE.out.merged_bam)

        SAMTOOLS_INDEX(bam_merged)
        bam_merged_index = bam_merged.join(SAMTOOLS_INDEX.out.bai)

    }else{
        bam_merged_index = bams
    }

    //If skip_markduplicates then QC tools are run on mapped bams,
    //if !skip_markduplicates, then QC tools are run on duplicate marked bams
    //TODO: Need to check if it is ok like this, but saves one unneccary merging step. We could look at it, like QC before BQSR and after or so maybe
    qualimap_bamqc = Channel.empty()
    if(!skip_bamqc){
        QUALIMAP_BAMQC(bam_merged_index, target_bed, params.target_bed)
        qualimap_bamqc = QUALIMAP_BAMQC.out
    }

    SAMTOOLS_BAM_TO_CRAM(bam_merged_index, fasta)

    samtools_stats = Channel.empty()
    if (!skip_samtools) {
        SAMTOOLS_STATS(SAMTOOLS_BAM_TO_CRAM.out.cram)
        samtools_stats = SAMTOOLS_STATS.out.stats
    }

    qc_reports = samtools_stats.mix(qualimap_bamqc)

    emit:
        cram = SAMTOOLS_BAM_TO_CRAM.out.cram // channel: [mandatory] meta, cram, crai
        qc   = qc_reports
}