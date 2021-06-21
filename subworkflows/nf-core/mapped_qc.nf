include { QUALIMAP_BAMQC }               from '../../modules/nf-core/software/qualimap/bamqc/main'   addParams(options: params.qualimap_bamqc_options)
include { SAMTOOLS_STATS }               from '../../modules/nf-core/software/samtools/stats/main'   addParams(options: params.samtools_stats_options)

workflow MAPPING {
    take:
        skip_bamqc      // boolean: true/false
        skip_samtools   // boolean: true/false
        bam             //


    qualimap_bamqc = Channel.empty()
    samtools_stats = Channel.empty()

     //If skip_markduplicates then QC tools are run on mapped bams,
    //if !skip_markduplicates, then QC tools are run on duplicate marked bams
    //TODO: Need to check if it is ok like this, but saves one unneccary merging step. We could look at it, like QC before BQSR and after or so maybe
    qualimap_bamqc = Channel.empty()
    if(!skip_bamqc){
        QUALIMAP_BAMQC(bam_markduplicates, target_bed, params.target_bed)
        qualimap_bamqc = QUALIMAP_BAMQC.out
    }

    SAMTOOLS_BAM_TO_CRAM(bam_markduplicates, fasta)

    samtools_stats = Channel.empty()
    if (!skip_samtools) {
        SAMTOOLS_STATS(SAMTOOLS_BAM_TO_CRAM.out.cram)
        samtools_stats = SAMTOOLS_STATS.out.stats
    }

    qc_reports = samtools_stats.mix(qualimap_bamqc)
    qc_reports = report_markduplicates.mix(qc_reports)