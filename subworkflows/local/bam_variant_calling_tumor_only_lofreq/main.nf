include { LOFREQ_CALLPARALLEL } from '../../../modules/nf-core/lofreq/callparallel/main.nf'

workflow BAM_VARIANT_CALLING_TUMOR_ONLY_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, tumor_bam, tumor_bai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals ]

    main:
    versions = Channel.empty()

    intervals_ch = intervals.map{ bed, num -> [bed]}

    input_intervals = input.combine(intervals_ch)
        // Sort channel elements for LOFREQ_SOMATIC module
        .map {meta, tumor_bam, tumor_bai, intervals -> [meta, tumor_bam, tumor_bai, intervals]}

    LOFREQ_CALLPARALLEL(input_intervals, fasta, fai)

    
    vcf = Channel.empty().mix(LOFREQ_CALLPARALLEL.out.vcf) 
    versions = versions.mix(LOFREQ_CALLPARALLEL.out.versions)
    
    emit:
    vcf
    versions
}
