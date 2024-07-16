include { LOFREQ_SOMATIC } from '../../../modules/nf-core/lofreq/somatic/main.nf'

workflow BAM_VARIANT_CALLING_SOMATIC_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals ]

    main:
    versions = Channel.empty()

    intervals_ch = intervals.map{ bed, num -> [bed]}

    input_intervals = input.combine(intervals_ch)
        // Sort channel elements for LOFREQ_SOMATIC module
        .map {meta, normal_bam, normal_bai, tumor_bam, tumor_bai, intervals -> [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, intervals]}
    
    LOFREQ_SOMATIC(input_intervals, fasta, fai)
    
    vcf = Channel.empty().mix(LOFREQ_SOMATIC.out.vcf) 
    versions = versions.mix(LOFREQ_SOMATIC.out.versions)

    emit:
    vcf
    versions
}
