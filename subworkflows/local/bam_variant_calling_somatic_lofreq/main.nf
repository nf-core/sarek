include { LOFREQ_SOMATIC } from '../../../modules/nf-core/lofreq/somatic/main '

workflow BAM_VARIANT_CALLING_SOMATIC_LOFREQ {
    take:
    input     // channel: [mandatory] [ meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai       // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals ]

    main:
    versions = Channel.empty()

    
    input_intervals = input.combine(intervals)
        // Sort channel elements for LOFREQ_SOMATIC module
        .map {meta, normal_bam, normal_bai, tumor_bam, tumor_bai, intervals -> [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, intervals]}

    
    fasta_ch = fasta ? fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] } : [[id: 'null'], []]  


    fasta_val_ch = fasta_ch.first() 


    fai_ch = fai ? fai.map{ fai -> [ [ id:fai.baseName ], fai ] } : [[id: 'null'], []]

   
    fai_val_ch = fai_ch.first()

    
    LOFREQ_SOMATIC(input_intervals, fasta_val_ch, fai_val_ch)


    
    vcf = Channel.empty().mix(LOFREQ_SOMATIC.out.vcf) 
    versions = versions.mix(LOFREQ_SOMATIC.out.versions)
    
    emit:
    vcf
    versions
}
