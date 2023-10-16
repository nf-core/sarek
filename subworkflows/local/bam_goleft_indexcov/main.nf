//
// INDEXCOV : call large CNVs using goleft indexcov
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX_FOR_INDEXCOV               } from '../../../modules/local/bcftools/mpileup/main'


workflow BAM_GOLEFT_INDEXCOV {
    take:
    bam      // channel: [mandatory] [ meta, bam, bai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai     // channel: [mandatory] [ fasta ]

    main:
    versions = Channel.empty()


    SAMTOOLS_INDEX_FOR_INDEXCOV(bam, fasta, fai)


    //GOLEFT_INDEXCOV();

    emit:
    versions
}
