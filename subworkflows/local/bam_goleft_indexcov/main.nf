//
// INDEXCOV : call large CNVs using goleft indexcov
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX_FOR_INDEXCOV  } from '../../../modules/local/samtools_index_for_indexcov/main'
include { GOLEFT_INDEXCOV              } from '../../../modules/local/goleft/indexcov/main'


workflow BAM_GOLEFT_INDEXCOV {
    take:
    bam      // channel: [mandatory] [ meta, bam, bai ]
    fasta     // channel: [mandatory] [ fasta ]
    fai     // channel: [mandatory] [ fai ]

    main:
    versions = Channel.empty()

    samtools_index = SAMTOOLS_INDEX_FOR_INDEXCOV(bam, fasta, fai)
    versions = versions.mix(samtools_index.versions.first())

    indexcov_ch = GOLEFT_INDEXCOV([id:"PREFIX"], samtools_index.output.flatten().collect(), fasta, fai)
    versions = versions.mix(indexcov_ch.versions)

    emit:
    versions
}

