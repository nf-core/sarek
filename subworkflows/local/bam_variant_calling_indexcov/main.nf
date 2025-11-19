//
// Indexcov calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_REINDEX_BAM } from '../../../modules/local/samtools/reindex_bam/main'
include { GOLEFT_INDEXCOV      } from '../../../modules/nf-core/goleft/indexcov/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_INDEXCOV {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] [ meta, fasta_fai ]

    main:
    versions = Channel.empty()

    // generate a cleaner bam index without duplicate, supplementary, etc. (Small workload because the bam itself is not re-generated)
    reindex_ch = SAMTOOLS_REINDEX_BAM(
        cram,
        fasta,
        fasta_fai
        )

    versions = versions.mix(reindex_ch.versions)

    // create [ [id:directory], bams, bais ]
    indexcov_input_ch = reindex_ch.output.map{[[id:"indexcov"], it[1], it[2]]}.groupTuple()

    goleft_ch = GOLEFT_INDEXCOV(
        indexcov_input_ch,
        fasta_fai
        )

    versions = versions.mix(goleft_ch.versions)


    emit:

    out_indexcov = goleft_ch.output
    versions
}
