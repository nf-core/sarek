//
// Indexcov calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_REINDEX_BAM } from '../../../modules/local/samtools/reindex_bam'
include { GOLEFT_INDEXCOV      } from '../../../modules/nf-core/goleft/indexcov'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_INDEXCOV {
    take:
    cram      // channel: [mandatory] [ meta, cram, crai ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]

    main:
    versions = Channel.empty()

    // generate a cleaner bam index without duplicate, supplementary, etc. (Small workload because the bam itself is not re-generated)
    SAMTOOLS_REINDEX_BAM(
        cram,
        fasta,
        fasta_fai,
    )

    // create [ [id:directory], bams, bais ] for input of GOLEFT_INDEXCOV
    GOLEFT_INDEXCOV(
        SAMTOOLS_REINDEX_BAM.out.output.map { _meta, bam, bai -> [[id: "indexcov"], bam, bai] }.groupTuple(),
        fasta_fai,
    )

    out_indexcov = GOLEFT_INDEXCOV.out.output

    versions = versions.mix(GOLEFT_INDEXCOV.out.versions)
    versions = versions.mix(SAMTOOLS_REINDEX_BAM.out.versions)

    emit:
    out_indexcov
    versions
}
