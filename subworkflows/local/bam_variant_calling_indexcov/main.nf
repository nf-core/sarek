//
// Indexcov calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_REINDEX_BAM } from '../../../modules/local/indexcov/main'
include { GOLEFT_INDEXCOV      } from '../../../modules/nf-core/goleft/indexcov/main'

// Seems to be the consensus on upstream modules implementation too
workflow BAM_VARIANT_CALLING_GERMLINE_INDEXCOV {
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

    goleft_ch = GOLEFT_INDEXCOV(
          reindex_ch.output.map{it[1]}.collect().combine(
              reindex_ch.output.map{it[2]}.collect()
              ).map{[[:], it[0], it[1] ]}
          ),
        fasta_fai
        )

    versions = versions.mix(goleft_ch.versions)


    emit:

    versions
}
