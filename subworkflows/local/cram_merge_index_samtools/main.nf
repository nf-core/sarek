//
// MERGE INDEX CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_CRAM } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_CRAM } from '../../../modules/nf-core/samtools/merge/main'

workflow CRAM_MERGE_INDEX_SAMTOOLS {
    take:
    cram      // channel: [mandatory] meta, cram
    fasta     // channel: [mandatory] fasta
    fasta_fai // channel: [mandatory] fai for fasta

    main:
    versions = Channel.empty()

    // Figuring out if there is one or more cram(s) from the same sample
    cram_to_merge = cram.branch { meta, cram ->
        single: cram.size() <= 1
        return [meta, cram[0]]
        multiple: cram.size() > 1
    }

    // Only when using intervals
    MERGE_CRAM(cram_to_merge.multiple, fasta.map { it -> [[id: 'fasta'], it] }, fasta_fai.map { it -> [[id: 'fasta_fai'], it] })

    // Mix intervals and no_intervals channels together
    cram_all = MERGE_CRAM.out.cram.mix(cram_to_merge.single)

    // Index cram
    INDEX_CRAM(cram_all)

    // Join with the crai file
    cram_crai = cram_all.join(INDEX_CRAM.out.crai, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    versions = versions.mix(INDEX_CRAM.out.versions)
    versions = versions.mix(MERGE_CRAM.out.versions)

    emit:
    cram_crai
    versions
}
