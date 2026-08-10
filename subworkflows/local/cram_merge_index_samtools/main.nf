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
    fasta     // channel: [mandatory] meta, fasta
    fasta_fai // channel: [mandatory] meta, fai

    main:

    // Figuring out if there is one or more cram(s) from the same sample
    cram_to_merge = cram.branch { meta, cram_files ->
        single: cram_files.size() <= 1
        return [meta, cram_files[0]]
        multiple: cram_files.size() > 1
    }

    // Only when using intervals
    MERGE_CRAM(cram_to_merge.multiple.map { meta, crams -> [ meta, crams, [] ] }, fasta.combine(fasta_fai).map { meta, fasta_, _meta_fai, fai -> [ meta, fasta_, fai, [] ] }.collect())

    // Mix intervals and no_intervals channels together
    cram_crai_merged = MERGE_CRAM.out.cram.join(MERGE_CRAM.out.index)

    // Index cram, multiple ones are indexed on the fly
    INDEX_CRAM(cram_to_merge.single)
    cram_crai_single = cram_to_merge.single.join(INDEX_CRAM.out.index)

    // Mix intervals and no_intervals channels together
    cram_crai = cram_crai_merged.mix(cram_crai_single)

    emit:
    cram_crai
}
