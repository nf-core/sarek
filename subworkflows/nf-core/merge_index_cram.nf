//
// MERGE INDEX CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_CRAM } from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_CRAM } from '../../modules/nf-core/modules/samtools/merge/main'

workflow MERGE_INDEX_CRAM {
    take:
        ch_cram       // channel: [mandatory] meta, cram
        fasta         // channel: [mandatory] fasta
        num_intervals

    main:
    ch_versions = Channel.empty()

    // Figuring out if there is one or more cram(s) from the same sample
    ch_cram = ch_cram
        .map{ meta, cram ->
            meta.id = meta.sample
            [meta, cram]
        }.groupTuple(size: num_intervals)
    .branch{
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{cram_to_merge}

    MERGE_CRAM(cram_to_merge.multiple, fasta)
    INDEX_CRAM(cram_to_merge.single.mix(MERGE_CRAM.out.cram))

    cram_crai = cram_to_merge.single
        .mix(MERGE_CRAM.out.cram)
        .join(INDEX_CRAM.out.crai)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(INDEX_CRAM.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_CRAM.out.versions.first())

    emit:
        cram_crai
        versions  = ch_versions
}
