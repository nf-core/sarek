//
// MERGE INDEX BAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX      as INDEX_CRAM } from '../../modules/local/samtools/index/main'
include { SAMTOOLS_MERGE_CRAM as MERGE_CRAM } from '../../modules/local/samtools/mergecram/main'

workflow MERGE_INDEX_CRAM {
    take:
        ch_cram       // channel: [mandatory] meta, cram
        fasta         // channel: [mandatory] fasta

    main:
    ch_versions = Channel.empty()

    // Figuring out if there is one or more cram(s) from the same sample
    cram_to_merge = ch_cram
        .map{ meta, cram ->
            meta.id = meta.sample
            def groupKey = groupKey(meta, meta.num_intervals)
            [meta, cram]
        }.groupTuple()
    .branch{
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }

    MERGE_CRAM(cram_to_merge.multiple, fasta)
    INDEX_CRAM(cram_to_merge.single.mix(MERGE_CRAM.out.cram))

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(INDEX_CRAM.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_CRAM.out.versions.first())

    emit:
        cram_crai = INDEX_CRAM.out.cram_crai
        versions  = ch_versions
}
