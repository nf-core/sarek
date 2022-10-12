//
// MERGE INDEX CRAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_CRAM } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_CRAM } from '../../../modules/nf-core/samtools/merge/main'

workflow CRAM_MERGE_INDEX_SAMTOOLS {
    take:
        ch_cram       // channel: [mandatory] meta, cram
        fasta         // channel: [mandatory] fasta
        fasta_fai     // channel: [mandatory] fai for fasta

    main:
    ch_versions = Channel.empty()

    // Figuring out if there is one or more cram(s) from the same sample
    ch_cram_to_merge = ch_cram.map{ meta, cram ->

        [groupKey([
                    data_type:      meta.data_type,
                    id:             meta.sample,
                    num_intervals:  meta.num_intervals,
                    patient:        meta.patient,
                    sample:         meta.sample,
                    sex:            meta.sex,
                    status:         meta.status,
                    ],
                meta.num_intervals),
        cram]
    }.groupTuple()
    .branch{
        //Warning: size() calculates file size not list length here, so use num_intervals instead
        single:   it[0].num_intervals <= 1
        multiple: it[0].num_intervals > 1
    }

    MERGE_CRAM(ch_cram_to_merge.multiple, fasta, fasta_fai)
    INDEX_CRAM(ch_cram_to_merge.single.mix(MERGE_CRAM.out.cram))

    cram_crai = ch_cram_to_merge.single.map{meta, cram -> [meta, cram[0]]}
        .mix(MERGE_CRAM.out.cram)
        .join(INDEX_CRAM.out.crai)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(INDEX_CRAM.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_CRAM.out.versions.first())

    emit:
        cram_crai
        versions  = ch_versions
}
