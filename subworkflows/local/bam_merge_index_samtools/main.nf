//
// MERGE INDEX BAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_BAM       } from '../../../modules/nf-core/samtools/merge/main'

workflow BAM_MERGE_INDEX_SAMTOOLS {
    take:
        bam // channel: [mandatory] meta, bam

    main:
    ch_versions = Channel.empty()

    // Figuring out if there is one or more bam(s) from the same sample
    bam.branch{
        //Here there actually is a list, so size() works
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{bam_to_merge}

    MERGE_BAM(bam_to_merge.multiple, [], [])
    INDEX_MERGE_BAM(bam_to_merge.single.mix(MERGE_BAM.out.bam))

    bam_bai = bam_to_merge.single.map{meta, bam -> [meta, bam[0]]}
        .mix(MERGE_BAM.out.bam)
        .join(INDEX_MERGE_BAM.out.bai)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(INDEX_MERGE_BAM.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BAM.out.versions.first())

    emit:
        bam_bai
        versions = ch_versions
}
