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
    versions = Channel.empty()

    // Figuring out if there is one or more bam(s) from the same sample
    bam_to_merge = bam.branch{
        // Here there actually is a list, so size() works
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }

    // Only when using intervals
    MERGE_BAM(bam_to_merge.multiple, [], [])

    // Mix intervals and no_intervals channels together
    INDEX_MERGE_BAM(MERGE_BAM.out.bam.mix(bam_to_merge.single))

    // Mix intervals and no_intervals channels together
    bam_bai = MERGE_BAM.out.bam.mix(bam_to_merge.single.map{ meta, bam -> [ meta, bam[0] ] })
        // Join with the bai file
        .join(INDEX_MERGE_BAM.out.bai)

    // Gather versions of all tools used
    versions = versions.mix(INDEX_MERGE_BAM.out.versions)
    versions = versions.mix(MERGE_BAM.out.versions)

    emit:
    bam_bai

    versions
}
