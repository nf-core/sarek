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
    bam_to_merge = bam.branch{ meta, bam ->
        // bam is a list, so use bam.size() to asses number of intervals
        single:   bam.size() <= 1
            return [ meta, bam[0] ]
        multiple: bam.size() > 1
    }

    bam_to_merge.multiple.view()
    // Only when using intervals
    MERGE_BAM(bam_to_merge.multiple, [], [])

    // Mix intervals and no_intervals channels together
    bam_all = MERGE_BAM.out.bam.mix(bam_to_merge.single)

    // Index bam
    INDEX_MERGE_BAM(bam_all)

    // Join with the bai file
    bam_bai = bam_all.join(INDEX_MERGE_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    versions = versions.mix(INDEX_MERGE_BAM.out.versions)
    versions = versions.mix(MERGE_BAM.out.versions)

    emit:
    bam_bai

    versions
}
