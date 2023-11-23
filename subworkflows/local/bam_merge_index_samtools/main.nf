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

    // Only when using intervals
    MERGE_BAM(bam_to_merge.multiple, [ [ id:'null' ], []], [ [ id:'null' ], []])

    // Mix intervals and no_intervals channels together
    bam_bai_merged = MERGE_BAM.out.bam.join(MERGE_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)

    // Index single bams, merged ones are indexed on the fly
    INDEX_MERGE_BAM(bam_to_merge.single)
    bam_bai_single = bam_to_merge.single.join(INDEX_MERGE_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)

    // Join with the bai file
    bam_bai = bam_bai_merged.mix(bam_bai_single)

    // Gather versions of all tools used
    versions = versions.mix(INDEX_MERGE_BAM.out.versions)
    versions = versions.mix(MERGE_BAM.out.versions)

    emit:
    bam_bai

    versions
}
