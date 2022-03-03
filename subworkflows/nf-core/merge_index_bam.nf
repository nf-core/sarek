//
// MERGE INDEX BAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from '../../modules/local/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_BAM       } from '../../modules/nf-core/modules/samtools/merge/main'

workflow MERGE_INDEX_BAM {
    take:
        bam // channel: [mandatory] meta, bam

    main:
    ch_versions = Channel.empty()

    // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
    // Except if and only if skipping markduplicates or saving mapped bams

    // Figuring out if there is one or more bam(s) from the same sample
    // Is done for all, but not blocking
    bam.branch{
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{bam_to_merge}

    // Only if the when clause is true
    MERGE_BAM(bam_to_merge.multiple, [])
    INDEX_MERGE_BAM(bam_to_merge.single.mix(MERGE_BAM.out.bam))

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(INDEX_MERGE_BAM.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_BAM.out.versions.first())

    emit:
        bam_bai  = INDEX_MERGE_BAM.out.bam_bai
        versions = ch_versions
}
