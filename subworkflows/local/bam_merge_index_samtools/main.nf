//
// MERGE INDEX BAM
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAMTOOLS_INDEX as INDEX_MERGE_BAM } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_BAM       } from '../../../modules/nf-core/samtools/merge/main'
include { FGBIO_COPYUMIFROMREADNAME         } from '../../../modules/nf-core/fgbio/copyumifromreadname/main'

workflow BAM_MERGE_INDEX_SAMTOOLS {
    take:
    bam // channel: [mandatory] meta, bam

    main:
    versions = Channel.empty()

    // If UMIs started in read header, copy to RX tag
    if (params.umi_in_read_header ) {
        FGBIO_COPYUMIFROMREADNAME(bam)
        bam = FGBIO_COPYUMIFROMREADNAME.out.bam
    }

    // Figuring out if there is one or more bam(s) from the same sample
    bam_to_merge = bam.branch{ meta, bam_ ->
        // bam is a list, so use bam.size() to asses number of intervals
        single:   bam_.size() <= 1
            return [ meta, bam_[0] ]
        multiple: bam_.size() > 1
    }

    // Only when using intervals
    MERGE_BAM(bam_to_merge.multiple, [ [ id:'null' ], []], [ [ id:'null' ], []])

    // Mix intervals and no_intervals channels together
    bam_all = MERGE_BAM.out.bam.mix(bam_to_merge.single)

    // Index bam
    INDEX_MERGE_BAM(bam_all)

    // Join with the bai file
    bam_bai = bam_all.join(INDEX_MERGE_BAM.out.bai, failOnDuplicate: true, failOnMismatch: true)

    // Gather versions of all tools used
    versions = versions.mix(INDEX_MERGE_BAM.out.versions)
    versions = versions.mix(MERGE_BAM.out.versions)
    versions = versions.mix(FGBIO_COPYUMIFROMREADNAME.out.versions)

    emit:
    bam_bai

    versions
}
