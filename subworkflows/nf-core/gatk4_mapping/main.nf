//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM                     } from '../../../modules/nf-core/modules/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM          } from '../../../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_INDEX as INDEX_MAPPING } from '../../../modules/local/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_MAPPING } from '../../../modules/nf-core/modules/samtools/merge/main'

workflow GATK4_MAPPING {
    take:
        reads     // channel: [mandatory] meta, reads
        bwa       // channel: [mandatory] bwa
        fasta     // channel: [mandatory] fasta
        fasta_fai // channel: [mandatory] fasta_fai

    main:
    ch_versions = Channel.empty()

    // Only one of the following will be run
    BWAMEM1_MEM(reads, bwa, true) // If aligner is bwa-mem
    BWAMEM2_MEM(reads, bwa, true) // If aligner is bwa-mem2

    // Grouping the bams from the same samples not to stall the workflow
    // Removing unneeded fields in the new_meta map
    bam_mapped = BWAMEM1_MEM.out.bam.mix(BWAMEM2_MEM.out.bam).map{ meta, bam ->
        new_meta = meta.clone()
        new_meta.remove('read_group')
        new_meta.remove('size')
        new_meta.id = meta.sample

        // groupKey is to makes sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        def groupKey = groupKey(new_meta, meta.numLanes * meta.size)

        //Returns the values we need
        tuple(groupKey, new_meta, bam)
    }.groupTuple(by:[0,1]).map{ groupKey, new_meta, bam -> [new_meta, bam] }

    // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
    // Except if and only if skipping markduplicates or saving mapped bams

    // Figuring out if there is one or more bam(s) from the same sample
    // Is done for all, but not blocking
    bam_mapped.branch{
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{bam_to_merge}

    // Only if the when clause is true
    MERGE_MAPPING(bam_to_merge.multiple, [])
    INDEX_MAPPING(bam_to_merge.single.mix(MERGE_MAPPING.out.bam))

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    ch_versions = ch_versions.mix(INDEX_MAPPING.out.versions.first())
    ch_versions = ch_versions.mix(MERGE_MAPPING.out.versions.first())

    emit:
        bam         = bam_mapped
        bam_indexed = INDEX_MAPPING.out.bam_bai
        versions    = ch_versions
}
