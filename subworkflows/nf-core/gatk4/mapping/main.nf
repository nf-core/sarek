//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../../modules/nf-core/modules/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../../modules/nf-core/modules/bwa/mem/main'

workflow GATK4_MAPPING {
    take:
        ch_reads // channel: [mandatory] meta, reads
        ch_bwa   // channel: [mandatory] bwa

    main:

    ch_versions = Channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(ch_reads, ch_bwa, true) // If aligner is bwa-mem
    BWAMEM2_MEM(ch_reads, ch_bwa, true) // If aligner is bwa-mem2

    // Grouping the bams from the same samples not to stall the workflow
    ch_bam_mapped = BWAMEM1_MEM.out.bam.mix(BWAMEM2_MEM.out.bam).map{ meta, bam ->
        new_meta = meta.clone()
        // Removing unneeded fields in the new_meta map
        new_meta.remove('read_group') // read_group is now in the bam header
        new_meta.remove('size')       // size is only needed for the groupKey
        new_meta.id = meta.sample

        // groupKey is to makes sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        def groupKey = groupKey(new_meta, meta.numLanes * meta.size)

        //Returns the values we need
        tuple(groupKey, new_meta, bam)
    }.groupTuple(by:[0,1]).map{ groupKey, new_meta, bam -> [new_meta, bam] }

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    emit:
        bam      = ch_bam_mapped // channel: [ [meta], bam ]
        versions = ch_versions   // channel: [ versions.yml ]
}
