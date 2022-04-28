//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../../modules/nf-core/modules/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../../modules/nf-core/modules/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../../modules/nf-core/modules/dragmap/align/main'

workflow GATK4_MAPPING {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index
        sort         // boolean: [mandatory] true -> sort, false -> don't sort

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(ch_reads,   ch_map_index, sort) // If aligner is bwa-mem
    BWAMEM2_MEM(ch_reads,   ch_map_index, sort) // If aligner is bwa-mem2
    DRAGMAP_ALIGN(ch_reads, ch_map_index, sort) // If aligner is dragmap

    // Grouping the bams from the same samples not to stall the workflow
    ch_bam_mapped = BWAMEM1_MEM.out.bam.mix(BWAMEM2_MEM.out.bam, DRAGMAP_ALIGN.out.bam).map{ meta, bam ->
        new_meta = meta.clone()

        numLanes = meta.numLanes ?: 1
        size     = meta.size     ?: 1

        // Use groupKey to make sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all reads from all channels are mapped
        def groupKey = groupKey(new_meta, numLanes * size)

        //Returns the values we need
        tuple(groupKey, new_meta, bam)
    }.groupTuple(by:[0,1]).map{ groupKey, new_meta, bam -> [new_meta, bam] }

    ch_bam_mapped.view()

    // Gather reports of all tools used
    ch_reports = ch_reports.mix(DRAGMAP_ALIGN.out.log)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())
    ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions.first())

    emit:
        bam      = ch_bam_mapped // channel: [ [meta], bam ]
        reports  = ch_reports
        versions = ch_versions   // channel: [ versions.yml ]
}
