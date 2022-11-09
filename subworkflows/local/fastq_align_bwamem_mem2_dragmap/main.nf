//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'

workflow FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP {
    take:
        ch_reads     // channel: [mandatory] meta, reads
        ch_map_index // channel: [mandatory] mapping index
        sort         // boolean: [mandatory] true -> sort, false -> don't sort

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(ch_reads,   ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort) // If aligner is bwa-mem
    BWAMEM2_MEM(ch_reads,   ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort) // If aligner is bwa-mem2
    DRAGMAP_ALIGN(ch_reads, ch_map_index.map{ it -> [[id:it[0].baseName], it] }, sort) // If aligner is dragmap

    // Get the bam files from the aligner
    // Only one aligner is run
    ch_bam_mapped = Channel.empty()
    ch_bam_mapped = ch_bam_mapped.mix(BWAMEM1_MEM.out.bam)
    ch_bam_mapped = ch_bam_mapped.mix(BWAMEM2_MEM.out.bam)
    ch_bam_mapped = ch_bam_mapped.mix(DRAGMAP_ALIGN.out.bam)

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
