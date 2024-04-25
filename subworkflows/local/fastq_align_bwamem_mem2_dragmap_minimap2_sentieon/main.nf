//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'
include { MINIMAP2_ALIGN         } from '../../../modules/nf-core/minimap2/align/main'
include { SENTIEON_BWAMEM        } from '../../../modules/nf-core/sentieon/bwamem/main'

workflow FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_MINIMAP2_SENTIEON {
    take:
    reads // channel: [mandatory] meta, reads
    index // channel: [mandatory] index
    sort  // boolean: [mandatory] true -> sort, false -> don't sort
    fasta
    fasta_fai

    main:

    versions = Channel.empty()
    reports = Channel.empty()

    cigar_paf_format = false
    cigar_bam = false

    // Only one of the following should be run
    BWAMEM1_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem
    BWAMEM2_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem2
    DRAGMAP_ALIGN(reads, index, [[id:'no_fasta'], []], sort) // If aligner is dragmap
    MINIMAP2_ALIGN(reads, index, [[id:'no_fasta'], []], sort, cigar_paf_format, cigar_bam) // If aligner is minimap2
    // The sentieon-bwamem-module does sorting as part of the conversion from sam to bam.
    SENTIEON_BWAMEM(reads, index, fasta, fasta_fai) // If aligner is sentieon-bwamem

    // Get the bam files from the aligner
    // Only one aligner is run
    bam = Channel.empty()
    bam = bam.mix(BWAMEM1_MEM.out.bam)
    bam = bam.mix(BWAMEM2_MEM.out.bam)
    bam = bam.mix(DRAGMAP_ALIGN.out.bam)
    bam = bam.mix(MINIMAP2_ALIGN.out.bam)
    bam = bam.mix(SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bam ] })

    bai = SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bai ] }

    // Gather reports of all tools used
    reports = reports.mix(DRAGMAP_ALIGN.out.log)

    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_MEM.out.versions)
    versions = versions.mix(BWAMEM2_MEM.out.versions)
    versions = versions.mix(DRAGMAP_ALIGN.out.versions)
    versions = versions.mix(MINIMAP2_ALIGN.out.versions)
    versions = versions.mix(SENTIEON_BWAMEM.out.versions)

    emit:
    bam      // channel: [ [meta], bam ]
    bai      // channel: [ [meta], bai ]
    reports
    versions // channel: [ versions.yml ]
}
