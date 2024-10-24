//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'
include { SENTIEON_BWAMEM        } from '../../../modules/nf-core/sentieon/bwamem/main'
include { PARABRICKS_FQ2BAM      } from '../../../modules/nf-core/parabricks/fq2bam/main'

workflow FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON {
    take:
    reads // channel: [mandatory] meta, reads
    // interval_file // channel: [optional for parabricks] intervals_bed_combined
    index // channel: [mandatory] index
    sort  // boolean: [mandatory] true -> sort, false -> don't sort
    fasta
    fasta_fai
    // known_sites // channel [optional for parabricks] known_sites_indels

    main:
    versions = Channel.empty()
    reports = Channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem
    BWAMEM2_MEM(reads, index, [[id:'no_fasta'], []], sort) // If aligner is bwa-mem2
    DRAGMAP_ALIGN(reads, index, [[id:'no_fasta'], []], sort) // If aligner is dragmap
    // The sentieon-bwamem-module does sorting as part of the conversion from sam to bam.
    SENTIEON_BWAMEM(reads, index, fasta, fasta_fai) // If aligner is sentieon-bwamem
    // The parabricks-fq2bam module performs alignment and sorting as part of the conversion from fastq to bam.
    // Additionally, it can perform mark duplicates and generate a bqsr table (input for parabricks-applybqsr module)
    PARABRICKS_FQ2BAM(reads.map{meta, reads -> [meta, reads, []] } , fasta, index, [])

    // Get the bam files from the aligner
    // Only one aligner is run
    bam = Channel.empty()
    bam = bam.mix(BWAMEM1_MEM.out.bam)
    bam = bam.mix(BWAMEM2_MEM.out.bam)
    bam = bam.mix(DRAGMAP_ALIGN.out.bam)
    bam = bam.mix(SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bam ] })
    bam = bam.mix(PARABRICKS_FQ2BAM.out.bam)

    bai = SENTIEON_BWAMEM.out.bam_and_bai.map{ meta, bam, bai -> [ meta, bai ] }
    bai = bai.mix(PARABRICKS_FQ2BAM.out.bai)

    // Optional outputs for parabricks_fq2bam
    qc_metrics = PARABRICKS_FQ2BAM.out.qc_metrics
    bqsr_table = PARABRICKS_FQ2BAM.out.bqsr_table
    duplicate_metrics = PARABRICKS_FQ2BAM.out.duplicate_metrics

    // Gather reports of all tools used
    reports = reports.mix(DRAGMAP_ALIGN.out.log)

    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_MEM.out.versions)
    versions = versions.mix(BWAMEM2_MEM.out.versions)
    versions = versions.mix(DRAGMAP_ALIGN.out.versions)
    versions = versions.mix(SENTIEON_BWAMEM.out.versions)
    versions = versions.mix(PARABRICKS_FQ2BAM.out.versions)

    emit:
    bam      // channel: [ [meta], bam ]
    bai      // channel: [ [meta], bai ]
    qc_metrics // file: qc_metrics/alignment.txt
    bqsr_table // file: *.table
    duplicate_metrics // file: duplicate-metrics.txt
    reports
    versions // channel: [ versions.yml ]
}
