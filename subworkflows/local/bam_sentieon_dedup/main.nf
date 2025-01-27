//
// SENTIEON DEDUP

include { CRAM_QC_MOSDEPTH_SAMTOOLS } from '../cram_qc_mosdepth_samtools/main'
include { SENTIEON_DEDUP            } from '../../../modules/nf-core/sentieon/dedup/main'

workflow BAM_SENTIEON_DEDUP {
    take:
    bam                    // channel: [mandatory] [ meta, bam ]  // Although the channel is named "bam", it may contain cram-files.
    bai
    fasta                  // channel: [mandatory] [ fasta ]
    fasta_fai              // channel: [mandatory] [ fasta_fai ]
    intervals_bed_combined // channel: [optional]  [ intervals_bed ]

    main:
    versions = Channel.empty()
    reports  = Channel.empty()

    bam = bam.map{ meta, bam -> [ meta - meta.subMap('data_type'), bam ] }
    bai = bai.map{ meta, bai -> [ meta - meta.subMap('data_type'), bai ] }
    bam_bai = bam.join(bai, failOnMismatch:true, failOnDuplicate:true)
    SENTIEON_DEDUP(bam_bai, fasta, fasta_fai)

    // Join with the crai file
    cram = SENTIEON_DEDUP.out.cram.join(SENTIEON_DEDUP.out.crai, failOnDuplicate: true, failOnMismatch: true)

    // QC on CRAM
    CRAM_QC_MOSDEPTH_SAMTOOLS(cram, fasta, intervals_bed_combined)

    // Gather all reports generated
    reports = reports.mix(SENTIEON_DEDUP.out.metrics)
    reports = reports.mix(SENTIEON_DEDUP.out.metrics_multiqc_tsv)
    reports = reports.mix(SENTIEON_DEDUP.out.score)
    reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)

    // Gather versions of all tools used
    versions = versions.mix(SENTIEON_DEDUP.out.versions)
    versions = versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
    cram
    reports

    versions    // channel: [ versions.yml ]
}
