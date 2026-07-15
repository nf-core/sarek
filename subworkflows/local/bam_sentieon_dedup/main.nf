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

    bam = bam.map{ meta, bam_ -> [ meta - meta.subMap('data_type'), bam_ ] }
    bai = bai.map{ meta, bai_ -> [ meta - meta.subMap('data_type'), bai_ ] }
    bam_bai = bam.join(bai, failOnMismatch:true, failOnDuplicate:true)
    SENTIEON_DEDUP(bam_bai, fasta, fasta_fai)

    // Unified alignment output — BAM or CRAM depending on save_output_as_bam.
    // Branched rather than mixed because SENTIEON_DEDUP.out.bai is non-optional
    // (sentieon driver always writes a .bai), so joining .out.bam (optional, empty in
    // CRAM mode) with .out.bai would fail on the empty side.
    alignment = params.save_output_as_bam
        ? SENTIEON_DEDUP.out.bam.join(SENTIEON_DEDUP.out.bai, failOnDuplicate: true, failOnMismatch: true)
        : SENTIEON_DEDUP.out.cram.join(SENTIEON_DEDUP.out.crai, failOnDuplicate: true, failOnMismatch: true)

    // QC on alignment
    CRAM_QC_MOSDEPTH_SAMTOOLS(alignment, fasta, intervals_bed_combined)

    // Gather all reports generated
    reports = reports.mix(SENTIEON_DEDUP.out.metrics)
    reports = reports.mix(SENTIEON_DEDUP.out.metrics_multiqc_tsv)
    reports = reports.mix(SENTIEON_DEDUP.out.score)
    reports = reports.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.reports)

    // Gather versions of all tools used
    versions = versions.mix(CRAM_QC_MOSDEPTH_SAMTOOLS.out.versions)

    emit:
    alignment   // channel: [ meta, file, index ] — BAM or CRAM
    reports

    versions    // channel: [ versions.yml ]
}
