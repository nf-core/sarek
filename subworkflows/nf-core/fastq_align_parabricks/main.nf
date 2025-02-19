//
// Alignment and BQSR with Nvidia CLARA Parabricks
//
include { PARABRICKS_FQ2BAM    } from '../../../modules/nf-core/parabricks/fq2bam/main'
include { PARABRICKS_APPLYBQSR } from '../../../modules/nf-core/parabricks/applybqsr/main'

workflow FASTQ_ALIGN_PARABRICKS {

    take:
    ch_reads // channel: [mandatory] meta, reads
    ch_fasta // channel: [mandatory] meta, fasta
    ch_index // channel: [mandatory] meta, index
    ch_interval_file // channel: [optional] meta, intervals_bed_combined
    ch_known_sites // channel [optional] known_sites_indels

    main:
    ch_versions          = Channel.empty()
    ch_bam               = Channel.empty()
    ch_bai               = Channel.empty()
    ch_bqsr_table        = Channel.empty()
    ch_qc_metrics        = Channel.empty()
    ch_duplicate_metrics = Channel.empty()

    PARABRICKS_FQ2BAM(
        ch_reads,
        ch_fasta,
        ch_index,
        ch_interval_file,
        ch_known_sites
    )

    // Collecting FQ2BAM outputs
    ch_bam               = ch_bam.mix(PARABRICKS_FQ2BAM.out.bam)
    ch_bai               = ch_bai.mix(PARABRICKS_FQ2BAM.out.bai)
    ch_qc_metrics        = ch_qc_metrics.mix(PARABRICKS_FQ2BAM.out.qc_metrics)
    ch_bqsr_table        = ch_bqsr_table.mix(PARABRICKS_FQ2BAM.out.bqsr_table)
    ch_duplicate_metrics = ch_duplicate_metrics.mix(PARABRICKS_FQ2BAM.out.duplicate_metrics)
    ch_versions          = ch_versions.mix(PARABRICKS_FQ2BAM.out.versions)

    // Apply BQSR
    PARABRICKS_APPLYBQSR(
        ch_bam,
        ch_bai,
        ch_bqsr_table.ifEmpty([]),
        ch_interval_file,
        ch_fasta
    )
    ch_versions = ch_versions.mix(PARABRICKS_APPLYBQSR.out.versions)

    emit:
    bam               = PARABRICKS_APPLYBQSR.out.bam // channel: [ [meta], bam ]
    bai               = PARABRICKS_APPLYBQSR.out.bai // channel: [ [meta], bai ]
    qc_metrics        = ch_qc_metrics                // channel: [ [meta], qc_metrics ]
    duplicate_metrics = ch_duplicate_metrics         // channel: [ [meta], duplicate_metrics ]
    bqsr_table        = ch_bqsr_table                // channel: [ [meta], bqsr_table ]
    versions          = ch_versions                  // channel: [ versions.yml ]
}
