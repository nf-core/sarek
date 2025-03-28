include { PARABRICKS_FQ2BAM        } from '../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { CRAM_SAMPLEQC            } from '../../../subworkflows/local/cram_sampleqc/main.nf'
include { CHANNEL_ALIGN_CREATE_CSV } from '../../../subworkflows/local/channel_align_create_csv/main'

workflow FASTQ_PREPROCESS_PARABRICKS {

    take:
    ch_reads                        // channel: [mandatory] meta, reads
    ch_fasta                        // channel: [mandatory] meta, fasta
    ch_index                        // channel: [mandatory] meta, index - bwa index
    ch_interval_file                // channel: [optional] meta, intervals_bed_combined
    ch_known_sites                  // channel [optional] known_sites_indels
    ch_ngscheckmate_bed             // channel: [mandatory] meta, ngscheckmate_bed
    ch_intervals_for_preprocessing  // channel: [optional] meta, intervals_for_preprocessing
    val_output_fmt                  // either bam or cram

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    PARABRICKS_FQ2BAM(
        ch_reads,           // channel: [ val(meta), reads ]
        ch_fasta,           // channel: [ val(meta), fasta ]
        ch_index,           // channel: [ val(meta), index ]
        ch_interval_file,   // channel: [ val(meta), interval_file ]
        ch_known_sites,     // channel: [ val(meta), known_sites ]
        val_output_fmt      // either bam or cram
    )

    ch_versions = ch_versions.mix(PARABRICKS_FQ2BAM.out.versions)

    cram_out = PARABRICKS_FQ2BAM.out.cram

    cram_out.view(tag: 'cram')

    CRAM_SAMPLEQC(
            cram_out,
            ch_ngscheckmate_bed,
            ch_fasta,
            params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'),
            ch_intervals_for_preprocessing
    )

    ch_versions = ch_versions.mix(CRAM_SAMPLEQC.out.versions)
    ch_reports = ch_reports.mix(CRAM_SAMPLEQC.out.reports.collect{ _meta, report -> [ report ] })

    cram_variant_calling =
        PARABRICKS_FQ2BAM.out.cram
        .join(PARABRICKS_FQ2BAM.out.crai, failOnDuplicate: true, failOnMismatch: true)

    CHANNEL_ALIGN_CREATE_CSV(
        cram_variant_calling,
        params.outdir,
        params.save_output_as_bam
    )

    emit:
    cram      = cram_variant_calling     // channel: [ val(meta), cram, crai ]
    versions  = ch_versions              // channel: [ versions.yml ]
    reports   = ch_reports
}

