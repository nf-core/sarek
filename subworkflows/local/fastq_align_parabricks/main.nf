include { PARABRICKS_FQ2BAM        } from '../../../modules/nf-core/parabricks/fq2bam/main.nf'
include { CRAM_SAMPLEQC            } from '../../../subworkflows/local/cram_sampleqc/main.nf'
include { CHANNEL_ALIGN_CREATE_CSV } from '../../../subworkflows/local/channel_align_create_csv/main'

workflow FASTQ_ALIGN_PARABRICKS {

    take:
    ch_reads // channel: [mandatory] meta, reads
    ch_fasta // channel: [mandatory] meta, fasta
    ch_index // channel: [mandatory] meta, index
    ch_interval_file // channel: [optional] meta, intervals_bed_combined
    ch_known_sites // channel [optional] known_sites_indels
    val_output_fmt // either bam or cram
    ch_ngscheckmate_bed // channel: [mandatory] meta, ngscheckmate_bed
    ch_intervals_for_preprocessing // channel: [optional] meta, intervals_for_preprocessing

    main:

    ch_versions = Channel.empty()

    PARABRICKS_FQ2BAM(
        ch_reads,
        ch_fasta,
        ch_index,
        ch_interval_file,
        ch_known_sites,
        val_output_fmt
    )

    ch_versions = ch_versions.mix(PARABRICKS_FQ2BAM.versions)

    cram_variant_calling = PARABRICKS_FQ2BAM.out.cram

    CRAM_SAMPLEQC(cram_variant_calling,
            ch_ngscheckmate_bed,
            ch_fasta,
            params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator'),
            ch_intervals_for_preprocessing
    )

    ch_versions = ch_versions.mix(CRAM_SAMPLEQC.versions)

    CHANNEL_ALIGN_CREATE_CSV(
        PARABRICKS_FQ2BAM.out.cram.join(PARABRICKS_FQ2BAM.out.crai, failOnDuplicate: true, failOnMismatch: true), 
        params.outdir, 
        params.save_output_as_bam
    )

    emit:
    // TODO nf-core: edit emitted channels
    cram      = cram          // channel: [ val(meta), [ bam ] ]
    crai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

