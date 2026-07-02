include { BAM_NGSCHECKMATE                           } from '../../../subworkflows/nf-core/bam_ngscheckmate'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL } from '../../../subworkflows/local/cram_qc_mosdepth_samtools'

workflow CRAM_SAMPLEQC {
    take:
    cram                        // channel: [ val(meta), cram, crai ]
    ngscheckmate_bed            // channel: [ ngscheckmate_bed ]
    fasta_fai                   // channel: [ fasta, fasta_fai ]
    skip_baserecalibration      // boolean:
    intervals_for_preprocessing // channel:

    main:
    versions = Channel.empty()
    reports = Channel.empty()
    // Split fasta channel into [fasta] for QC and [meta, fasta, fai] for BAM_NGSCHECKMATE
    fasta_ch  = fasta_fai.map { f, f_ -> [f] }
    fasta_fai_ch = fasta_fai.map { _f, f_ -> [[id: 'fasta_fai'], f_] }

    if (!skip_baserecalibration) {

        CRAM_QC_RECAL(
            cram,
            fasta_ch,
            intervals_for_preprocessing,
        )

        // Gather QC reports
        reports = CRAM_QC_RECAL.out.reports.collect { _meta, report -> report }

        // Gather used softwares versions
        versions = versions.mix(CRAM_QC_RECAL.out.versions)
    }

    BAM_NGSCHECKMATE(cram.map { meta, cram_, _crai -> [meta, cram_] }, ngscheckmate_bed.map { bed -> [[id: "ngscheckmate"], bed] }, fasta_fai_ch)

    emit:
    corr_matrix = BAM_NGSCHECKMATE.out.corr_matrix // channel: [ meta, corr_matrix ]
    matched     = BAM_NGSCHECKMATE.out.matched // channel: [ meta, matched ]
    all         = BAM_NGSCHECKMATE.out.all // channel: [ meta, all ]
    vcf         = BAM_NGSCHECKMATE.out.vcf // channel: [ meta, vcf ]
    pdf         = BAM_NGSCHECKMATE.out.pdf // channel: [ meta, pdf ]
    reports
    versions    // channel: [ versions.yml ]
}
