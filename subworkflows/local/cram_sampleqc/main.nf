include { BAM_NGSCHECKMATE                           } from '../../../subworkflows/nf-core/bam_ngscheckmate'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL } from '../../../subworkflows/local/cram_qc_mosdepth_samtools'

workflow CRAM_SAMPLEQC {
    take:
    cram                        // channel: [ val(meta), cram, crai ]
    ngscheckmate_bed            // channel: [ ngscheckmate_bed ]
    fasta_fai                   // channel: [ val(meta), fasta, fasta_fai ]
    intervals_for_preprocessing // channel:
    skip_tools                  // list:    lowercased --skip_tools
    tools                       // list:    lowercased --tools

    main:
    reports = channel.empty()

    if ('baserecalibrator' !in skip_tools) {

        CRAM_QC_RECAL(
            cram,
            fasta_fai.map{meta, fasta, _fai -> [meta, fasta]},
            fasta_fai.map{meta, _fasta, fai -> [meta, fai]},
            intervals_for_preprocessing,
            skip_tools,
        )

        // Gather QC reports
        reports = CRAM_QC_RECAL.out.reports.collect { _meta, report -> report }
    }

    BAM_NGSCHECKMATE(cram.filter { 'ngscheckmate' in tools }.map { meta, cram_, _crai -> [meta, cram_] }, ngscheckmate_bed.map { bed -> [[id: "ngscheckmate"], bed] }, fasta_fai)

    emit:
    corr_matrix = BAM_NGSCHECKMATE.out.corr_matrix // channel: [ meta, corr_matrix ]
    matched     = BAM_NGSCHECKMATE.out.matched // channel: [ meta, matched ]
    all         = BAM_NGSCHECKMATE.out.all // channel: [ meta, all ]
    vcf         = BAM_NGSCHECKMATE.out.vcf // channel: [ meta, vcf ]
    pdf         = BAM_NGSCHECKMATE.out.pdf // channel: [ meta, pdf ]
    reports
}
