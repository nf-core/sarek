include { BAM_NGSCHECKMATE                           } from '../../../subworkflows/nf-core/bam_ngscheckmate'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL } from '../../../subworkflows/local/cram_qc_mosdepth_samtools'

workflow CRAM_SAMPLEQC {
    take:
    cram                        // channel: [ val(meta), cram, crai ]
    ngscheckmate_bed            // channel: [ ngscheckmate_bed ]
    fasta_fai                   // channel: [ val(meta), fasta, fasta_fai ]
    skip_baserecalibration      // boolean:
    intervals_for_preprocessing // channel:
    skip_tools                  // string:  --skip_tools
    tools                       // string:  --tools

    main:
    reports = channel.empty()

    if (!skip_baserecalibration) {

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

    corr_matrix = channel.empty()
    matched     = channel.empty()
    all         = channel.empty()
    vcf         = channel.empty()
    pdf         = channel.empty()

    if (tools && tools.split(',').contains('ngscheckmate')) {
        BAM_NGSCHECKMATE(cram.map { meta, cram_, _crai -> [meta, cram_] }, ngscheckmate_bed.map { bed -> [[id: "ngscheckmate"], bed] }, fasta_fai)
        corr_matrix = BAM_NGSCHECKMATE.out.corr_matrix
        matched     = BAM_NGSCHECKMATE.out.matched
        all         = BAM_NGSCHECKMATE.out.all
        vcf         = BAM_NGSCHECKMATE.out.vcf
        pdf         = BAM_NGSCHECKMATE.out.pdf
    }

    emit:
    corr_matrix // channel: [ meta, corr_matrix ]
    matched     // channel: [ meta, matched ]
    all         // channel: [ meta, all ]
    vcf         // channel: [ meta, vcf ]
    pdf         // channel: [ meta, pdf ]
    reports
}
