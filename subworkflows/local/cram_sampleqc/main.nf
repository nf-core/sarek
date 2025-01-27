include { BAM_NGSCHECKMATE                           } from '../../../subworkflows/nf-core/bam_ngscheckmate/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL } from '../../../subworkflows/local/cram_qc_mosdepth_samtools/main'

workflow CRAM_SAMPLEQC {
    take:
    cram                        // channel: [ meta, cram, crai ]
    ngscheckmate_bed            // channel: [ meta, ngscheckmate_bed ]
    fasta                       // channel: [ meta, fasta ]
    skip_baserecalibration      // boolean
    intervals_for_preprocessing // channel: [ meta, intervals ]

    main:

    versions = Channel.empty()
    reports = Channel.empty()

    if (!skip_baserecalibration) {

        CRAM_QC_RECAL(
            cram,
            fasta,
            intervals_for_preprocessing,
        )

        // Gather QC reports
        reports = CRAM_QC_RECAL.out.reports.collect { _meta, report -> report }

        // Gather used softwares versions
        versions = versions.mix(CRAM_QC_RECAL.out.versions)
    }

    BAM_NGSCHECKMATE(cram.map { meta, cram_, _crai -> [meta, cram_] }, ngscheckmate_bed, fasta)
    versions = versions.mix(BAM_NGSCHECKMATE.out.versions.first())

    emit:
    corr_matrix = BAM_NGSCHECKMATE.out.corr_matrix // channel: [ meta, corr_matrix ]
    matched     = BAM_NGSCHECKMATE.out.matched // channel: [ meta, matched ]
    all         = BAM_NGSCHECKMATE.out.all // channel: [ meta, all ]
    vcf         = BAM_NGSCHECKMATE.out.vcf // channel: [ meta, vcf ]
    pdf         = BAM_NGSCHECKMATE.out.pdf // channel: [ meta, pdf ]
    reports
    versions    // channel: [ versions.yml ]
}
