include { BAM_NGSCHECKMATE } from '../../../subworkflows/nf-core/bam_ngscheckmate/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL  } from '../../../subworkflows/local/cram_qc_mosdepth_samtools/main'

workflow CRAM_SAMPLEQC {

    take:
    ch_cram                     // channel: [ val(meta), cram, crai ]
    ngscheckmate_bed            // channel: [ ngscheckmate_bed ]
    fasta                       // channel: [ fasta ]
    skip_baserecalibration      // boolean:
    intervals_for_preprocessing // channel:

    main:

    ch_versions = Channel.empty()
    reports     = Channel.empty()

    if(!skip_baserecalibration){

        CRAM_QC_RECAL(
            ch_cram,
            fasta,
            intervals_for_preprocessing)

        // Gather QC reports
        reports = CRAM_QC_RECAL.out.reports.collect{ meta, report -> report }

        // Gather used softwares versions
        ch_versions = ch_versions.mix(CRAM_QC_RECAL.out.versions)
    }

    ch_ngscheckmate_bed = ngscheckmate_bed.map{bed -> [[id: "ngscheckmate"], bed]}

    ch_fasta            = fasta.map{fasta -> [[id: "genome"], fasta]}

    BAM_NGSCHECKMATE ( ch_cram.map{meta, cram, crai -> [meta, cram]}, ch_ngscheckmate_bed, ch_fasta)
    ch_versions = ch_versions.mix(BAM_NGSCHECKMATE.out.versions.first())

    emit:
    corr_matrix  = BAM_NGSCHECKMATE.out.corr_matrix  // channel: [ meta, corr_matrix ]
    matched      = BAM_NGSCHECKMATE.out.matched      // channel: [ meta, matched ]
    all          = BAM_NGSCHECKMATE.out.all          // channel: [ meta, all ]
    vcf          = BAM_NGSCHECKMATE.out.vcf          // channel: [ meta, vcf ]
    pdf          = BAM_NGSCHECKMATE.out.pdf          // channel: [ meta, pdf ]
    reports

    versions = ch_versions                     // channel: [ versions.yml ]
}

