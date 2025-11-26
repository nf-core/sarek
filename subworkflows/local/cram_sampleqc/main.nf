include { BAM_NGSCHECKMATE                           } from '../../../subworkflows/nf-core/bam_ngscheckmate'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL } from '../../../subworkflows/local/cram_qc_mosdepth_samtools'

workflow CRAM_SAMPLEQC {
    take:
    cram                        // channel: [ val(meta), cram, crai ]
    ngscheckmate_bed            // channel: [ ngscheckmate_bed ]
    fasta                       // channel: [ fasta ]
    skip_baserecalibration      // boolean:
    intervals_for_preprocessing // channel:

    main:

    versions = Channel.empty()
    reports = Channel.empty()
    samtools_stats = Channel.empty()
    mosdepth_global = Channel.empty()
    mosdepth_region = Channel.empty()
    mosdepth_summary = Channel.empty()
    mosdepth_regions_bed = Channel.empty()
    mosdepth_regions_csi = Channel.empty()

    if (!skip_baserecalibration) {

        CRAM_QC_RECAL(
            cram,
            fasta,
            intervals_for_preprocessing,
        )

        // Gather QC reports
        reports = CRAM_QC_RECAL.out.reports.collect { _meta, report -> report }

        // Capture individual QC channels with meta
        samtools_stats = CRAM_QC_RECAL.out.samtools_stats
        mosdepth_global = CRAM_QC_RECAL.out.mosdepth_global
        mosdepth_region = CRAM_QC_RECAL.out.mosdepth_region
        mosdepth_summary = CRAM_QC_RECAL.out.mosdepth_summary
        mosdepth_regions_bed = CRAM_QC_RECAL.out.mosdepth_regions_bed
        mosdepth_regions_csi = CRAM_QC_RECAL.out.mosdepth_regions_csi

        // Gather used softwares versions
        versions = versions.mix(CRAM_QC_RECAL.out.versions)
    }

    BAM_NGSCHECKMATE(cram.map { meta, cram_, _crai -> [meta, cram_] }, ngscheckmate_bed.map { bed -> [[id: "ngscheckmate"], bed] }, fasta)
    versions = versions.mix(BAM_NGSCHECKMATE.out.versions)

    emit:
    corr_matrix = BAM_NGSCHECKMATE.out.corr_matrix // channel: [ meta, corr_matrix ]
    matched     = BAM_NGSCHECKMATE.out.matched // channel: [ meta, matched ]
    all         = BAM_NGSCHECKMATE.out.all // channel: [ meta, all ]
    vcf         = BAM_NGSCHECKMATE.out.vcf // channel: [ meta, vcf ]
    pdf         = BAM_NGSCHECKMATE.out.pdf // channel: [ meta, pdf ]
    reports
    samtools_stats                                 // channel: [ meta, stats ]
    mosdepth_global                                // channel: [ meta, txt ]
    mosdepth_region                                // channel: [ meta, txt ]
    mosdepth_summary                               // channel: [ meta, txt ]
    mosdepth_regions_bed                           // channel: [ meta, bed.gz ]
    mosdepth_regions_csi                           // channel: [ meta, csi ]
    versions    // channel: [ versions.yml ]
}
