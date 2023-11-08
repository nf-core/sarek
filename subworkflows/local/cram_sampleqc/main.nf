include { BAM_NGSCHECKMATE } from '../../../subworkflows/nf-core/bam_ngscheckmate/main'

workflow CRAM_SAMPLEQC {

    take:
    ch_cram             // channel: [ val(meta), cram, crai ]
    ngscheckmate_bed    // channel: [ ngscheckmate_bed ]
    fasta               // channel: [ fasta ]

    main:

    ch_versions = Channel.empty()

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

    versions = ch_versions                     // channel: [ versions.yml ]
}

