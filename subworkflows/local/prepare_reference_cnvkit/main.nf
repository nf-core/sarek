include { CNVKIT_ANTITARGET } from '../../../modules/nf-core/cnvkit/antitarget'
include { CNVKIT_REFERENCE  } from '../../../modules/nf-core/cnvkit/reference'

workflow PREPARE_REFERENCE_CNVKIT {
    take:
    fasta                  // channel: [mandatory] fasta
    intervals_bed_combined // channel: []

    main:
    versions = Channel.empty()

    // prepare a antitarget reference files for tumor_only mode of cnvkit
    CNVKIT_ANTITARGET(intervals_bed_combined.flatten().map { it -> [[id: 'intervals'], it] })
    CNVKIT_REFERENCE(fasta.map { _meta, fasta_ -> [fasta_] }, intervals_bed_combined, CNVKIT_ANTITARGET.out.bed.map { _meta, bed -> [bed] })

    versions = versions.mix(CNVKIT_ANTITARGET.out.versions)
    versions = versions.mix(CNVKIT_REFERENCE.out.versions)

    emit:
    cnvkit_reference = CNVKIT_REFERENCE.out.cnn.collect()
    versions
}
