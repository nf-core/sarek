include { CNVKIT_ANTITARGET } from '../../../modules/nf-core/cnvkit/antitarget'
include { CNVKIT_REFERENCE  } from '../../../modules/nf-core/cnvkit/reference'

workflow PREPARE_REFERENCE_CNVKIT {
    take:
    fasta                  // channel: [mandatory] fasta
    intervals_bed_combined // channel: []

    main:
    // prepare a antitarget reference files for tumor_only mode of cnvkit
    CNVKIT_ANTITARGET(intervals_bed_combined.flatten().map { bed -> [[id: 'intervals'], bed] })
    CNVKIT_REFERENCE(fasta.map { _meta, fasta_ -> [fasta_] }, intervals_bed_combined, CNVKIT_ANTITARGET.out.bed.map { _meta, bed -> [bed] })

    emit:
    cnvkit_reference = CNVKIT_REFERENCE.out.cnn.collect()
}
