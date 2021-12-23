//
// ANNOTATION
//

include { BAM2FQ                                  } from '../nf-core/annotation_snpeff/main'

workflow ANNOTATE {
    take:
    input  // bam or cram meta, bam/cram
    fasta  // [optional] for cram files, if no header is provided in cram files

    main:
    ch_versions = Channel.empty()



    emit:
    reads
    versions = ch_versions     //    path: versions.yml
}
