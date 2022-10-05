//
// Read QC
//

include { FASTQC     } from '../../modules/nf-core/fastqc/main'

workflow RUN_FASTQC {
    take:
    reads         // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    FASTQC(reads)

    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    emit:

    fastqc_html = FASTQC.out.html   // channel: [ val(meta), [ html ] ]
    fastqc_zip  = FASTQC.out.zip    // channel: [ val(meta), [ zip ] ]

    versions = ch_versions          // channel: [ versions.yml ]
}
