//
// Read QC and trimming
//

include { TRIMGALORE } from '../../modules/nf-core/modules/trimgalore/main'

workflow RUN_TRIMGALORE {
    take:
    reads         // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    TRIMGALORE(reads)

    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    emit:
    reads     = TRIMGALORE.out.reads // channel: [ val(meta), [ reads ] ]

    trim_html = TRIMGALORE.out.html // channel: [ val(meta), [ html ] ]
    trim_zip  = TRIMGALORE.out.zip  // channel: [ val(meta), [ zip ] ]
    trim_log  = TRIMGALORE.out.log  // channel: [ val(meta), [ txt ] ]

    versions = ch_versions          // channel: [ versions.yml ]
}
