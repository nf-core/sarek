//
// Read QC and trimming
//

include { FASTQC     } from '../../modules/nf-core/modules/fastqc/main'
include { TRIMGALORE } from '../../modules/nf-core/modules/trimgalore/main'

workflow FASTQC_TRIMGALORE {
    take:
    reads         // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = Channel.empty()

    FASTQC(reads)
    TRIMGALORE(reads)

    out_reads = reads.mix(TRIMGALORE.out.reads)

    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    emit:
    reads = out_reads               // channel: [ val(meta), [ reads ] ]

    fastqc_html = FASTQC.out.html   // channel: [ val(meta), [ html ] ]
    fastqc_zip  = FASTQC.out.zip    // channel: [ val(meta), [ zip ] ]

    trim_html = TRIMGALORE.out.html // channel: [ val(meta), [ html ] ]
    trim_zip  = TRIMGALORE.out.zip  // channel: [ val(meta), [ zip ] ]
    trim_log  = TRIMGALORE.out.log  // channel: [ val(meta), [ txt ] ]

    versions = ch_versions          // channel: [ versions.yml ]
}
