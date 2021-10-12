//
// Read QC and trimming
//

params.fastqc_options     = [:]
params.trimgalore_options = [:]

include { FASTQC }     from '../../modules/nf-core/modules/fastqc/main'     addParams( options: params.fastqc_options )
include { TRIMGALORE } from '../../modules/nf-core/modules/trimgalore/main' addParams( options: params.trimgalore_options )

workflow FASTQC_TRIMGALORE {
    take:
    reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc   // boolean: true/false
    skip_trimming // boolean: true/false

    main:
    fastqc_html    = Channel.empty()
    fastqc_zip     = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc) {
        FASTQC(reads)
        fastqc_html    = FASTQC.out.html
        fastqc_zip     = FASTQC.out.zip
        fastqc_version = FASTQC.out.versions.first()
    }

    trim_reads = reads
    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()
    trimgalore_version = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE(reads)
        trim_reads = TRIMGALORE.out.reads
        trim_html  = TRIMGALORE.out.html
        trim_zip   = TRIMGALORE.out.zip
        trim_log   = TRIMGALORE.out.log
        trimgalore_version = TRIMGALORE.out.versions.first()
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]

    fastqc_html        // channel: [ val(meta), [ html ] ]
    fastqc_zip         // channel: [ val(meta), [ zip ] ]
    fastqc_version     //    path: versions.yml

    trim_html          // channel: [ val(meta), [ html ] ]
    trim_zip           // channel: [ val(meta), [ zip ] ]
    trim_log           // channel: [ val(meta), [ txt ] ]
    trimgalore_version //    path: versions.yml
}
