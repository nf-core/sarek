/*
 * Read QC and trimming
 */

include { FASTQC     } from '../software/fastqc'
include { TRIMGALORE } from '../software/trimgalore'

workflow QC_TRIM {
    take:

    reads           // channel: [ val(meta), [ reads ] ]
    skip_fastqc     // boolean: true/false
    skip_trimming   // boolean: true/false
    fastqc_opts     //     map: options for FastQC module
    trimgalore_opts //     map: options for TrimGalore! module

    main:

    fastqc_html    = Channel.empty()
    fastqc_version = Channel.empty()
    fastqc_zip     = Channel.empty()
    if (!skip_fastqc) {
        FASTQC(reads, fastqc_opts)
        fastqc_html    = FASTQC.out.html
        fastqc_version = FASTQC.out.version
        fastqc_zip     = FASTQC.out.zip
    }

    trim_reads         = reads
    trimgalore_html    = Channel.empty()
    trimgalore_zip     = Channel.empty()
    trimgalore_log     = Channel.empty()
    trimgalore_version = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE(reads, trimgalore_opts)
        trim_reads         = TRIMGALORE.out.reads
        trimgalore_html    = TRIMGALORE.out.html
        trimgalore_zip     = TRIMGALORE.out.zip
        trimgalore_log     = TRIMGALORE.out.log
        trimgalore_version = TRIMGALORE.out.version
    }

    emit:

    fastqc_html        //    path: *.html
    fastqc_zip         //    path: *.zip
    fastqc_version     //    path: *.version.txt
    reads = trim_reads // channel: [ val(meta), [ reads ] ]
    trimgalore_html    //    path: *.html
    trimgalore_zip     //    path: *.zip
    trimgalore_log     //    path: *.txt
    trimgalore_version //    path: *.version.txt
}
