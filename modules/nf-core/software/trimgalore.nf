include { initOptions; saveFiles; getSoftwareName } from './functions'

process TRIMGALORE {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (options.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/trim-galore:0.6.5--0"

    conda (params.conda ? "trim-galore=0.6.5" : null)

    input:
      tuple val(meta), path(reads)
      val options

    output:
      tuple val(meta), path("*_1.fq.gz"), path("*_2.fq.gz"), emit: reads
      path "*.html" ,                                        emit: html optional true
      path "*.txt" ,                                         emit: log
      path "*.version.txt",                                  emit: version
      path "*.zip" ,                                         emit: zip optional true

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Clipping presets have to be evaluated in the context of SE/PE
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''

    // Added soft-links to original fastqs for consistent naming in MultiQC
    prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

    trim_galore \\
        ${options.args} \\
        --cores ${cores} \\
        --paired \\
        --gzip \\
        ${c_r1} \\
        ${c_r2} \\
        ${tpc_r1} \\
        ${tpc_r2} \\
        ${prefix}_1.fastq.gz \\
        ${prefix}_2.fastq.gz

    trim_galore --version > trim_galore.version.txt
    """
}
