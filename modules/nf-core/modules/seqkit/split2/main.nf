// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::seqkit=0.16.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*${prefix}/*.gz"), emit: reads
    path "versions.yml"                     , emit: versions

    script:
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if(meta.single_end){
        """
        seqkit \\
            split2 \\
            $options.args \\
            --threads $task.cpus \\
            -1 $reads \\
            --out-dir $prefix

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        seqkit \\
            split2 \\
            $options.args \\
            --threads $task.cpus \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --out-dir $prefix

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
}
