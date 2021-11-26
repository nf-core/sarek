process SEQKIT_SPLIT2 {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqkit=0.16.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*${prefix}/*.gz"), emit: reads
    path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    if(meta.single_end){
        """
        seqkit \\
            split2 \\
            $args \\
            --threads $task.cpus \\
            -1 $reads \\
            --out-dir $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        seqkit \\
            split2 \\
            $args \\
            --threads $task.cpus \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            --out-dir $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$(echo \$(seqkit 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
}
