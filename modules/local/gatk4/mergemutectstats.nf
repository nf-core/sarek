process GATK4_MERGEMUTECTSTATS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.4.1--hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.2.4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(stats)


    output:
    tuple val(meta), path("*.stats"), emit: stats
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = stats.collect{ " -stats ${it} "}.join()

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK MergeMutectStats] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" MergeMutectStats \\
        ${input} \\
        -O ${prefix}.stats \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
