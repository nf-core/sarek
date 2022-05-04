process DEEPTOOLS_BAMCOVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::deeptools=3.5.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0':
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(meta), path(input), path(input_index)

    output:
    tuple val(meta), path("*.bigWig")   , emit: bigwig, optional: true
    tuple val(meta), path("*.bedgraph") , emit: bedgraph, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.bigWig"

    """
    bamCoverage \\
        --bam $input \\
        $args \\
        --numberOfProcessors ${task.cpus} \\
        --outFileName ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """
}
