process SAMBAMBA_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity//sambamba:1.0--h98b6b92_0':
        'biocontainers/sambamba:1.0--h98b6b92_0' }"

    // [[patient:test, sample:test, sex:XX, status:0, id:test, data_type:bam], /workspace/nextflow_dir/work/5f/bd451fdce7965b839d600fcc122530/test.bam]
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sambamba \\
        markdup \\
        $args \\
        -t $task.cpus \\
        $bam \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sambamba: \$(echo \$(sambamba --version 2>&1) | awk '{print \$2}' )
    END_VERSIONS
    """
}
