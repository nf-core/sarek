process DRAGMAP_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:7eed251370ac7f3537c3d9472cdb2f9f5d8da1c5-0':
        'biocontainers/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:7eed251370ac7f3537c3d9472cdb2f9f5d8da1c5-0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(hashmap)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path('*.log'), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-1 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def samtools_command = sort_bam ? 'sort' : 'view'

    """
    dragen-os \\
        -r $hashmap \\
        $args \\
        --num-threads $task.cpus \\
        $reads_command \\
        2> >(tee ${prefix}.dragmap.log >&2) \\
        | samtools $samtools_command $args2 --threads $task.cpus -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
