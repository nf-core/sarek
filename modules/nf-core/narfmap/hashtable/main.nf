process NARFMAP_HASHTABLE {
    tag "$fasta"
    label 'process_high'

    // TODO Add a singularity image
    // TODO Make a seperate image
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nf-core/modules/narfmap_align:narfmap--8a04bcf8bd9b6242':
        'nf-core/modules/narfmap_align:narfmap--8a04bcf8bd9b6242' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("narfmap")    , emit: hashmap
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir narfmap
    dragen-os \\
        --build-hash-table true \\
        --ht-reference $fasta \\
        --output-directory narfmap \\
        $args \\
        --ht-num-threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        narfmap: \$(echo \$(dragen-os --version 2>&1))
    END_VERSIONS
    """

    stub:
    """
    mkdir narfmap

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        narfmap: \$(echo \$(dragen-os --version 2>&1))
    END_VERSIONS
    """

}
