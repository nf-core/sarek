process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_low'

    conda     (params.enable_conda ? "bioconda::minimap2=2.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.17--hed695b0_3' :
        'quay.io/biocontainers/minimap2:2.17--hed695b0_3' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)


    output:

    tuple val(meta), path("*.sam"), emit: align_sam
    path  "versions.yml"          , emit: versions

    script:
    def preset    = "-ax map-ont"
    def kmer      = ""
    def stranded  = ""
    def junctions = ""
    def md        = "--MD"
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    minimap2 \\
        $preset \\
        $kmer \\
        $stranded \\
        $junctions \\
        $md \\
        -t $task.cpus \\
        $index \\
        $reads > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
