process MINIMAP2_INDEX {
    tag "$fasta"
    label 'process_high'

    conda     (params.enable_conda ? "bioconda::minimap2=2.17" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2:2.17--hed695b0_3' :
        'quay.io/biocontainers/minimap2:2.17--hed695b0_3' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi")  , emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def preset    = "-ax map-ont"
    def kmer      = ""
    def stranded  = ""
    def junctions = ""
    """
    minimap2 \\
        $preset \\
        $kmer \\
        $stranded \\
        $junctions \\
        -t $task.cpus \\
        -d ${fasta}.mmi \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
