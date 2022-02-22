process FGBIO_CALLMOLECULARCONSENSUSREADS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::fgbio=1.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:1.3.0--0' :
        'quay.io/biocontainers/fgbio:1.3.0--0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgbio \\
        CallMolecularConsensusReads \\
        -i $bam \\
        $args \\
        -o ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
