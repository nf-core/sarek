process FGBIO_GROUPREADSBYUMI {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::fgbio=2.0.2"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/fgbio:2.0.2--hdfd78af_0',
        singularity: 'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(taggedbam)
    val(strategy)

    output:
    tuple val(meta), path("*_umi-grouped.bam")  , emit: bam
    tuple val(meta), path("*_umi_histogram.txt"), emit: histogram
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    fgbio \\
        --tmp-dir=. \\
        GroupReadsByUmi \\
        -s $strategy \\
        $args \\
        -i $taggedbam \\
        -o ${prefix}_umi-grouped.bam \\
        -f ${prefix}_umi_histogram.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
