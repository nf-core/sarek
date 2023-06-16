process TABIX_TABIX {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::tabix=1.11"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/tabix:1.11--hdfd78af_0',
        singularity: 'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tbi"), optional:true, emit: tbi
    tuple val(meta), path("*.csi"), optional:true, emit: csi
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    tabix $args $tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${tab}.tbi
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
