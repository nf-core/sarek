process CONTROLFREEC_FREEC2BED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::control-freec=11.6"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/control-freec:11.6--h1b792b2_1',
        singularity: 'https://depot.galaxyproject.org/singularity/control-freec:11.6--h1b792b2_1',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    tuple val(meta), path(ratio)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    freec2bed.pl -f ${ratio} ${args} > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: \$(echo \$(freec -version 2>&1) | sed 's/^.*Control-FREEC  //; s/:.*\$//' | sed -e "s/Control-FREEC v//g" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        controlfreec: \$(echo \$(freec -version 2>&1) | sed 's/^.*Control-FREEC  //; s/:.*\$//' | sed -e "s/Control-FREEC v//g" )
    END_VERSIONS
    """
}
