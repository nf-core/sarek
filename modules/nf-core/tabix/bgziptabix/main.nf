process TABIX_BGZIPTABIX {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz"), path("*.tbi"), optional: true, emit: gz_tbi
    tuple val(meta), path("*.gz"), path("*.csi"), optional: true, emit: gz_csi
    path  "versions.yml" ,                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bgzip  --threads ${task.cpus} -c $args $input > ${prefix}.${input.getExtension()}.gz
    tabix $args2 ${prefix}.${input.getExtension()}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${input.getExtension()}.gz
    touch ${prefix}.${input.getExtension()}.gz.tbi
    touch ${prefix}.${input.getExtension()}.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
