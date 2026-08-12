process VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/varlociraptor%3A8.9.5--h24073b4_0'
        : 'quay.io/biocontainers/varlociraptor:8.9.5--h24073b4_0'}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.alignment-properties.json"), emit: alignment_properties_json
    tuple val("${task.process}"), val('varlociraptor'), eval("varlociraptor --version | sed 's/^varlociraptor //'"), topic: versions, emit: versions_varlociraptor

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    varlociraptor estimate alignment-properties \\
        ${fasta} \\
        --bams ${bam} \\
        ${args} \\
        > ${prefix}.alignment-properties.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment-properties.json
    """
}
