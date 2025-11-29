process VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/834c886cf862aade1d5d28a4a750b0676dfbced7300b10a795d5ebb993fa2586/data'
        : 'community.wave.seqera.io/library/varlociraptor:8.9.0--7199b2ed2f0e184f'}"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.alignment-properties.json"), emit: alignment_properties_json
    path "versions.yml",                                  emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment-properties.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """
}
