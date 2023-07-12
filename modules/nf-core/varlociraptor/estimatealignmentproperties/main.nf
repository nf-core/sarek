process VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES {
    tag "$meta.id"
    label 'process_single'
    conda "bioconda::varlociraptor=8.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varlociraptor:8.1.1--hc349b7f_0':
        'biocontainers/varlociraptor:8.1.1--hc349b7f_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.alignment-properties.json"), emit: alignment_properties_json
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    varlociraptor estimate alignment-properties \\
        $fasta \\
        --bam $bam \\
        $args \\
        > ${prefix}.alignment-properties.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """
}
