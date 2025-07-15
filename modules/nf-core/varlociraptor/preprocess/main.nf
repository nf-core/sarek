process VARLOCIRAPTOR_PREPROCESS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varlociraptor:8.7.3--ha8ac579_2':
        'biocontainers/varlociraptor:8.7.3--ha8ac579_2' }"

    input:
    tuple val(meta), path(bam), path(bai), path(candidates), path(alignment_json)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bcf"), emit: bcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def alignment_properties_json = alignment_json ? "--alignment-properties ${alignment_json}" : ""
    """
    varlociraptor preprocess variants \\
        ${fasta} \\
        ${alignment_properties_json} \\
        --bam ${bam} \\
        --candidates ${candidates} \\
        ${args} \\
        --output ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """
}
