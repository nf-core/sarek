process VARLOCIRAPTOR_PREPROCESS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::varlociraptor=8.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varlociraptor:8.1.1--hc349b7f_0':
        'biocontainers/varlociraptor:8.1.1--hc349b7f_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(candidates), path(alignment_json)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bcf.gz"), emit: bcf_gz, optional: true
    tuple val(meta), path("*.vcf.gz"), emit: vcf_gz, optional: true
    tuple val(meta), path("*.bcf")   , emit: bcf   , optional: true
    tuple val(meta), path("*.vcf")   , emit: vcf   , optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.vcf.gz"
    def alignment_properties_json = alignment_json ? "--alignment-properties ${alignment_json}" : ""
    """
    varlociraptor preprocess variants \\
        $fasta \\
        $alignment_properties_json \\
        --bam $bam \\
        --candidates $candidates \\
        ${args} \\
        > ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varlociraptor: \$(echo \$(varlociraptor --version 2>&1) | sed 's/^.*varlociraptor //; s/:.*\$//' )
    END_VERSIONS
    """
}
