process SNPSIFT_ANNMEM_CREATE_DB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/3570381a93c22449d48bdaa85097c5e8a075e90437565546acb2e40a29171bca/data'
        : 'community.wave.seqera.io/library/snpsift:5.3.0a--67d3871d6f67ac2b'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), val(fields)

    output:
    tuple val(meta), path("*.snpsift.vardb"), emit: database
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fields_arg = fields ? "-fields ${fields}" : ""

    """
    SnpSift \\
        annmem \\
        -create \\
        ${args} \\
        -dbfile ${vcf} \\
        ${fields_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${vcf}.snpsift.vardb
    touch ${vcf}.snpsift.vardb/chr1.snpsift.df

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """
}
