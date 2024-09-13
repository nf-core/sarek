process TILEDBVCF_CREATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::tiledbvcf-py=0.20.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiledbvcf-py:0.20.0--py39h9197a36_0' :
        'biocontainers/tiledbvcf-py:0.20.0--py39h9197a36_0' }"

    input:
    tuple val(meta), val(db_name)

    output:
    tuple val(meta), path("${params.tiledb_store_name}"), emit: tiledb_db
    path "versions.yml"           , emit: versions

    when:
    params.tiledb_create_store

    script:
    """
    tiledbvcf create --uri ${params.tiledb_store_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiledbvcf: \$(tiledbvcf --version 2>&1 | sed 's/^.*version //; s/Using.*\$//')
    END_VERSIONS
    """
}