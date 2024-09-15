process TILEDBVCF_CREATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::tiledbvcf-py=0.35.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'tiledb/tiledbvcf-py:latest' :
        'tiledb/tiledbvcf-py:latest' }"

    input:
    tuple val(meta), val(db_name)

    output:
    tuple val(meta), path("${params.tiledb_dataset_name}"), emit: tiledb_db
    path "versions.yml"           , emit: versions

    when:
    params.tiledb_create_dataset

    script:
    """
    tiledbvcf create --uri ${params.tiledb_dataset_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiledbvcf: \$(tiledbvcf --version 2>&1 | sed 's/^.*version //; s/Using.*\$//')
    END_VERSIONS
    """
}