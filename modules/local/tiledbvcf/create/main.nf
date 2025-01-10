process TILEDBVCF_CREATE {
    tag "$meta.id"
    label 'process_low'

    conda "tiledb::tiledbvcf-py>=0.34.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/tiledb/tiledbvcf-py:latest' :
        'docker.io/tiledb/tiledbvcf-py:latest' }"

    input:
    tuple val(meta), val(db_name)

    output:
    tuple val(meta), val(db_name), emit: tiledb_db
    path "versions.yml", emit: versions
    
    when:
    params.tiledb_create_dataset

    script:
    def datasetExists = { db_name ->
        try {
            // Execute the command and capture the exit status
            def process = "tiledbvcf stat --uri ${db_name}".execute()
            process.waitFor() // Wait for the process to complete

            // Check the exit status
            if (process.exitValue() == 0) {
                return true // Dataset exists
            } else {
                return false // Dataset does not exist
            }
        } catch (Exception e) {
            return false // Handle any exceptions
        }
    }

    if (datasetExists(db_name)) {
        println "Dataset already exists. Skipping creation."
    } else {
        println "Dataset does not exist. Creating dataset."
        "tiledbvcf create --uri ${db_name}".execute()
    }

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiledbvcf: \$(tiledbvcf --version 2>&1 | sed 's/^.*version //; s/Using.*\$//')
    END_VERSIONS
    """
}