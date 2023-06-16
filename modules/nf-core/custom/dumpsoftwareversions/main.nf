process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_single'

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    conda "bioconda::multiqc=1.14"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/multiqc:1.14--pyhdfd78af_0',
        singularity: 'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
