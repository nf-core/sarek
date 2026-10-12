process SNPSIFT_ANNMEMCREATE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c5/c5ee8f24fb1bcd6e98ecb91fcf013f35a8db193dc9af4883c468e623aaebd7cf/data'
        : 'community.wave.seqera.io/library/htslib_snpsift:6167cd8f20036c55'}"

    input:
    tuple val(meta), path(db_vcf), path(db_vcf_tbi), val(db_fields)

    output:
    tuple val(meta), path("*.snpsift.vardb"), emit: database
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift -version 2>&1 | grep -oE '[0-9]+\\.[0-9]+[a-z]?'"), emit: versions_snpsift, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def fields = db_fields instanceof List ? db_fields.join(',') : db_fields

    """
    SnpSift \\
        annmem \\
        -create \\
        ${args} \\
        -dbfile ${db_vcf} \\
        ${fields ? "-fields ${fields}" : ""}
    """

    stub:
    """
    mkdir -p ${db_vcf}.snpsift.vardb
    touch ${db_vcf}.snpsift.vardb/chr1.snpsift.df
    """
}
