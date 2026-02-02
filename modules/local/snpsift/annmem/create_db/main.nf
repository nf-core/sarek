process SNPSIFT_ANNMEM_CREATE_DB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/56f40c9dc9094c1e74d82a7120d64260073b123976503d8943ec19f5c0627a3a/data' :
        'community.wave.seqera.io/library/snpsift:5.4.0a--6680e6faf23ef488' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), val(fields)

    output:
    tuple val(meta), path("*.snpsift.vardb"), emit: database
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift split -h 2>&1 | head -1 | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g'"), topic: 'versions'

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Convert semicolons to commas for SnpSift (semicolons used in CSV to avoid delimiter conflict)
    def fields_csv = fields ? fields.replace(';', ',') : ''
    def fields_arg = fields_csv ? "-fields ${fields_csv}" : ""

    """
    SnpSift \\
        annmem \\
        -create \\
        ${args} \\
        -dbfile ${vcf} \\
        ${fields_arg}
    """

    stub:
    """
    mkdir -p ${vcf}.snpsift.vardb
    touch ${vcf}.snpsift.vardb/chr1.snpsift.df
    """
}
