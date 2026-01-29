process SNPSIFT_ANNMEM_CREATE_DB {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/3570381a93c22449d48bdaa85097c5e8a075e90437565546acb2e40a29171bca/data'
        : 'community.wave.seqera.io/library/htslib_snpsift:4df051c7ff79f7f9'}"

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
