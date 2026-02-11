process BCFTOOLS_ISEC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(vcfs), path(tbis), path(file_list), path(targets_file), path(regions_file)

    output:
    tuple val(meta), path("${prefix}", type: "dir"), emit: results
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    targets_file_args = targets_file ? "-T ${targets_file}" : ''
    regions_file_args = regions_file ? "-R ${regions_file}" : ''
    vcf_files = file_list ? "-l ${file_list}" : "${vcfs}"

    """
    bcftools isec  \\
        ${args} \\
        ${targets_file_args} \\
        ${regions_file_args} \\
        -p ${prefix} \\
        ${vcf_files} \\
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/README.txt
    touch ${prefix}/sites.txt
    echo "" | gzip > ${prefix}/0000.vcf.gz
    touch ${prefix}/0000.vcf.gz.tbi
    echo "" | gzip > ${prefix}/0001.vcf.gz
    touch ${prefix}/0001.vcf.gz.tbi
    """
}
