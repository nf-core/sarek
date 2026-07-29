process BCFTOOLS_ANNOTATE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(input), path(index), path(annotations), path(annotations_index), path(columns), path(header_lines), path(rename_chrs)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: vcf
    tuple val(meta), path("${prefix}.${extension}.{tbi,csi}"), emit: index, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def annotations_file = annotations ? "--annotations ${annotations}" : ''
    def columns_file = columns ? "--columns-file ${columns}" : ''
    def header_file = header_lines ? "--header-lines ${header_lines}" : ''
    def rename_chrs_file = rename_chrs ? "--rename-chrs ${rename_chrs}" : ''
    extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov") ? "vcf" : "vcf"
    def index_command = !index ? "bcftools index ${input}" : ''

    if ("${input}" == "${prefix}.${extension}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    ${index_command}

    bcftools \\
        annotate \\
        ${args} \\
        ${annotations_file} \\
        ${columns_file} \\
        ${header_file} \\
        ${rename_chrs_file} \\
        --output ${prefix}.${extension} \\
        --threads ${task.cpus} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov") ? "vcf" : "vcf"
    def index_extension = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W") ? "csi" : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index_extension.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index_extension}" : ""

    if ("${input}" == "${prefix}.${extension}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
