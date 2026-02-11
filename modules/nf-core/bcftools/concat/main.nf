process BCFTOOLS_CONCAT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    tuple val(meta), path(vcfs), path(tbi)

    output:
    tuple val(meta), path("${prefix}.${extension}"), emit: vcf
    tuple val(meta), path("${prefix}.${extension}.tbi"), emit: tbi, optional: true
    tuple val(meta), path("${prefix}.${extension}.csi"), emit: csi, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def tbi_names = tbi.findAll { file -> !(file instanceof List) }.collect { file -> file.name }
    def create_input_index = vcfs.collect { vcf -> tbi_names.contains(vcf.name + ".tbi") || tbi_names.contains(vcf.name + ".csi") ? "" : "tabix ${vcf}" }.join("\n    ")
    extension = args.contains("--output-type b") || args.contains("-Ob")
        ? "bcf.gz"
        : args.contains("--output-type u") || args.contains("-Ou")
            ? "bcf"
            : args.contains("--output-type z") || args.contains("-Oz")
                ? "vcf.gz"
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def input = vcfs.sort { vcf -> vcf.toString() }.join(" ")
    """
    ${create_input_index}

    bcftools concat \\
        --output ${prefix}.${extension} \\
        ${args} \\
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
                : args.contains("--output-type v") || args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def index_extension = args.contains("--write-index=tbi") || args.contains("-W=tbi")
        ? "tbi"
        : args.contains("--write-index=csi") || args.contains("-W=csi")
            ? "csi"
            : args.contains("--write-index") || args.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index_extension.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index_extension}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
