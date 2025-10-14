process VCFLIB_VCFFILTER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fc/fc33d59c090cef123aca26ae17fbddbd596640304d8325cbd5816229fa2c05ee/data'
        : 'community.wave.seqera.io/library/vcflib:1.0.14--cc8ffb2c1a080797'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.filter"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '1.0.14'
    if (!(args.contains("-f") || args.contains("--info-filter") || args.contains("-g") || args.contains("--genotype-filter"))) {
        error("VCFLIB_VCFFILTER requires either the -f/--info-filter or -g/--genotype-filter arguments to be supplied using ext.args.")
    }
    if ("${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    vcffilter \\
        ${args} \\
        ${vcf} \\
        | bgzip -c ${args2} > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: ${VERSION}
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.filter"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '1.0.14'
    if ("${vcf}" == "${prefix}.vcf.gz") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }
    """
    echo | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcflib: ${VERSION}
    END_VERSIONS
    """
}
