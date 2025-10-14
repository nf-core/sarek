process SENTIEON_HAPLOTYPER {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f1dfe59ef66d7326b43db9ab1f39ce6220b358a311078c949a208f9c9815d4e/data'
        : 'community.wave.seqera.io/library/sentieon:202503.01--1863def31ed8e4d5'}"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals), path(recal_table)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dbsnp)
    tuple val(meta5), path(dbsnp_tbi)
    val emit_vcf
    val emit_gvcf

    output:
    // added the substring ".unfiltered" in the filename of the vcf-files since without that the g.vcf.gz-files were ending up in the vcf-channel
    tuple val(meta), path("*.unfiltered.vcf.gz"),     emit: vcf,      optional: true
    tuple val(meta), path("*.unfiltered.vcf.gz.tbi"), emit: vcf_tbi,  optional: true
    // these output-files have to have the extension ".vcf.gz", otherwise the subsequent GATK-MergeVCFs will fail.
    tuple val(meta), path("*.g.vcf.gz"),              emit: gvcf,     optional: true
    tuple val(meta), path("*.g.vcf.gz.tbi"),          emit: gvcf_tbi, optional: true
    path "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '' // options for the driver
    def args2 = task.ext.args2 ?: '' // options for the vcf generation
    def args3 = task.ext.args3 ?: '' // options for the gvcf generation
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = input instanceof List ? input.collect { "-i ${it}" }.join(' ') : "-i ${input}"
    def dbsnp_command = dbsnp ? "-d ${dbsnp} " : ""
    def interval_command = intervals ? "--interval ${intervals}" : ""
    def recal_table_command = recal_table ? "-q ${recal_table}" : ""
    def base_cmd = '--algo Haplotyper ' + dbsnp_command

    // The Sentieon --algo Haplotyper can create a VCF or gVCF but not both
    // Luckily, we can run it twice while reading the BAM once, therefore we construct the two separate commands
    // and run them twice while using the sentieon driver once. This allows us to create both types of VCF indels
    // one process

    // Create VCF command to export a VCF
    def vcf_cmd = emit_vcf
        ? base_cmd + args2 + ' --emit_mode ' + emit_vcf + ' ' + prefix + '.unfiltered.vcf.gz'
        : ""

    // Create a gVCF command to export a gVCF
    def gvcf_cmd = emit_gvcf
        ? base_cmd + args3 + ' --emit_mode gvcf ' + prefix + '.g.vcf.gz'
        : ""

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver \\
        ${args} \\
        -r ${fasta} \\
        -t ${task.cpus} \\
        ${interval_command} \\
        ${input_list} \\
        ${recal_table_command} \\
        ${vcf_cmd} \\
        ${gvcf_cmd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.unfiltered.vcf.gz
    touch ${prefix}.unfiltered.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sentieon: \$(echo \$(sentieon driver --version 2>&1) | sed -e "s/sentieon-genomics-//g")
    END_VERSIONS
    """
}
