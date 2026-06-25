process SENTIEON_DEDUP {
    tag "${meta.id}"
    label 'process_medium'
    label 'sentieon'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.cram"),                emit: cram, optional: true
    tuple val(meta), path("*.crai"),                emit: crai, optional: true
    tuple val(meta), path("*.bam"),                 emit: bam,  optional: true
    tuple val(meta), path("*.bai"),                 emit: bai
    tuple val(meta), path("*.score"),               emit: score
    tuple val(meta), path("*.metrics"),             emit: metrics
    tuple val(meta), path("*.metrics.multiqc.tsv"), emit: metrics_multiqc_tsv
    tuple val("${task.process}"), val('sentieon'), eval('sentieon driver --version | sed "s/.*-//g"'), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.cram"
    def metrics = "${prefix}.metrics"
    def input_list = bam.collect {input -> "-i ${input}" }.join(' ')
    def prefix_basename = prefix.substring(0, prefix.lastIndexOf("."))
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon driver ${args} -t ${task.cpus} ${input_list} -r ${fasta} --algo LocusCollector ${args2} --fun score_info ${prefix_basename}.score
    sentieon driver ${args3} -t ${task.cpus} ${input_list} -r ${fasta} --algo Dedup ${args4} --score_info ${prefix_basename}.score --metrics ${metrics} ${prefix}

    # This following tsv-file is produced in order to get a proper tsv-file with Dedup-metrics for importing in MultiQC as "custom content".
    # It should be removed once MultiQC has a module for displaying Dedup-metrics.
    head -3 ${metrics} > ${metrics}.multiqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.cram"
    def metrics = "${prefix}.metrics"
    def prefix_basename = prefix.substring(0, prefix.lastIndexOf("."))

    """
    touch "${prefix}"
    touch "${prefix}.crai"
    touch "${prefix}.bai"
    touch "${metrics}"
    touch "${metrics}.multiqc.tsv"
    touch "${prefix_basename}.score"
    """
}
