process GATK4_FILTERMUTECTCALLS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi), path(stats), path(orientationbias), path(segmentation), path(table), val(estimate)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.filteringStats.tsv"), emit: stats
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def orientationbias_command = orientationbias ? orientationbias.collect { orientationbias_ -> "--orientation-bias-artifact-priors ${orientationbias_}" }.join(' ') : ''
    def segmentation_command = segmentation ? segmentation.collect { segmentation_ -> "--tumor-segmentation ${segmentation_}" }.join(' ') : ''
    def estimate_command = estimate ? " --contamination-estimate ${estimate} " : ''
    def table_command = table ? table.collect { table_ -> "--contamination-table ${table_}" }.join(' ') : ''

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK FilterMutectCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        FilterMutectCalls \\
        --variant ${vcf} \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        ${orientationbias_command} \\
        ${segmentation_command} \\
        ${estimate_command} \\
        ${table_command} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.filteringStats.tsv
    """
}
