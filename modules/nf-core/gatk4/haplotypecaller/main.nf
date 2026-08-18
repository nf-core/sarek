process GATK4_HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    tuple val(meta5), path(dbsnp)
    tuple val(meta6), path(dbsnp_tbi)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.realigned.bam"), emit: bam, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbsnp_command = dbsnp ? "--dbsnp ${dbsnp}" : ""
    def interval_command = intervals ? "--intervals ${intervals}" : ""
    def dragstr_command = dragstr_model ? "--dragstr-params-path ${dragstr_model}" : ""
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        HaplotypeCaller \\
        --input ${input} \\
        --output ${prefix}.vcf.gz \\
        --reference ${fasta} \\
        --native-pair-hmm-threads ${task.cpus} \\
        ${dbsnp_command} \\
        ${interval_command} \\
        ${dragstr_command} \\
        ${bamout_command} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${stub_realigned_bam}
    """
}
