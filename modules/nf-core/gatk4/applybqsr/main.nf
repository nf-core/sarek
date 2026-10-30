process GATK4_APPLYBQSR {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b33544fcc9268dad01801b928e566f1beaacd078b1d58d7904d7c6aa68eba9a/data'
        : 'community.wave.seqera.io/library/gatk4-lite_samtools:9b6806e26dd32afe'}"

    input:
    tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    tuple val(meta2), path(fasta), path(fai), path(dict)
    val output_format

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: reads
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // default to cram output with crai index
    def file_type = output_format ? output_format[0] : 'cram'
    def index_format = output_format ? output_format[1] : 'crai'

    if ((file_type == 'bam' && index_format == 'crai') || (file_type == 'cram' && ['bai', 'csi'].contains(index_format))) {
        error("GATK4_APPLYBQSR: invalid index format '${index_format}' for '${file_type}' output.")
    }

    // GATK/htsjdk writes CRAM 3.1 by default (or 3.0), and only htsjdk 4.1+/5 can read those files back,
    // so GATK streams BAM to stdout and samtools does the cram conversion and indexing
    def samtools_view_args = file_type == 'cram' ? "--reference ${fasta} -C" : "-b"
    def output_file = "${prefix}.${file_type}##idx##${prefix}.${file_type}.${index_format}"

    def interval_command = intervals ? "--intervals ${intervals}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        --input ${input} \\
        --output /dev/stdout \\
        --reference ${fasta} \\
        --bqsr-recal-file ${bqsr_table} \\
        ${interval_command} \\
        -OBI false \\
        --tmp-dir . \\
        ${args} | samtools view --write-index ${samtools_view_args} -o ${output_file} ${args2} -
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    // default to cram output with crai index
    def file_type = output_format ? output_format[0] : 'cram'
    def index_format = output_format ? output_format[1] : 'crai'

    if ((file_type == 'bam' && index_format == 'crai') || (file_type == 'cram' && ['bai', 'csi'].contains(index_format))) {
        error("GATK4_APPLYBQSR: invalid index format '${index_format}' for '${file_type}' output.")
    }

    """
    touch ${prefix}.${file_type}
    touch ${prefix}.${file_type}.${index_format}
    """
}
