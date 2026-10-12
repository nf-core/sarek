process LOFREQ_CALLPARALLEL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4'
        : 'quay.io/biocontainers/lofreq:2.1.5--py38h588ecb2_4'}"

    input:
    tuple val(meta), path(bam), path(bai), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"),     emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('lofreq'), eval("lofreq version | sed -n '1s/^.* //p'"), emit: versions_lofreq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_intervals = intervals ? "-l ${intervals}" : ""

    def alignment_cram = bam.Extension == "cram" ? true : false
    def alignment_out = alignment_cram ? bam.BaseName + ".bam" : "${bam}"

    def samtools_cram_convert = ''
    samtools_cram_convert += alignment_cram ? "    samtools view -T ${fasta} ${bam} -@ ${task.cpus} -o ${alignment_out}\n" : ''
    samtools_cram_convert += alignment_cram ? "    samtools index ${alignment_out}\n" : ''

    def samtools_cram_remove = ''
    samtools_cram_remove += alignment_cram ? "    rm ${alignment_out}\n" : ''
    samtools_cram_remove += alignment_cram ? "    rm ${alignment_out}.bai\n " : ''
    """
    ${samtools_cram_convert}

    lofreq \\
        call-parallel \\
        --pp-threads ${task.cpus} \\
        ${args} \\
        ${options_intervals} \\
        -f ${fasta} \\
        -o ${prefix}.vcf.gz \\
        ${alignment_out}

    ${samtools_cram_remove}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.vcf.gz.tbi
    """
}
