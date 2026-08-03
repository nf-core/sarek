process SAMTOOLS_BAM2FQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/31/315d2445cd42b0f5512fa37965a9c59bc93ae8614b7d105150caece6c61e2e71/data'
        : 'community.wave.seqera.io/library/htslib_samtools_xz:1595ae0727655963'}"

    input:
    tuple val(meta), path(inputbam)
    val split

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val("bgzip"), eval('bgzip --version | head -1 | sed "s/bgzip (htslib) //"'), emit: versions_bgzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (split) {
        """
        samtools \\
            bam2fq \\
            ${args} \\
            -@ ${task.cpus} \\
            -1 ${prefix}_1.fq.gz \\
            -2 ${prefix}_2.fq.gz \\
            -0 ${prefix}_other.fq.gz \\
            -s ${prefix}_singleton.fq.gz \\
            ${inputbam}
        """
    }
    else {
        """
        samtools \\
            bam2fq \\
            ${args} \\
            -@ ${task.cpus} \\
            ${inputbam} | bgzip > ${prefix}_interleaved.fq.gz
        """
    }

    stub:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def bgzip_command_1 = "echo | bgzip -c > ${prefix}_1.fq.gz"
    def bgzip_command_2 = "echo | bgzip -c > ${prefix}_2.fq.gz"
    def bgzip_command_other = "echo | bgzip -c > ${prefix}_other.fq.gz"
    def bgzip_command_singleton = "echo | bgzip -c > ${prefix}_singleton.fq.gz"
    def bgzip_command_interleaved = "echo | bgzip -c > ${prefix}_interleaved.fq.gz"

    if (split) {
        """
        ${bgzip_command_1}
        ${bgzip_command_2}
        ${bgzip_command_other}
        ${bgzip_command_singleton}
        """
    }
    else {
        """
        ${bgzip_command_interleaved}
        """
    }
}
