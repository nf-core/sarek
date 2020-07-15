params.bwa_options = "-M -B 2"
params.sequencer = "ILLUMINA"

process BWAMEM2_MEM {
    tag {id}

    publishDir "${params.outdir}/bwamem2_mem", mode: 'copy'

    input:
        tuple val(id), path(reads)
        path genomeindex
        val indexprefix

    output:
        tuple path("*.bam"), path("*.bai")

    script:
    CN = params.sequencing_center ? "\\tCN:${params.sequencing_center}\\t" : ""
    """
    bwa-mem2 mem -t ${task.cpus} -R "@RG\\tID:${id}${CN}\\tLB:${id}\\tSM:${id}\\tPL:${params.sequencer}" \\
    ${params.bwa_options} ${indexprefix} ${reads} | samtools sort -@8 -O BAM -o ${id}.bam -
    samtools index ${id}.bam
    """
}