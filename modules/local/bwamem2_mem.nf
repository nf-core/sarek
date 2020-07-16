params.bwa_options = "-K 100000000 -M"
params.sequencer = "ILLUMINA"

process BWAMEM2_MEM {
    label 'CPUS_MAX'

    tag "${sample}_${run}"

    publishDir "${params.outdir}/bwamem2_mem", mode: 'copy'

    input:
        tuple val(patient), val(sample), val(run), path(read1), path(read2)
        path bwa
        path fasta
        path fai

    output:
        tuple val(patient), val(sample), val(run), path("*.bam"), path("*.bai")

    script:
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${run}\\t${CN}PU:${run}\\tSM:${sample}\\tLB:${sample}\\tPL:${params.sequencer}"
    """
    bwa-mem2 mem ${params.bwa_options} -R \"${readGroup}\" -t ${task.cpus} \
    ${fasta} ${read1} ${read2} | \
    samtools sort --threads ${task.cpus} -m 2G - > ${sample}_${run}.bam
    samtools index ${sample}_${run}.bam
    """
}