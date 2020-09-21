process BWA_MEM {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/bwa/${meta.sample}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"
    //container "https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:eabfac3657eda5818bae4090db989e3d41b01542-0"

    conda (params.conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.10" : null)

    input:
        tuple val(meta), path(reads)
        path bwa
        path fasta
        path fai
        val options

    output:
        tuple val(meta), path("*.bam"), emit: bam
        path  "*.version.txt"         , emit: version

    script:
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${meta.run}\\t${CN}PU:${meta.run}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA"
    extra = meta.status == 1 ? "-B 3" : ""
    """
    bwa mem \
        ${options.args} \
        -R \"${readGroup}\" \
        ${extra} \
        -t ${task.cpus} \
        ${fasta} ${reads} | \
    samtools sort --threads ${task.cpus} -m 2G - > ${meta.id}.bam

    # samtools index ${meta.id}.bam

    echo \$(bwa version 2>&1) > bwa.version.txt
    """
}
