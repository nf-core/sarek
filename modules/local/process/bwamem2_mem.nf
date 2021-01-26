include { initOptions; saveFiles; getSoftwareName } from './../../nf-core/software/functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::bwa-mem2=2.0 bioconda::samtools=1.10" : null
container = "quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:876eb6f1d38fbf578296ea94e5aede4e317939e7-0"
if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) container = "https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:876eb6f1d38fbf578296ea94e5aede4e317939e7-0"

process BWAMEM2_MEM {
    label 'process_high'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda environment
    container container

    input:
        tuple val(meta), path(reads)
        path bwa
        path fasta
        path fai

    output:
        tuple val(meta), path("*.bam")

    script:
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${meta.run}\\t${CN}PU:${meta.run}\\tSM:${meta.sample}\\tLB:${meta.sample}\\tPL:ILLUMINA"
    extra = meta.status == 1 ? "-B 3" : ""
    """
    bwa-mem2 mem \
        ${options.args} \
        -R \"${readGroup}\" \
        ${extra} \
        -t ${task.cpus} \
        ${fasta} ${reads} | \
    samtools sort --threads ${task.cpus} -m 2G - > ${meta.id}.bam

    # samtools index ${meta.id}.bam

    echo \$(bwa-mem2 version 2>&1) > bwa-mem2.version.txt
    """
}