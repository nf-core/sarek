include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.1--he513fc3_0 bioconda::samtools=1.11--h6270b1f_0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:e6f0d20c9d78572ddbbf00d8767ee6ff865edd4e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:e6f0d20c9d78572ddbbf00d8767ee6ff865edd4e-0"
    }

    input:
    tuple val(meta), path(reads)
    path  index
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "*.version.txt"         , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    """
    bwa-mem2 mem \
        $options.args \
        $read_group \
        -t $task.cpus \
        $fasta \
        $reads \
        | samtools $options.args2 --threads $task.cpus -o ${prefix}.bam -

    echo \$(bwa-mem2 version 2>&1) > ${software}.version.txt
    """
}