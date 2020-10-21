include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::adam=0.32.0" : null
container = "quay.io/biocontainers/adam:0.32.0--0"
if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) container = "https://depot.galaxyproject.org/singularity/adam:0.32.0--0"

process ADAM_TRANSFORMALIGNMENTS {
    label 'cpus_16'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda environment
    container container

    input:
        tuple val(meta), path("${meta.sample}.bam"), path("${meta.sample}.bam.bai")

    output:
        tuple val(meta), path("${meta.sample}.md.bam"), emit: bam
        val meta,                                       emit: tsv
          
    script:
    """
    adam-submit \
      --master local[${task.cpus}] \
      --driver-memory ${task.memory.toGiga()} \
      -- \
      transformAlignments \
      -mark_duplicate_reads \
      -single \
      ${meta.sample}.bam \
      ${meta.sample}.md.bam
    """
}