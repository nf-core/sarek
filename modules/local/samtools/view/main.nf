//There is a -L option to only output alignments in interval, might be an option for exons/panel data?
process SAMTOOLS_VIEW {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::samtools=1.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.14--hb421002_0"
    } else {
        container "quay.io/biocontainers/samtools:1.14--hb421002_0"
    }

    input:
    tuple val(meta), path(input), path(index)
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bam"), path("*.bai")  , optional: true, emit: bam_bai
    tuple val(meta), path("*.cram"), path("*.crai"), optional: true, emit: cram_crai
    path  "versions.yml"                                           , emit: versions

    script:
    def prefix    = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def reference = fasta ? "--reference ${fasta} -C" : ""
    """
    samtools view --threads ${task.cpus-1} ${reference} $options.args $input > ${prefix}.cram
    samtools index -@${task.cpus} ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
