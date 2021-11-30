process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'
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
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bam", includeInputs:true), path("*.bai")  , optional:true, emit: bam_bai
    tuple val(meta), path("*.bam", includeInputs:true), path("*.csi")  , optional:true, emit: bam_csi
    tuple val(meta), path("*.cram", includeInputs:true), path("*.crai"), optional:true, emit: cram_crai
    path  "versions.yml"                                                              , emit: version

    script:
    """
    samtools index -@ ${task.cpus-1} $options.args $input

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
