// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MSISENSORPRO_MSI {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::msisensor-pro=1.1.a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/msisensor-pro:1.1.a--hb3646a4_0"
    } else {
        container "quay.io/biocontainers/msisensor-pro:1.1.a--hb3646a4_0"
    }

    input:
    tuple val(meta), path(cram_normal), path(crai_normal), path(cram_tumor), path(crai_tumor)
    path msisensorpro_scan

    output:
    tuple val(meta), path("${prefix}.list")         , emit: output
    tuple val(meta), path("${prefix}_dis.list")     , emit: output_dis
    tuple val(meta), path("${prefix}_germline.list"), emit: output_germline
    tuple val(meta), path("${prefix}_somatic.list") , emit: output_somatic
    path "versions.yml"                             , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    msisensor-pro msi \\
        -d $msisensorpro_scan \\
        -n $bam_normal \\
        -t $bam_tumor \\
        -o $prefix \\
        -b $task.cpus \\
        $options.args

    mv ${prefix}          ${prefix}.list
    mv ${prefix}_dis      ${prefix}_dis.list
    mv ${prefix}_germline ${prefix}_germline.list
    mv ${prefix}_somatic  ${prefix}_somatic.list

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(msisensor 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
