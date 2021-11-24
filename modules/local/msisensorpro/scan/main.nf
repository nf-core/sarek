// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MSISENSORPRO_SCAN {
    tag "$fasta"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::msisensor-pro=1.1.a" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/msisensor-pro:1.1.a--hb3646a4_0"
    } else {
        container "quay.io/biocontainers/msisensor-pro:1.1.a--hb3646a4_0"
    }

    input:
    path fasta

    output:
    tuple val(meta), path("*.list"), emit: list
    path "versions.yml"            , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${options.suffix}" : ""
    """
    msisensor-pro scan \\
        -d $fasta \\
        -o ${fasta.baseName}.list \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(msisensor 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p')
    END_VERSIONS
    """
}
