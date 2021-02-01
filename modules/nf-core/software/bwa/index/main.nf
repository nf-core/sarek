include { initOptions; saveFiles; getSoftwareName } from './../../functions'

params.options = [:]
def options    = initOptions(params.options)

process BWA_INDEX {
    label 'process_high'

    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"false") }

    conda (params.enable_conda ? "bioconda::bwa=0.7.17" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7"
    } else {
        container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
    }

    input:
        path fasta

    output:
        path "${fasta}.*"   , emit: index
        path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    bwa index ${ioptions.args} ${fasta}
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    """
}
