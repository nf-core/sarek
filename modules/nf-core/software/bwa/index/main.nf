include { initOptions; saveFiles; getSoftwareName } from './../../functions'

process BWA_INDEX {
    tag "${fasta}"

    label 'process_high'

    publishDir params.outdir,
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"

    conda (params.conda ? "bioconda::bwa=0.7.17" : null)

    input:
        path fasta
        val options

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
