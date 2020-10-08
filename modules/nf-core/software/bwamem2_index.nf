include { initOptions; saveFiles; getSoftwareName } from './functions'

process BWAMEM2_INDEX {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/bwa-mem2:2.0--he513fc3_1"

    conda (params.conda ? "bioconda::bwa-mem2=2.0" : null)

    input:
        path fasta
        val options

    output:
        path "${fasta}.*"

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    bwa-mem2 index ${ioptions.args} ${fasta}

    echo \$(bwa-mem2 version 2>&1) > ${software}.version.txt
    """
}
