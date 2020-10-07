include { initOptions; saveFiles; getSoftwareName } from './../functions'

process SAMTOOLS_FAIDX {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
        path fasta
        val options

    output:
        path "${fasta}.fai"

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    samtools faidx ${fasta}

    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/ Using.*\$//' > ${software}.version.txt
    """
}
