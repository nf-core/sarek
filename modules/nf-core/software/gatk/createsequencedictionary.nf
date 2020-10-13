include { initOptions; saveFiles; getSoftwareName } from './../functions'

environment = params.conda ? "bioconda::gatk4-spark=4.1.8.1" : null
container = "quay.io/biocontainers/gatk4-spark:4.1.8.1--0"
if (workflow.containerEngine == 'singularity') container = "https://depot.galaxyproject.org/singularity/gatk4-spark:4.1.8.1--0"

process GATK_CREATESEQUENCEDICTIONARY {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda environment
    container container

    input:
        path fasta
        val options

    output:
        path "${fasta.baseName}.dict"

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict

    echo \$(gatk CreateSequenceDictionary --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}