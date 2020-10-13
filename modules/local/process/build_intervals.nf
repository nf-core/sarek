include { initOptions; saveFiles; getSoftwareName } from './../../nf-core/software/functions'

environment = params.conda ? "anaconda::gawk=5.1.0" : null
container = "quay.io/biocontainers/gawk:5.1.0"
if (workflow.containerEngine == 'singularity') container = "https://depot.galaxyproject.org/singularity/gawk:5.1.0"

process BUILD_INTERVALS {
    tag fai

    publishDir params.outdir, mode: params.publish_dir_mode,
    saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    conda environment
    container container

    input:
        path fai

    output:
        path "${fai.baseName}.bed"

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fai} > ${fai.baseName}.bed
    """
}