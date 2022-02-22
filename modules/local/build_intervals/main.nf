process BUILD_INTERVALS {
    tag "$fasta_fai"
    label 'process_medium'

    conda (params.enable_conda ? "anaconda::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    path fasta_fai

    output:
    path "*.bed", emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fasta_fai} > ${fasta_fai.baseName}.bed
    """
}
