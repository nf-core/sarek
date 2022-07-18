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
    path("*.bed")       , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fasta_fai} > ${fasta_fai.baseName}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>/dev/null | head -n 1 | awk '{split(\$0,a,","); print a[1];}' | egrep -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
