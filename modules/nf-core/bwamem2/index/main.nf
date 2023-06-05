process BWAMEM2_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::bwa-mem2=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--he513fc3_0' :
        'biocontainers/bwa-mem2:2.2.1--he513fc3_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwamem2"), emit: index
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bwamem2
    bwa-mem2 \\
        index \\
        $args \\
        $fasta -p bwamem2/${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir bwamem2
    touch bwamem2/${fasta}.0123
    touch bwamem2/${fasta}.ann
    touch bwamem2/${fasta}.pac
    touch bwamem2/${fasta}.amb
    touch bwamem2/${fasta}.bwt.2bit.64

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
