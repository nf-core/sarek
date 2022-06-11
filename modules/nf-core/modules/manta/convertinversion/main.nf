process MANTA_CONVERTINVERSION {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1 bioconda::manta=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-40295ae41112676b05b649e513fe7000675e9b84:0b4be2c719f99f44df34be7b447b287bb7f86e01-0':
        'quay.io/biocontainers/mulled-v2-40295ae41112676b05b649e513fe7000675e9b84:0b4be2c719f99f44df34be7b447b287bb7f86e01-0' }"

    input:
    tuple val(meta), path(vcf)
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convertInversion.py \$(which samtools) $fasta $vcf | bgzip --threads $task.cpus > ${prefix}.vcf.gz
    tabix ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
