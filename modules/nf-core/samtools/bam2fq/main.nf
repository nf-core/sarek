process SAMTOOLS_BAM2FQ {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
        'quay.io/biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(meta), path(inputbam)
    val split

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (split){
        """
        samtools \\
            bam2fq \\
            $args \\
            -@ $task.cpus \\
            -1 ${prefix}_1.fq.gz \\
            -2 ${prefix}_2.fq.gz \\
            -0 ${prefix}_other.fq.gz \\
            -s ${prefix}_singleton.fq.gz \\
            $inputbam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else {
        """
        samtools \\
            bam2fq \\
            $args \\
            -@ $task.cpus \\
            $inputbam | gzip --no-name > ${prefix}_interleaved.fq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
