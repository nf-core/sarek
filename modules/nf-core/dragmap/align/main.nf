process DRAGMAP_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::dragmap=1.2.1 bioconda::samtools=1.15.1 conda-forge::pigz=2.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:5ebebbc128cd624282eaa37d2c7fe01505a91a69-0':
        'quay.io/biocontainers/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:5ebebbc128cd624282eaa37d2c7fe01505a91a69-0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(hashmap)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path('*.log'), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-1 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def samtools_command = sort_bam ? 'sort' : 'view'

    """
    dragen-os \\
        -r $hashmap \\
        $args \\
        --num-threads $task.cpus \\
        $reads_command \\
        2> ${prefix}.dragmap.log \\
        | samtools $samtools_command $args2 --threads $task.cpus -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
