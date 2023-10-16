process SAMTOOLS_INDEX_FOR_INDEXCOV {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    output:
    tuple val(meta), path("*.bam"),  emit: bam
    tuple val(meta), path("*.bai") , emit: bai
    path  "versions.yml"           , emit: versions

    when:
    def args = task.ext.args ?: ''

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    # write BAM header only
    samtools view --header-only -O BAM  \
        --threads ${task.cpus} \
        -o "${prefix}.bam" \
        ${reference} \
		"${input}"

    # create index without writing BAM (redirecting to /dev/null)
	samtools view ${args} --uncompressed \
		--threads ${task.cpus} \
		-o "/dev/null##idx##${prefix}.bam.bai" \
		--write-index \
		-O BAM  \
		${reference} \
		"${input}"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch "${prefix}.bam"
    touch "${prefix}.bai"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
