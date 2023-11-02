process SAMTOOLS_INDEX_FOR_INDEXCOV {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val(meta), path(input), path(bai)
    path(fasta)
    path(fai)

    output:
    tuple path("${meta.id}.indexcov.bam"),path("${meta.id}.indexcov.bam.bai"), emit: output
    path  "versions.yml"                 , emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    # write BAM header only
    samtools view --header-only -O BAM  \
        --threads ${task.cpus} \
        -o "${meta.id}.indexcov.bam" \
        --reference "${fasta}" \
		"${input}"

    # create index without writing BAM (redirecting to /dev/null)
    samtools view ${args} --uncompressed \
		--threads ${task.cpus} \
		-o "/dev/null##idx##${meta.id}.indexcov.bam.bai" \
		--write-index \
		-O BAM  \
	        --reference "${fasta}" \
		"${input}"



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.indexcov.bam"
    touch "${meta.id}.indexcov.bam.bai"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
