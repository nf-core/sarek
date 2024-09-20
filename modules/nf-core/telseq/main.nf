process TELSEQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5ce2a0c04652b0d0cc87f012a2240e1e5a90bc90:eb084e8aa92146f3987d00af8d38b36214d1f39f-0':
        'biocontainers/mulled-v2-5ce2a0c04652b0d0cc87f012a2240e1e5a90bc90:eb084e8aa92146f3987d00af8d38b36214d1f39f-0' }"

    input:
    tuple val(meta ), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(bed)

    output:
    tuple val(meta), path("*.telseq.tsv"), emit: output
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exome  = bed ? " --exomebed=${bed}" : ""
    """
    # telseq doesn't support CRAM. See https://github.com/zd1/telseq/issues/26
    if ${bam.name.endsWith(".cram")}
    then
        samtools view -T ${fasta} -O BAM --uncompressed ${bam} |\\
        telseq ${args} ${exome} - > tmp.tsv
    else
        telseq ${args} ${exome} ${bam} > tmp.tsv
    fi

    #
    # 'bug' in telseq, messages that should be printed on stderr are printed on stdout
    # We remove them with awk
    #
    awk '/^ReadGroup/ {ok=1;} {if(ok) print;}' tmp.tsv > ${prefix}.telseq.tsv
    rm tmp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telseq: \$(telseq --help 2>&1 | grep "^Version" -m1 | cut -d ' ' -f2)
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.telseq.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telseq: \$(telseq --help 2>&1 | grep "^Version" -m1 | cut -d ' ' -f2)
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
