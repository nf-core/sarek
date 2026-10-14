process XENGSORT_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xengsort:2.2.1--pyhdfd78af_0':
        'quay.io/biocontainers/xengsort:2.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(hash), path(info)
    val mode

    output:
    tuple val(meta), path("*-host*.fq.gz"),                                     emit: host, optional: true
    tuple val(meta), path("*-graft*.fq.gz"),                                    emit: graft, optional: true
    tuple val(meta), path("*-both*.fq.gz"),                                     emit: both, optional: true
    tuple val(meta), path("*-neither*.fq.gz"),                                  emit: neither, optional: true
    tuple val(meta), path("*-ambiguous*.fq.gz"),                                emit: ambiguous, optional: true
    tuple val("${task.process}"), val('xengsort'), eval("xengsort --version"),  emit: versions_xengsort, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "--fastq ${reads}" : "--fastq ${reads[0]} --pairs ${reads[1]}"
    def index_base = hash.baseName
    """
    xengsort classify \
        --index ${index_base} \
        ${input_reads} \
        --out ${prefix} \
        --mode ${mode} \
        --compression gz \
        --threads ${task.cpus} \
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        echo "" | gzip > ${prefix}-host.fq.gz
        echo "" | gzip > ${prefix}-graft.fq.gz
        echo "" | gzip > ${prefix}-both.fq.gz
        echo "" | gzip > ${prefix}-neither.fq.gz
        echo "" | gzip > ${prefix}-ambiguous.fq.gz
        """
    }
    else {
        """
        echo "" | gzip > ${prefix}-host.1.fq.gz
        echo "" | gzip > ${prefix}-host.2.fq.gz
        echo "" | gzip > ${prefix}-graft.1.fq.gz
        echo "" | gzip > ${prefix}-graft.2.fq.gz
        echo "" | gzip > ${prefix}-both.1.fq.gz
        echo "" | gzip > ${prefix}-both.2.fq.gz
        echo "" | gzip > ${prefix}-neither.1.fq.gz
        echo "" | gzip > ${prefix}-neither.2.fq.gz
        echo "" | gzip > ${prefix}-ambiguous.1.fq.gz
        echo "" | gzip > ${prefix}-ambiguous.2.fq.gz
        """
    }
}
