process XENGSORT_CLASSIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/xengsort:2.1.0--pyhdfd78af_1':
        'community.wave.seqera.io/library/htslib_pip_xengsort:39f2d076bf655602'}"

    input:
    tuple val(meta), path(reads)
    path  index_folder

    output:
    tuple val(meta), path('*-graft.*fq.gz')    , optional:true, emit: graft_fastq
    tuple val(meta), path('*-host.*fq.gz')     , optional:true, emit: host_fastq
    tuple val(meta), path('*-both.*fq.gz')     , optional:true, emit: both_fastq
    tuple val(meta), path('*-neither.*fq.gz')  , optional:true, emit: neither_fastq
    tuple val(meta), path('*-ambiguous.*fq.gz'), optional:true, emit: ambiguous_fastq
    tuple val("${task.process}"), val('xengsort'), eval("xengsort --version"), topic: versions, emit: versions_xengsort_classify

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mode = task.ext.mode ?: 'count'

    def fastq_args = meta.single_end ? "--fastq ${reads}" : "--fastq ${reads[0]} --pairs ${reads[1]}"

    """
    # Check if the index folder contains the expected index files
    HASH_COUNT=\$(ls "${index_folder}"/*.hash 2>/dev/null | wc -l || true)
    INFO_COUNT=\$(ls "${index_folder}"/*.info 2>/dev/null | wc -l || true)

    if [ "\$HASH_COUNT" -ne 1 ] || [ "\$INFO_COUNT" -ne 1 ]; then
        echo "ERROR: The input index directory must contain exactly one .hash and one .info file." >&2
        echo "Found \$HASH_COUNT .hash files and \$INFO_COUNT .info files in ${index_folder}." >&2
        exit 1
    fi

    # Extract the index prefix (the .hash and .info have the same basename)
    HASH_FILE=\$(ls "${index_folder}"/*.hash)
    INDEX_PREFIX=\$(basename "\$HASH_FILE" .hash)

    # Run xengsort classify
    xengsort classify \\
        --index "${index_folder}/\$INDEX_PREFIX" \\
        $fastq_args \\
        --mode $mode \\
        --prefix $prefix \\
        --threads ${task.cpus} \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
        """
        echo $args

        echo '' | gzip > ${prefix}-graft.fq.gz
        echo '' | gzip > ${prefix}-host.fq.gz
        echo '' | gzip > ${prefix}-both.fq.gz
        echo '' | gzip > ${prefix}-neither.fq.gz
        echo '' | gzip > ${prefix}-ambiguous.fq.gz
        """
    } else {
        """
        echo $args

        echo '' | gzip > ${prefix}-graft.1.fq.gz
        echo '' | gzip > ${prefix}-graft.2.fq.gz
        echo '' | gzip > ${prefix}-host.1.fq.gz
        echo '' | gzip > ${prefix}-host.2.fq.gz
        echo '' | gzip > ${prefix}-both.1.fq.gz
        echo '' | gzip > ${prefix}-both.2.fq.gz
        echo '' | gzip > ${prefix}-neither.1.fq.gz
        echo '' | gzip > ${prefix}-neither.2.fq.gz
        echo '' | gzip > ${prefix}-ambiguous.1.fq.gz
        echo '' | gzip > ${prefix}-ambiguous.2.fq.gz
        """
    }
}
