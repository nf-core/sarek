process GOLEFT_INDEXCOV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/goleft:0.2.4--h9ee0642_1'
        : 'quay.io/biocontainers/goleft:0.2.4--h9ee0642_1'}"

    input:
    tuple val(meta), path(bams), path(indexes)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("${prefix}/*"),           emit: output
    tuple val(meta), path("${prefix}/*ped"),        emit: ped,       optional: true
    tuple val(meta), path("${prefix}/*bed.gz"),     emit: bed,       optional: true
    tuple val(meta), path("${prefix}/*bed.gz.tbi"), emit: bed_index, optional: true
    tuple val(meta), path("${prefix}/*roc"),        emit: roc,       optional: true
    tuple val(meta), path("${prefix}/*html"),       emit: html,      optional: true
    tuple val(meta), path("${prefix}/*png"),        emit: png,       optional: true
    tuple val("${task.process}"), val('goleft'), eval("goleft --version |& sed '1!d;s/^.*goleft Version: //'"), topic: versions, emit: versions_goleft
    tuple val("${task.process}"), val('tabix'), eval("tabix -h |& sed -n 's/^.*Version: //p'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // indexcov uses BAM files or CRAI
    def input_files = bams.findAll {bam_file -> bam_file.name.endsWith(".bam") } + indexes.findAll {index_file -> index_file.name.endsWith(".crai") }
    def extranormalize = input_files.any {input_file -> input_file.name.endsWith(".crai") } ? " --extranormalize " : ""
    """
    goleft indexcov \\
        --fai ${fai}  \\
        --directory ${prefix} \\
        ${extranormalize} \\
        ${args} \\
        ${input_files.join(" ")}

    if [ -f "${prefix}/${prefix}-indexcov.bed.gz" ] ; then
        tabix -p bed ${prefix}/${prefix}-indexcov.bed.gz
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    echo "" | gzip > ${prefix}/${prefix}-indexcov.bed.gz
    touch ${prefix}/${prefix}-indexcov.bed.gz.tbi
    """
}
