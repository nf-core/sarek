process CAT_CAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.8' :
        'quay.io/biocontainers/pigz:2.8' }"

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    tuple val("${task.process}"), val("pigz"), eval("pigz --version 2>&1 | sed 's/pigz //g'"),  topic: versions, emit: versions_cat

    when:
    task.ext.when == null || task.ext.when

    script:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/find/concatenate

    Reason:
    This module passes all input files as shell arguments, which can exceed the UNIX ARG_MAX
    limit when concatenating large numbers of files. The find/concatenate module resolves this
    by staging files into a directory and enumerating them with `find`, and also enforces
    consistent input compression (all gzipped or all uncompressed). Also enables faster, parallel
    decompression with pigz.
    """
    assert false: deprecation_message

    def file_list = files_in.collect { file -> file.toString() }
    prefix   = task.ext.prefix ?: "${meta.id}${getFileSuffix(file_list[0])}"

    stub:
    def deprecation_message = """
    WARNING: This module has been deprecated. Please use nf-core/modules/find/concatenate

    Reason:
    This module passes all input files as shell arguments, which can exceed the UNIX ARG_MAX
    limit when concatenating large numbers of files. The find/concatenate module resolves this
    by staging files into a directory and enumerating them with `find`, and also enforces
    consistent input compression (all gzipped or all uncompressed). Also enables faster, parallel
    decompression with pigz.
    """
    assert false: deprecation_message
}

// for .gz files also include the second to last extension if it is present. E.g., .fasta.gz
def getFileSuffix(filename) {
    def match = filename =~ /^.*?((\.\w{1,5})?(\.\w{1,5}\.gz$))/
    return match ? match[0][1] : filename.substring(filename.lastIndexOf('.'))
}
