process SAMTOOLS_MERGE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(input_files, stageAs: "?/*"), path(index_files, stageAs: "?/*")
    tuple val(meta2), path(fasta), path(fai), path(gzi)
    val index_format

    output:
    tuple val(meta), path("${prefix}.bam"), optional: true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional: true, emit: cram
    tuple val(meta), path("${prefix}.${file_type}.{bai,crai,csi}"), optional: true, emit: index
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    //setting default values
    def write_index = ""
    def output_file = "${prefix}.${file_type}"

    // Update if index is requested
    if (index_format != '' && index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error("Index format not one of bai, csi, crai.")
        }
        write_index = "--write-index"
        output_file = "${prefix}.${file_type}##idx##${prefix}.${file_type}.${index_format}"
    }
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        merge \\
        --threads ${task.cpus - 1} \\
        ${write_index} \\
        ${args} \\
        ${reference} \\
        ${output_file} \\
        ${input_files}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()

    if (index_format) {
        if (!index_format.matches('bai|csi|crai')) {
            error("Index format not one of bai, csi, crai.")
        }
    }

    def index = index_format ? "touch ${prefix}.${file_type}.${index_format}" : ""

    """
    touch ${prefix}.${file_type}
    ${index}
    """
}
