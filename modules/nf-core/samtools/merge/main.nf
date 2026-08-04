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

    output:
    tuple val(meta), path("${prefix}.bam"), optional: true, emit: bam
    tuple val(meta), path("${prefix}.cram"), optional: true, emit: cram
    tuple val(meta), path("*.{bai,crai,csi}"), optional: true, emit: index
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    def reference = fasta ? "--reference ${fasta}" : ""
    // samtools merge --write-index defaults to .csi for BAM; force the legacy .bai/.crai extensions via the ##idx## naming trick so index files match what downstream sample-sheet validation expects
    def index_type = file_type == "bam" ? "bai" : "crai"
    def output = args.contains("--write-index") ? "${prefix}.${file_type}##idx##${prefix}.${file_type}.${index_type}" : "${prefix}.${file_type}"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        merge \\
        --threads ${task.cpus - 1} \\
        ${args} \\
        ${reference} \\
        ${output} \\
        ${input_files}
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input_files instanceof List ? input_files[0].getExtension() : input_files.getExtension()
    def index_type = file_type == "bam" ? "bai" : "crai"
    def index = args.contains("--write-index") ? "touch ${prefix}.${file_type}.${index_type}" : ""
    """
    touch ${prefix}.${file_type}
    ${index}
    """
}
