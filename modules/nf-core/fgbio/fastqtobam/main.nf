process FGBIO_FASTQTOBAM {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4d/4d1150a2e123f49f8c268f0ab429847afae642376fa52af713b846b084df4a9f/data'
        : 'community.wave.seqera.io/library/fgbio:3.1.2--6e9400d507a9dc55'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: bam
    tuple val("${task.process}"), val('fgbio'), eval('fgbio --version 2>&1 | tr -d "[:cntrl:]" | sed -e "s/^.*Version: //;s/\\[.*$//"'), topic: versions, emit: versions_fgbio

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bam"
    def sample_name = args.contains("--sample") ? "" : "--sample ${prefix}"
    def library_name = args.contains("--library") ? "" : "--library ${prefix}"

    def mem_gb = 8
    if (!task.memory) {
        log.info('[fgbio FastqToBam] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.')
    }
    else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        }
        else {
            mem_gb = task.memory.giga - 1
        }
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        FastqToBam \\
        ${args} \\
        --input ${reads} \\
        --output ${prefix}.${suffix} \\
        ${sample_name} \\
        ${library_name}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = task.ext.suffix ?: "bam"
    """
    touch ${prefix}.${suffix}
    """
}
