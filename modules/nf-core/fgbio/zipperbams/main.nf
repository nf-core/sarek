process FGBIO_ZIPPERBAMS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::fgbio=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0' :
        'biocontainers/fgbio:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(unmapped_bam)
    tuple val(meta), path(mapped_bam)
    path(fasta)
    path(dict)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def compression = task.ext.compression ?: '0'
    prefix = task.ext.prefix ?: "${meta.id}_zipped"
    def fgbio_mem_gb = 4

    if (!task.memory) {
        log.info '[fgbio ZipperBams] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else if (fgbio_mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            fgbio_mem_gb = 1
        } else {
            fgbio_mem_gb = task.memory.giga - 1
        }
    }

    """
    fgbio -Xmx${fgbio_mem_gb}g \\
        --compression ${compression} \\
        --async-io=true \\
        ZipperBams \\
        --unmapped ${unmapped_bam} \\
        --input ${mapped_bam} \\
        --ref ${fasta} \\
        ${args} \\
        --output ${prefix}.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
