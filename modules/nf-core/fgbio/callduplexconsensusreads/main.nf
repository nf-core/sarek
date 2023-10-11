process FGBIO_CALLDUPLEXCONSENSUSREADS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::fgbio=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0' :
        'biocontainers/fgbio:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    // please note:
    // --min-reads is a required argument with no default
    // --min-input-base-quality is a required argument with no default
    // make sure they are specified via ext.args in your config

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_consensus"

    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio CallDuplexConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else if (mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            mem_gb = 1
        } else {
            mem_gb = task.memory.giga - 1
        }
    }

    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CallDuplexConsensusReads \\
        --input $bam \\
        --output ${prefix}.bam \\
        --threads ${task.cpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
