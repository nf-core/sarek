process SNPEFF_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2' :
        'biocontainers/snpeff:5.1--hdfd78af_2' }"

    input:
    tuple val(meta), val(snpeff_db)

    output:
    tuple val(meta), path('snpeff_cache'), emit: cache
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6144
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    snpEff \\
        -Xmx${avail_mem}M \\
        download ${snpeff_db} \\
        -dataDir \${PWD}/snpeff_cache \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p snpeff_cache/${snpeff_db}

    touch snpeff_cache/${snpeff_db}/sequence.I.bin
    touch snpeff_cache/${snpeff_db}/sequence.bin

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
