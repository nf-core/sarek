process TMB {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::tmb=1.3.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tmb:1.3.0--pyh5e36f6f_0':
        'quay.io/biocontainers/tmb:1.3.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(vcf), path (dbconfig), path (varconfig)
    path (target_bed)

    output:
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def target_bed = target_bed ? "--bed ${target_bed}" : ""
    """
    pyTMB.py -i $vcf \
        --dbConfig ${dbconfig} \
        --varConfig ${varconfig} \
        ${target_bed} \
        --sample TUMOR
        $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tmb: \$(echo \$(pyTMB.py --version 2>&1) | sed 's/^.*pyTMB.py //; s/.*\$//' | sed 's|[()]||g')
    END_VERSIONS
    """
}