process STRELKA_GERMLINE {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1'
        : 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'}"

    input:
    tuple val(meta), path(input), path(input_index), path(target_bed), path(target_bed_index)
    path fasta
    path fai

    output:
    tuple val(meta), path("*variants.vcf.gz"),     emit: vcf
    tuple val(meta), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*genome.vcf.gz"),       emit: genome_vcf
    tuple val(meta), path("*genome.vcf.gz.tbi"),   emit: genome_vcf_tbi
    tuple val("${task.process}"), val('strelka'), eval("configureStrelkaGermlineWorkflow.py --version"), emit: versions_strelka, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = target_bed ? "--callRegions ${target_bed}" : ""
    """
    configureStrelkaGermlineWorkflow.py \\
        --bam ${input} \\
        --referenceFasta ${fasta} \\
        ${regions} \\
        ${args} \\
        --runDir strelka

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g strelka/runWorkflow.py

    python strelka/runWorkflow.py -m local -j ${task.cpus}
    mv strelka/results/variants/genome.*.vcf.gz     ${prefix}.genome.vcf.gz
    mv strelka/results/variants/genome.*.vcf.gz.tbi ${prefix}.genome.vcf.gz.tbi
    mv strelka/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
    mv strelka/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.genome.vcf.gz
    touch ${prefix}.genome.vcf.gz.tbi
    echo "" | gzip > ${prefix}.variants.vcf.gz
    touch ${prefix}.variants.vcf.gz.tbi
    """
}
