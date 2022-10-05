process STRELKA_GERMLINE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1' }"

    input:
    tuple val(meta), path(input), path(input_index), path (target_bed), path (target_bed_tbi)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*variants.vcf.gz")    , emit: vcf
    tuple val(meta), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*genome.vcf.gz")      , emit: genome_vcf
    tuple val(meta), path("*genome.vcf.gz.tbi")  , emit: genome_vcf_tbi
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions  = target_bed ? "--callRegions ${target_bed}" : ""
    """
    configureStrelkaGermlineWorkflow.py \\
        --bam $input \\
        --referenceFasta $fasta \\
        $regions \\
        $args \\
        --runDir strelka

    python strelka/runWorkflow.py -m local -j $task.cpus
    mv strelka/results/variants/genome.*.vcf.gz     ${prefix}.genome.vcf.gz
    mv strelka/results/variants/genome.*.vcf.gz.tbi ${prefix}.genome.vcf.gz.tbi
    mv strelka/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
    mv strelka/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaGermlineWorkflow.py --version )
    END_VERSIONS
    """
}
