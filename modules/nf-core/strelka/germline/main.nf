process STRELKA_GERMLINE {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    conda "bioconda::strelka=2.9.10"
    container { NfcoreTemplate.getContainer(
        docker: 'biocontainers/strelka:2.9.10--h9ee0642_1',
        singularity: 'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1',
        registry: 'quay.io',
        use_full_uri: task.ext.container_full_uri ?: false
    )}

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

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g strelka/runWorkflow.py

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
