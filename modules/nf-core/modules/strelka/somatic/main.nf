process STRELKA_SOMATIC {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1' }"

    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor),  path(manta_candidate_small_indels), path(manta_candidate_small_indels_tbi), path(target_bed), path(target_bed_index)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.somatic_indels.vcf.gz")    , emit: vcf_indels
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    tuple val(meta), path("*.somatic_snvs.vcf.gz")      , emit: vcf_snvs
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi")  , emit: vcf_snvs_tbi
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_target_bed = target_bed ? "--callRegions ${target_bed}" : ""
    def options_manta = manta_candidate_small_indels ? "--indelCandidates ${manta_candidate_small_indels}" : ""
    """

    configureStrelkaSomaticWorkflow.py \\
        --tumor $input_tumor \\
        --normal $input_normal \\
        --referenceFasta $fasta \\
        ${options_target_bed} \\
        ${options_manta} \\
        $args \\
        --runDir strelka

    python strelka/runWorkflow.py -m local -j $task.cpus

    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """
}
