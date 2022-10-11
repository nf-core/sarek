process MANTA_GERMLINE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::manta=1.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1' :
        'quay.io/biocontainers/manta:1.6.0--h9ee0642_1' }"

    input:
    //Matching the target bed with the input sample allows to parallelize the same sample run across different intervals or a single bed file
    tuple val(meta), path(input), path(index), path(target_bed), path(target_bed_tbi)
    path fasta
    path fasta_fai

    output:
    tuple val(meta), path("*candidate_small_indels.vcf.gz")    , emit: candidate_small_indels_vcf
    tuple val(meta), path("*candidate_small_indels.vcf.gz.tbi"), emit: candidate_small_indels_vcf_tbi
    tuple val(meta), path("*candidate_sv.vcf.gz")              , emit: candidate_sv_vcf
    tuple val(meta), path("*candidate_sv.vcf.gz.tbi")          , emit: candidate_sv_vcf_tbi
    tuple val(meta), path("*diploid_sv.vcf.gz")                , emit: diploid_sv_vcf
    tuple val(meta), path("*diploid_sv.vcf.gz.tbi")            , emit: diploid_sv_vcf_tbi
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_files = input.collect{"--bam ${it}"}.join(' ')
    def options_manta = target_bed ? "--callRegions $target_bed" : ""
    """
    configManta.py \
        ${input_files} \
        --reference $fasta \
        --runDir manta \
        $options_manta \
        $args

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz \
        ${prefix}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \
        ${prefix}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \
        ${prefix}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz \
        ${prefix}.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi \
        ${prefix}.diploid_sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.candidate_small_indels.vcf.gz
    touch ${prefix}.candidate_small_indels.vcf.gz.tbi
    touch ${prefix}.candidate_sv.vcf.gz
    touch ${prefix}.candidate_sv.vcf.gz.tbi
    touch ${prefix}.diploid_sv.vcf.gz
    touch ${prefix}.diploid_sv.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS
    """
}
