// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRELKA_SOMATIC {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1"
    } else {
        container "quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
    }

    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor),  path(manta_candidate_small_indels), path(manta_candidate_small_indels_tbi)
    path  fasta
    path  fai
    path  target_bed
    path  target_bed_tbi

    output:
    tuple val(meta), path("*.somatic_indels.vcf.gz")    , emit: vcf_indels
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    tuple val(meta), path("*.somatic_snvs.vcf.gz")      , emit: vcf_snvs
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi")  , emit: vcf_snvs_tbi
    path "versions.yml"                                 , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def options_target_bed = target_bed ? "--exome --callRegions ${target_bed}" : ""
    def options_manta = manta_candidate_small_indels ? "--indelCandidates ${manta_candidate_small_indels}" : ""
    """
    configureStrelkaSomaticWorkflow.py \\
        --tumor $input_tumor \\
        --normal $input_normal \\
        --referenceFasta $fasta \\
        $options_target_bed \\
        $options_manta \\
        $options.args \\
        --runDir strelka

    python strelka/runWorkflow.py -m local -j $task.cpus

    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """
}
