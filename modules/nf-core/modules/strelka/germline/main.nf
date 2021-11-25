// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRELKA_GERMLINE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
    } else {
        container "quay.io/biocontainers/strelka:2.9.10--0"
    }

    input:
    tuple val(meta), path(input), path(input_index)
    path  fasta
    path  fai
    path  target_bed
    path  target_bed_tbi


    output:
    tuple val(meta), path("*variants.vcf.gz")    , emit: vcf
    tuple val(meta), path("*variants.vcf.gz.tbi"), emit: vcf_tbi
    tuple val(meta), path("*genome.vcf.gz")      , emit: genome_vcf
    tuple val(meta), path("*genome.vcf.gz.tbi")  , emit: genome_vcf_tbi
    path "versions.yml"                          , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def regions  = target_bed ? "--exome --callRegions ${target_bed}" : ""
    """
    configureStrelkaGermlineWorkflow.py \\
        --bam $input \\
        --referenceFasta $fasta \\
        $regions \\
        $options.args \\
        --runDir strelka

    python strelka/runWorkflow.py -m local -j $task.cpus
    mv strelka/results/variants/genome.*.vcf.gz     ${prefix}.genome.vcf.gz
    mv strelka/results/variants/genome.*.vcf.gz.tbi ${prefix}.genome.vcf.gz.tbi
    mv strelka/results/variants/variants.vcf.gz     ${prefix}.variants.vcf.gz
    mv strelka/results/variants/variants.vcf.gz.tbi ${prefix}.variants.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( configureStrelkaGermlineWorkflow.py --version )
    END_VERSIONS
    """
}
