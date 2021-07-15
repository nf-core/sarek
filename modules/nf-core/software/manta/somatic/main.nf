// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MANTA_SOMATIC {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::manta=1.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/manta:1.6.0--py27_0"
    } else {
        container "quay.io/biocontainers/manta:1.6.0--py27_0"
    }

    input:
        tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
        path fasta
        path fai
        tuple path(target_bed), path(target_bed_tbi)

    output:
        tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor), path("*.candidateSmallIndels.vcf.gz"), path("*.candidateSmallIndels.vcf.gz.tbi"), emit: manta_csi_for_strelka_bp
        tuple val(meta), path("*.candidateSmallIndels.vcf.gz"), path("*.candidateSmallIndels.vcf.gz.tbi"), emit: candidate_small_indels_vcf
        tuple val(meta), path("*.candidateSV.vcf.gz"), path("*.candidateSV.vcf.gz.tbi"),                   emit: candidate_sv_vcf
        tuple val(meta), path("*.diploidSV.vcf.gz"), path("*.diploidSV.vcf.gz.tbi"),                       emit: diploid_sv_vcf
        tuple val(meta), path("*.somaticSV.vcf.gz"), path("*.somaticSV.vcf.gz.tbi"),                       emit: somatic_sv_vcf
        path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "manta_${meta.id}${ioptions.suffix}" : "manta_${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    options_manta = params.target_bed ? "--exome --callRegions $target_bed" : ""
    """
    configManta.py \
        --tumorBam $bam_tumor \
        --normalBam $bam_normal \
        --reference $fasta \
        $options_manta \
        --runDir manta

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz     ${prefix}.candidateSmallIndels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi ${prefix}.candidateSmallIndels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz              ${prefix}.candidateSV.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi          ${prefix}.candidateSV.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz                ${prefix}.diploidSV.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi            ${prefix}.diploidSV.vcf.gz.tbi
    mv manta/results/variants/somaticSV.vcf.gz                ${prefix}.somaticSV.vcf.gz
    mv manta/results/variants/somaticSV.vcf.gz.tbi            ${prefix}.somaticSV.vcf.gz.tbi

    echo configManta.py --version &> ${software}.version.txt #2>&1
    """
}