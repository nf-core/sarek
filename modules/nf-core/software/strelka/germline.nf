include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process STRELKA_GERMLINE {
    tag "$meta.id"
    label 'process_high'
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
    } else {
        container "quay.io/biocontainers/strelka:2.9.10--0"
    }

    input:
        tuple val(meta), path(bam), path(bai)
        path fasta
        path fai
        tuple path(target_bed), path(target_bed_tbi)

    output:
        tuple val(meta), path("*_variants.vcf.gz"), path("*_variants.vcf.gz.tbi"), emit: vcf
        tuple val(meta), path("*_genome.vcf.gz"), path("*_genome.vcf.gz.tbi"), emit: genome_vcf
        path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "strelka_${meta.id}${ioptions.suffix}" : "strelka_${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    options_strelka = params.target_bed ? "--exome --callRegions ${target_bed}" : ""
    """
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options_strelka} \
        --runDir strelka

    python strelka/runWorkflow.py -m local -j ${task.cpus}

    mv strelka/results/variants/genome.*.vcf.gz     ${prefix}_genome.vcf.gz
    mv strelka/results/variants/genome.*.vcf.gz.tbi ${prefix}_genome.vcf.gz.tbi
    mv strelka/results/variants/variants.vcf.gz     ${prefix}_variants.vcf.gz
    mv strelka/results/variants/variants.vcf.gz.tbi ${prefix}_variants.vcf.gz.tbi

    echo configureStrelkaGermlineWorkflow.py --version &> ${software}.version.txt #2>&1
    """
}