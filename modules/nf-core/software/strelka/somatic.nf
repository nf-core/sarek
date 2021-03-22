include { initOptions; saveFiles; getSoftwareName } from './../functions'

params.options = [:]
def options    = initOptions(params.options)

process STRELKA_SOMATIC {
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
        tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
        path fasta
        path fai
        tuple path(target_bed), path(target_bed_tbi)

    output:
        tuple val(meta), path("*_somatic_indels.vcf.gz"), path("*_somatic_indels.vcf.gz.tbi"), emit: indels_vcf
        tuple val(meta), path("*_somatic_snvs.vcf.gz"), path("*_somatic_snvs.vcf.gz.tbi"),     emit: snvs_vcf
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
    configureStrelkaSomaticWorkflow.py \
        --tumor ${bam_tumor} \
        --normal ${bam_normal} \
        --referenceFasta ${fasta} \
        ${options_strelka} \
        --runDir strelka

    python strelka/runWorkflow.py -m local -j ${task.cpus}

    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}_somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}_somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}_somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}_somatic_snvs.vcf.gz.tbi

    echo configureStrelkaSomaticWorkflow.py --version &> ${software}.version.txt #2>&1
    """
}

process STRELKA_SOMATIC_BEST_PRACTICES {
    tag "${meta.id}"
    
    label 'CPUS_MAX'
    label 'MEMORY_MAX'
    
    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::strelka=2.9.10" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/strelka:2.9.10--0"
    } else {
        container "quay.io/biocontainers/strelka:2.9.10--0"
    }

    input:
        tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor), path(manta_csi), path(manta_csi_tbi)
        path fasta
        path fai
        tuple path(target_bed), path(target_bed_tbi)

    output:
        tuple val(meta), path("*_somatic_indels.vcf.gz"), path("*_somatic_indels.vcf.gz.tbi"), emit: indels_vcf
        tuple val(meta), path("*_somatic_snvs.vcf.gz"), path("*_somatic_snvs.vcf.gz.tbi"),     emit: snvs_vcf
        path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "strelka_bp_${meta.id}${ioptions.suffix}" : "strelka_bp_${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    options_strelka = params.target_bed ? "--exome --callRegions ${target_bed}" : ""
    """
    configureStrelkaSomaticWorkflow.py \
        --tumor ${bam_tumor} \
        --normal ${bam_normal} \
        --referenceFasta ${fasta} \
        --indelCandidates ${manta_csi} \
        ${options_strelka} \
        --runDir strelka

    python strelka/runWorkflow.py -m local -j ${task.cpus}

    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}_somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}_somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}_somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}_somatic_snvs.vcf.gz.tbi

    echo configureStrelkaSomaticWorkflow.py --version &> ${software}.version.txt #2>&1
    """
}