include { initOptions; saveFiles; getSoftwareName } from './../functions'

process STRELKA_GERMLINE {
    tag "$meta.id"
    
    label 'CPUS_MAX'
    label 'MEMORY_MAX'
    
    publishDir "${params.outdir}",
         mode: params.publish_dir_mode,
         saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "quay.io/biocontainers/strelka:2.9.10--0"
    
    conda (params.conda ? "bioconda::strelka=2.9.10" : null)

    input:
    tuple val(meta), path(bam), path (bai)
    path fasta
    path fai
    path target_bed
    val options

    output:
    tuple val(meta), path("*_variants.vcf.gz"), path("*_variants.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*_genome.vcf.gz"), path("*_genome.vcf.gz.tbi"), emit: genome_vcf
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "Strelka_${meta.id}" : "Strelka_${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${target_bed} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options_strelka = params.target_bed ? ioptions.args : ""
    """
    ${beforeScript}
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options_strelka} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/genome.*.vcf.gz ${prefix}_genome.vcf.gz

    mv Strelka/results/variants/genome.*.vcf.gz.tbi ${prefix}_genome.vcf.gz.tbi

    mv Strelka/results/variants/variants.vcf.gz ${prefix}_variants.vcf.gz

    mv Strelka/results/variants/variants.vcf.gz.tbi ${prefix}_variants.vcf.gz.tbi

    echo configureStrelkaGermlineWorkflow.py --version &> ${software}.version.txt #2>&1
    """
}