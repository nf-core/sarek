// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process STRELKA {
    tag "$meta.id"
    
    label 'cpus_max'
    label 'memory_max'
    
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id)
    //     }

    container "quay.io/biocontainers/strelka:2.9.10--0"
    //conda (params.conda ? "bioconda::strelka=2.9.10" : null)

    input:
    tuple val(meta), path(bam), path (bai)
    path(fasta)
    path(fai)
    path(target_bed)
    val options

    output:
    tuple val("Strelka"), val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcfStrelkaSingle
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$ioptions.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${target_bed} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? ioptions.args : ""
    """
    ${beforeScript}
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options} \
        --runDir Strelka
    python Strelka/runWorkflow.py -m local -j ${task.cpus}
    mv Strelka/results/variants/genome.*.vcf.gz \
        Strelka_${meta.sample}_genome.vcf.gz
    mv Strelka/results/variants/genome.*.vcf.gz.tbi \
        Strelka_${meta.sample}_genome.vcf.gz.tbi
    mv Strelka/results/variants/variants.vcf.gz \
        Strelka_${meta.sample}_variants.vcf.gz
    mv Strelka/results/variants/variants.vcf.gz.tbi \
        Strelka_${meta.sample}_variants.vcf.gz.tbi

    echo configureStrelkaGermlineWorkflow.py --version &> ${software}.version.txt #2>&1 || true
    #echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}