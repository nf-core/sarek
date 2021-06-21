include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process GATK4_MARKDUPLICATES_SPARK {
    label 'process_high'
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4==4.1.9.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.1.9.0--py39_0"
    } else {
        container "broadinstitute/gatk:4.2.0.0"
    }

    input:
        tuple val(meta), path(bam)
        path(reference)
        path(dict) //need to be present in the path
        path(fai)  //need to be present in the path
        val(format) //either "bam" or "cram"

    output:
        tuple val(meta), path('*.${format}'), emit: output
        path("*.version.txt"),           emit: version

    script:
    def software = getSoftwareName(task.process)
    def bams = bam.collect(){ x -> "-I ".concat(x.toString()) }.join(" ")
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """

    gatk  \
        MarkDuplicatesSpark \
        ${bams} \
        -O ${prefix}.${format} \
        --reference ${reference} \
        --tmp-dir . \
        --spark-master local[${task.cpus}] \\
        $options.args

    echo \$(gatk MarkDuplicatesSpark --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}