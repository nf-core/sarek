// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_ESTIMATELIBRARYCOMPLEXITY {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gatk4:4.2.0.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/gatk4:4.2.0.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(bam)
    path(reference)
    path(dict) //need to be present in the path
    path(fai)  //need to be present in the path

    output:
    path('*.md.metrics'), emit: report
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def bams = bam.collect(){ x -> "-I ".concat(x.toString()) }.join(" ")
    """
    gatk EstimateLibraryComplexity \
        ${bams} \
        -O ${prefix}.metrics \
        --REFERENCE_SEQUENCE ${reference} \
        --VALIDATION_STRINGENCY SILENT \
        --TMP_DIR . $options.args

    echo \$(gatk EstimateLibraryComplexity --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
