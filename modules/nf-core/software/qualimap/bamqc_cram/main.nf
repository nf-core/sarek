// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process QUALIMAP_BAMQC_CRAM{ // Import generic module functions {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::qualimap=2.2.2d bioconda::samtools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-d3934ca6bb4e61334891ffa2e9a4c87a530e3188:4bf11d12f2c3eccf1eb585097c0b6fd31c18c418-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-d3934ca6bb4e61334891ffa2e9a4c87a530e3188:4bf11d12f2c3eccf1eb585097c0b6fd31c18c418-0"
    }

    input:
    tuple val(meta), path(cram), path(crai)
    path gff
    val use_gff
    path reference
    path fai

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "*.version.txt"             , emit: version

    script:
    def software   = getSoftwareName(task.process)
    prefix         = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'
    def memory     = task.memory.toGiga() + "G"
    def regions = use_gff ? "--gff $gff" : ''

    def strandedness = 'non-strand-specific'
    if (meta.strandedness == 'forward') {
        strandedness = 'strand-specific-forward'
    } else if (meta.strandedness == 'reverse') {
        strandedness = 'strand-specific-reverse'
    }
    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp

    samtools view -hb -T ${reference} ${cram} |
    qualimap \\
        --java-mem-size=$memory \\
        bamqc \\
        $options.args \\
        -bam /dev/stdin \\
        $regions \\
        -p $strandedness \\
        $collect_pairs \\
        -outdir $prefix \\
        -nt $task.cpus

    echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//' > ${software}.version.txt
    """
}
