process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mosdepth=0.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.3--hdfd78af_1' :
        'quay.io/biocontainers/mosdepth:0.3.3--hdfd78af_1'}"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(bed)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path('*.global.dist.txt')      , emit: global_txt
    tuple val(meta), path('*.summary.txt')          , emit: summary_txt
    tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regions_txt
    tuple val(meta), path('*.per-base.d4')          , optional:true, emit: per_base_d4
    tuple val(meta), path('*.per-base.bed.gz')      , optional:true, emit: per_base_bed
    tuple val(meta), path('*.per-base.bed.gz.csi')  , optional:true, emit: per_base_csi
    tuple val(meta), path('*.regions.bed.gz')       , optional:true, emit: regions_bed
    tuple val(meta), path('*.regions.bed.gz.csi')   , optional:true, emit: regions_csi
    tuple val(meta), path('*.quantized.bed.gz')     , optional:true, emit: quantized_bed
    tuple val(meta), path('*.quantized.bed.gz.csi') , optional:true, emit: quantized_csi
    tuple val(meta), path('*.thresholds.bed.gz')    , optional:true, emit: thresholds_bed
    tuple val(meta), path('*.thresholds.bed.gz.csi'), optional:true, emit: thresholds_csi
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        exit 1, "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        exit 1, "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    mosdepth \\
        --threads $task.cpus \\
        $interval \\
        $reference \\
        $args \\
        $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.global.dist.txt
    touch ${prefix}.region.dist.txt
    touch ${prefix}.summary.txt
    touch ${prefix}.per-base.d4
    touch ${prefix}.per-base.bed.gz
    touch ${prefix}.per-base.bed.gz.csi
    touch ${prefix}.regions.bed.gz
    touch ${prefix}.regions.bed.gz.csi
    touch ${prefix}.quantized.bed.gz
    touch ${prefix}.quantized.bed.gz.csi
    touch ${prefix}.thresholds.bed.gz
    touch ${prefix}.thresholds.bed.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    END_VERSIONS
    """
}
