process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a3/a3dc5ea2ce788c24079d24d1721ed28086874152c43b5e7dde3f638dcf64336a/data' :
        'community.wave.seqera.io/library/htslib_mosdepth_gzip:4108dd38be84e40a'}"

    input:
    tuple val(meta),  path(bam), path(bai), path(bed)
    tuple val(meta2), path(fasta)
    val(quantize_labels)

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
    tuple val("${task.process}"), val('mosdepth'), eval("mosdepth --version | sed 's/mosdepth //g'"), topic: versions, emit: versions_mosdepth
    tuple val("${task.process}"), val('gzip'), eval("gzip -V 2>&1 | sed 's/gzip \\([0-9.]*\\).*/\\1/;q'"), topic: versions, emit: versions_gzip
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    def interval = bed ? "--by ${bed}" : ""
    def quantize_env_vars = []
    if (quantize_labels instanceof List && quantize_labels.size() > 0) {
        quantize_labels.eachWithIndex { label, index ->
            quantize_env_vars << "MOSDEPTH_Q${index}=${label}"
        }
    }
    if (bed && (args.contains("--by") || args.contains("-b "))) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (args.contains("--thresholds") && !(bed || args.contains("--by") || args.contains("-b "))) {
        error "'--thresholds' can only be specified in conjunction with '--by' or an input bed file"
    }

    """
    ${quantize_env_vars.join(" ")} mosdepth \\
        --threads $task.cpus \\
        $interval \\
        $reference \\
        $args \\
        $prefix \\
        $bam
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (bed && (args.contains("--by") || args.contains("-b "))) {
        error "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (args.contains("--thresholds") && !(bed || args.contains("--by") || args.contains("-b "))) {
        error "'--thresholds' can only be specified in conjunction with '--by' or an input bed file"
    }
    """
    touch ${prefix}.global.dist.txt
    touch ${prefix}.region.dist.txt
    touch ${prefix}.summary.txt
    touch ${prefix}.per-base.d4
    echo "" | gzip > ${prefix}.per-base.bed.gz
    touch ${prefix}.per-base.bed.gz.csi
    echo "" | gzip > ${prefix}.regions.bed.gz
    touch ${prefix}.regions.bed.gz.csi
    echo "" | gzip > ${prefix}.quantized.bed.gz
    touch ${prefix}.quantized.bed.gz.csi
    echo "" | gzip > ${prefix}.thresholds.bed.gz
    touch ${prefix}.thresholds.bed.gz.csi
    """
}
