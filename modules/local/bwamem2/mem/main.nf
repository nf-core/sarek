// Import generic module functions
 include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:cf603b12db30ec91daa04ba45a8ee0f35bbcd1e2-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:cf603b12db30ec91daa04ba45a8ee0f35bbcd1e2-0"
    }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def split_cpus = Math.floor(task.cpus/2)
    def software   = getSoftwareName(task.process)
    def part       = params.split_fastq > 1 ? reads.get(0).name.findAll(/part_([0-9]+)?/).last().concat('.') : ""
    def prefix     = options.suffix ? "${meta.id}${options.suffix}.${part}" : "${meta.id}.${part}"
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    //MD Spark NEEDS name sorted reads or runtime goes through the roof.
    //However, if duplicate marking is skipped, reads need to be coordinate sorted.
    //Spark can be used also for BQSR, therefore check for both: only name sort if spark + duplicate marking is done
    def sort_order = ('markduplicates' in params.use_gatk_spark) & !params.skip_markduplicates ? "-n" : ""
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa-mem2 mem \\
        $options.args \\
        $read_group \\
        -t ${split_cpus} \\
        \$INDEX \\
        $reads \\
        | samtools $options.args2 $sort_order -@ ${split_cpus} -o ${prefix}bam -

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
