process BWA_MEM {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0' :
        'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0' }"

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    script:
    def split_cpus = Math.floor(task.cpus/2)
    def args       = task.ext.args  ?: ''
    def args2      = task.ext.args2 ?: ''
    def part       = params.split_fastq > 1 ? reads.get(0).name.findAll(/part_([0-9]+)?/).last().concat('.') : ""
    def prefix     = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""    //MD Spark NEEDS name sorted reads or runtime goes through the roof.
    //However, if duplicate marking is skipped, reads need to be coordinate sorted.
    //Spark can be used also for BQSR, therefore check for both: only name sort if spark + duplicate marking is done
    def sort_order = ('markduplicates' in params.use_gatk_spark) & !params.skip_markduplicates ? "-n" : ""
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
        $args \\
        $read_group \\
        -t ${split_cpus} \\
        \$INDEX \\
        $reads \\
        | samtools $args2 $sort_order --threads ${split_cpus} -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
