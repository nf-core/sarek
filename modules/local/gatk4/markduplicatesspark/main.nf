process GATK4_MARKDUPLICATES_SPARK {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::gatk4==4.1.9.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
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
    tuple val(meta), path("*.${format}"), emit: output
    path "versions.yml"                 , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def bams = bam.collect(){ x -> "-I ".concat(x.toString()) }.join(" ")
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    if (!task.memory) {
        log.info '[GATK MarkDuplicatesSpark] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk  \
        MarkDuplicatesSpark \
        $bams \
        -O ${prefix}.${format} \
        --reference ${reference} \
        --tmp-dir . \
        --spark-master local[${task.cpus}] \\
        $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
