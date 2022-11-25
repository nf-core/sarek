process GATK4_CNNSCOREVARIANTS {
    tag "$meta.id"
    label 'process_low'

    //Conda is not supported at the moment: https://github.com/broadinstitute/gatk/issues/7811
    if (params.enable_conda) {
        exit 1, "Conda environments cannot be used for GATK4/CNNScoreVariants at the moment. Please use docker or singularity containers."
    }
    container "broadinstitute/gatk:4.3.0.0" //Biocontainers is missing a package

    input:
    tuple val(meta), path(vcf), path(tbi), path(aligned_input), path(intervals)
    path fasta
    path fai
    path dict
    path architecture
    path weights

    output:
    tuple val(meta), path("*cnn.vcf.gz")    , emit: vcf
    tuple val(meta), path("*cnn.vcf.gz.tbi"), emit: tbi
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligned_input = aligned_input ? "--input $aligned_input" : ""
    def interval_command = intervals ? "--intervals $intervals" : ""
    def architecture = architecture ? "--architecture $architecture" : ""
    def weights = weights ? "--weights $weights" : ""

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK CnnScoreVariants] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" CNNScoreVariants \\
        --variant $vcf \\
        --output ${prefix}.cnn.vcf.gz \\
        --reference $fasta \\
        $interval_command \\
        $aligned_input \\
        $architecture \\
        $weights \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
