process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals), path(dragstr_model)
    path  fasta
    path  fai
    path  dict
    path  dbsnp
    path  dbsnp_tbi

    output:
    tuple val(meta), path("*.vcf.gz")       , emit: vcf
    tuple val(meta), path("*.tbi")          , optional:true, emit: tbi
    tuple val(meta), path("*.realigned.bam"), optional:true, emit: bam
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dbsnp_command = dbsnp ? "--dbsnp $dbsnp" : ""
    def interval_command = intervals ? "--intervals $intervals" : ""
    def dragstr_command = dragstr_model ? "--dragstr-params-path $dragstr_model" : ""
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        HaplotypeCaller \\
        --input $input \\
        --output ${prefix}.vcf.gz \\
        --reference $fasta \\
        $dbsnp_command \\
        $interval_command \\
        $dragstr_command \\
        $bamout_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bamout_command = args.contains("--bam-writer-type") ? "--bam-output ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""

    def stub_realigned_bam = bamout_command ? "touch ${prefix.replaceAll('.g\\s*$', '')}.realigned.bam" : ""
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    ${stub_realigned_bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
