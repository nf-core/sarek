process GATK_BASERECALIBRATOR {
    label 'cpus_1'

    tag "${meta.id}-${interval.baseName}"

    input:
        tuple val(meta), path(bam), path(bai), path(interval)
        path dbsnp
        path dbsnp_tbi
        path dict
        path fai
        path fasta
        path known_indels
        path known_indels_tbi
        
    output:
        tuple val(meta), path("${prefix}${meta.sample}.recal.table"), emit: report
        val meta, emit: tsv

    //when: params.known_indels

    script:
    options_dbsnp = params.dbsnp ? "--known-sites ${dbsnp}" : ""
    options_intervals = params.no_intervals ? "" : "-L ${interval}"
    options_known_indels = params.known_indels ? known_indels.collect{"--known-sites ${it}"}.join(' ') : ""
    prefix = params.no_intervals ? "" : "${interval.baseName}_"
    // TODO: --use-original-qualities ???
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        -I ${bam} \
        -O ${prefix}${meta.sample}.recal.table \
        --tmp-dir . \
        -R ${fasta} \
        ${options_dbsnp} \
        ${options_known_indels} \
        ${options_intervals} \
        --verbosity INFO
    """
}