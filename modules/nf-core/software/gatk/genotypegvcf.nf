process GATK_GENOTYPEVCF {
    tag "${meta.id}-${interval.baseName}"

    container "quay.io/biocontainers/gatk4-spark:4.1.8.1--0"

    conda (params.conda ? "bioconda::gatk4-spark=4.1.8.1" : null)

    input:
        tuple val(meta), path(interval), path(gvcf)
        path dbsnp
        path dbsnpIndex
        path dict
        path fasta
        path fai

    output:
     tuple val("HaplotypeCaller"), val(meta), path("${interval.baseName}_${meta.id}.vcf")

    script:
    // Using -L is important for speed and we have to index the interval files also
    intervalsOptions = params.no_intervals ? "" : "-L ${interval}"
    dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        IndexFeatureFile \
        -I ${gvcf}

    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -V ${gvcf} \
        -O ${interval.baseName}_${meta.id}.vcf
    """
}