process GATK_HAPLOTYPECALLER {
    label 'MEMORY_SINGLECPU_TASK_SQ'
    label 'CPUS_2'

    tag "${meta.id}-${interval.baseName}"

    input:
        tuple val(meta), path(bam), path(bai), file(interval)
        path dbsnp
        path dbsnpIndex
        path dict
        path fasta
        path fai

    output:
        tuple val("HaplotypeCallerGVCF"), val(meta), path("${interval.baseName}_${meta.id}.g.vcf")  emit: gvcfHaplotypeCaller
        tuple val(meta), path(interval), path("${intervalBed.baseName}_${meta.id}.g.vcf")           emit: gvcfGenotypeGVCFs

   

    script:
    intervalsOptions = params.no_intervals ? "" : "-L ${interval}"
    dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -O ${interval.baseName}_${meta.id}.g.vcf \
        -ERC GVCF
    """
}
