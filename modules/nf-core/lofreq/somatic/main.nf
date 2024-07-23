process LOFREQ_SOMATIC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py38h588ecb2_4' :
        'biocontainers/lofreq:2.1.5--py38h588ecb2_4' }"
    input:
    tuple val(meta), path(tumor), path(tumor_index), path(normal), path(normal_index), path(target_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*somatic_final.snvs.vcf.gz")      , emit: vcf_snvs
    tuple val(meta), path("*somatic_final.snvs.vcf.gz.tbi")  , emit: vcf_snvs_tbi
    tuple val(meta), path("*somatic_final.indels.vcf.gz")    , emit: vcf_indels
    tuple val(meta), path("*somatic_final.indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}_"
    def options_target_bed = target_bed ? "-l ${target_bed}" : ""

    def tumor_cram =  tumor.Extension == "cram" ? true : false
    def normal_cram =  normal.Extension == "cram" ? true : false

    def tumor_out = tumor_cram ? tumor.BaseName + ".bam" : "${tumor}"
    def normal_out = normal_cram ? normal.BaseName + ".bam" : "${normal}"

    def samtools_cram_convert = ''
    samtools_cram_convert += normal_cram ? "    samtools view -T $fasta $normal -@ $task.cpus -o $normal_out\n" : ''
    samtools_cram_convert += normal_cram ? "    samtools index $normal_out\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools view -T ${fasta} $tumor -@ $task.cpus -o $tumor_out\n" : ''
    samtools_cram_convert += tumor_cram ? "    samtools index ${tumor_out}\n" : ''

    def samtools_cram_remove = ''
    samtools_cram_remove += tumor_cram ? "     rm $tumor_out\n" : ''
    samtools_cram_remove += tumor_cram ? "     rm ${tumor_out}.bai\n " : ''
    samtools_cram_remove += normal_cram ? "     rm $normal_out\n" : ''
    samtools_cram_remove += normal_cram ? "     rm ${normal_out}.bai\n " : ''
    """
    $samtools_cram_convert

    lofreq \\
        somatic \\
        --threads $task.cpus \\
        $args \\
        -f $fasta \\
        -t $tumor_out \\
        -n $normal_out \\
        ${options_target_bed} \\
        -o ${prefix}

    $samtools_cram_remove

    if [ -f "*_minus-dbsnp.snvs.vcf.gz" ] && [ -f "*_minus-dbsnp.indels.vcf.gz" ]; then

        mv *somatic_final.snvs.vcf.gz  ${prefix}somatic.snvs.vcf.gz 
        mv *somatic_final.indels.vcf.gz ${prefix}somatic.indels.vcf.gz

        mv *_minus-dbsnp.snvs.vcf.gz ${prefix}somatic_final.snvs.vcf.gz 
        mv *_minus-dbsnp.indels.vcf.gz ${prefix}somatic_final.indels.vcf.gz 
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_"
    """
    echo "" | gzip > ${prefix}somatic_final.snvs.vcf.gz
    echo "" | gzip > ${prefix}somatic_final.snvs.vcf.gz.tbi
    echo "" | gzip > ${prefix}somatic_final.indels.vcf.gz
    echo "" | gzip > ${prefix}somatic_final.indels.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(echo \$(lofreq version 2>&1) | sed 's/^version: //; s/ *commit.*\$//')
    END_VERSIONS
    """
}