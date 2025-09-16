process VALIDATE_READGROUPS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.validated.bam"), path("*.validated.bai"), emit: bam
    path "*.validation_report.txt"                                    , emit: report
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # First, check current read groups
    echo "Current read groups in BAM:" > ${prefix}.validation_report.txt
    samtools view -H ${bam} | grep '^@RG' >> ${prefix}.validation_report.txt || echo "No read groups found!" >> ${prefix}.validation_report.txt

    # Validate with GATK's ValidateSamFile
    gatk ValidateSamFile \\
        -I ${bam} \\
        -R ${fasta} \\
        -MODE SUMMARY \\
        --VALIDATE_INDEX true \\
        --MAX_OUTPUT 100 \\
        --IGNORE_WARNINGS false \\
        -O ${prefix}.gatk_validation.txt || true

    echo "\nGATK Validation Results:" >> ${prefix}.validation_report.txt
    cat ${prefix}.gatk_validation.txt >> ${prefix}.validation_report.txt

    # Check for missing required read group fields
    RG_LINE=\$(samtools view -H ${bam} | grep '^@RG' | head -1)
    MISSING_FIELDS=""
    
    for field in ID SM LB PL PU; do
        if ! echo "\$RG_LINE" | grep -q "\${field}:"; then
            MISSING_FIELDS="\$MISSING_FIELDS \$field"
        fi
    done

    if [ -n "\$MISSING_FIELDS" ]; then
        echo "\nWARNING: Missing required read group fields:\$MISSING_FIELDS" >> ${prefix}.validation_report.txt
        echo "Adding missing fields..." >> ${prefix}.validation_report.txt
        
        # Use GATK's AddOrReplaceReadGroups to fix
        gatk AddOrReplaceReadGroups \\
            -I ${bam} \\
            -O ${prefix}.validated.bam \\
            -ID ${meta.id}_fixed \\
            -SM ${meta.sample ?: meta.patient ?: meta.id} \\
            -LB ${meta.library ?: meta.id} \\
            -PL ${meta.platform ?: 'ILLUMINA'} \\
            -PU ${meta.platform_unit ?: 'UNKNOWN.1.0'} \\
            --CREATE_INDEX true \\
            --SORT_ORDER coordinate \\
            --VALIDATION_STRINGENCY LENIENT
        
        echo "Read groups fixed!" >> ${prefix}.validation_report.txt
    else
        echo "\nRead groups are GATK-compliant!" >> ${prefix}.validation_report.txt
        # Just copy the files
        cp ${bam} ${prefix}.validated.bam
        cp ${bai} ${prefix}.validated.bai
    fi

    # Final validation
    echo "\nFinal read groups:" >> ${prefix}.validation_report.txt
    samtools view -H ${prefix}.validated.bam | grep '^@RG' >> ${prefix}.validation_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
