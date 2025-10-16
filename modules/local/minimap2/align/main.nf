process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::minimap2=2.26 bioconda::samtools=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val sort_bam
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample ?: meta.id}_${task.index}"
    
    // Extract read group - remove quotes if present
    def rg = meta.read_group ?: "@RG\\\\tID:${meta.id}\\\\tSM:${meta.sample}\\\\tLB:${meta.sample}\\\\tPL:ILLUMINA\\\\tPU:${meta.lane ?: '1'}"
    def read_group_string = rg.toString().replaceAll('^"|"$', '')
    
    def reference = index[0].toString().endsWith('.mmi') ? index[0] : fasta
    """
    set -euo pipefail
    
    echo "DEBUG: Starting minimap2 alignment for ${meta.id}"
    echo "DEBUG: Read group string: ${read_group_string}"
    
    # For Illumina short reads, use -x sr preset
    PRESET="-x sr"
    echo "DEBUG: Using short-read preset"
    
    # Find the appropriate index
    INDEX=\$( find -L ./ -name "*.mmi" | head -n 1 )
    if [ -z "\$INDEX" ]; then
        INDEX=${fasta}
        echo "DEBUG: No .mmi index found, using fasta: \$INDEX"
    else
        echo "DEBUG: Using .mmi index: \$INDEX"
    fi
    
    echo "DEBUG: Starting alignment and sorting..."
    minimap2 \\
        \$PRESET \\
        -a \\
        -t $task.cpus \\
        -R '${read_group_string}' \\
        $args \\
        \$INDEX \\
        ${reads[0]} \\
        ${reads[1]} \\
        | samtools sort \\
            $args2 \\
            -@ $task.cpus \\
            -m \$(echo $task.memory | sed 's/ GB/G/g' | sed 's/.GB/G/g' | awk '{print int(\$1/8)"G"}') \\
            -O BAM \\
            -T ${prefix}.tmp \\
            -o ${prefix}.bam -
    
    echo "DEBUG: Indexing BAM..."
    samtools index -@ $task.cpus ${prefix}.bam
    
    echo "DEBUG: Validating read groups..."
    samtools view -H ${prefix}.bam | grep '^@RG' || (echo "ERROR: No read groups found in BAM!" && exit 1)
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
