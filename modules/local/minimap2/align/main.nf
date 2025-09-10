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
    path fasta
    path fai

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_group = meta.read_group ?: "${meta.id}.${meta.lane}"
    
    // Construct GATK-compliant read group string
    // Required fields: ID, LB, PL, SM, PU
    def rg_id = meta.read_group ?: "${meta.id}_T${task.index}"
    def rg_sm = meta.sample ?: meta.patient ?: meta.id  // Sample name (required)
    def rg_lb = meta.library ?: meta.id                  // Library (required)
    def rg_pl = meta.platform ?: "ILLUMINA"              // Platform (required)
    def rg_pu = meta.platform_unit ?: "${meta.flowcell ?: 'UNKNOWN'}.${meta.lane ?: '1'}.${meta.index ?: '0'}" // Platform unit (required)
    def rg_cn = meta.sequencing_center ?: ""             // Sequencing center (optional)
    def rg_dt = meta.date ?: ""                          // Date (optional)
    def rg_pm = meta.platform_model ?: ""                // Platform model (optional)
    
    // Build read group string with all required and optional fields
    def read_group_string = "@RG\\\\tID:${rg_id}\\\\tSM:${rg_sm}\\\\tLB:${rg_lb}\\\\tPL:${rg_pl}\\\\tPU:${rg_pu}"
    if (rg_cn) read_group_string += "\\\\tCN:${rg_cn}"
    if (rg_dt) read_group_string += "\\\\tDT:${rg_dt}"
    if (rg_pm) read_group_string += "\\\\tPM:${rg_pm}"
    
    def reference = index[0].toString().endsWith('.mmi') ? index[0] : fasta
    """
    # Set minimap2 preset based on read length
    READ_LENGTH=\$(gunzip -c ${reads[0]} | head -n 4000 | awk 'NR%4==2 {sum+=length(\$0); count++} END {print int(sum/count)}')
    
    if [[ \$READ_LENGTH -gt 1000 ]]; then
        PRESET="-x map-ont"  # For long reads
    else
        PRESET="-x sr"       # For short reads
    fi

    # Find the appropriate index
    INDEX=\$( find -L ./ -name "*.mmi" | head -n 1 )
    if [ -z "\$INDEX" ]; then
        # If no .mmi index exists, use the reference fasta
        INDEX=${fasta}
    fi

    # Align with minimap2 and pipe to samtools for sorting
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

    # Index the BAM file (required for GATK)
    samtools index -@ $task.cpus ${prefix}.bam

    # Validate the BAM file has proper read groups
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
