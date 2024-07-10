process NARFMAP_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/narfmap_align:8bad41386eab9997':
        'community.wave.seqera.io/library/narfmap_align:517a1fed8e4e84c1' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(hashmap)
    tuple val(meta3), path(fasta)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path('*.log'), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command       = meta.single_end ? "-1 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def samtools_command    = sort_bam ? 'sort' : 'view'
    def extension_pattern   = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher   =  (args2 =~ extension_pattern)
    def extension           = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference           = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    """
    dragen-os \\
        -r $hashmap \\
        $args \\
        --num-threads $task.cpus \\
        $reads_command \\
        2> >(tee ${prefix}.narfmap.log >&2) \\
        | samtools $samtools_command $args2 --threads $task.cpus ${reference} -o ${prefix}.${extension} -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        narfmap: \$(echo \$(dragen-os --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_command = meta.single_end ? "-1 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension_pattern = /(--output-fmt|-O)+\s+(\S+)/
    def extension_matcher =  (args2 =~ extension_pattern)
    def extension = extension_matcher.getCount() > 0 ? extension_matcher[0][2].toLowerCase() : "bam"
    def reference = fasta && extension=="cram"  ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"

    def create_index = ""
    if (extension == "cram") {
        create_index = "touch ${prefix}.crai"
    } else if (extension == "bam") {
        create_index = "touch ${prefix}.csi"
    }

    """
    touch ${prefix}.${extension}
    ${create_index}
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        narfmap: \$(echo \$(dragen-os --version 2>&1))
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
