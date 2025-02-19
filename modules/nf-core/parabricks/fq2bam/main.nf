process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.4.0-1"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)
    tuple val(meta4), path(interval_file)
    path(known_sites)

    output:
    tuple val(meta), path("*.bam")  , emit: bam
    tuple val(meta), path("*.bai")  , emit: bai
    tuple val(meta), path("*.table"), emit: bqsr_table         , optional:true
    path("versions.yml")            , emit: versions
    path("qc_metrics")              , emit: qc_metrics          , optional:true
    path("duplicate-metrics.txt")   , emit: duplicate_metrics  , optional:true

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def in_fq_command = meta.single_end ? "--in-se-fq $reads" : "--in-fq $reads"
    def known_sites_command = known_sites ? known_sites.collect{"--knownSites $it"}.join(' ') : ""
    def known_sites_output = known_sites ? "--out-recal-file ${prefix}.table" : ""
    def interval_file_command = interval_file ? interval_file.collect{"--interval-file $it"}.join(' ') : ""
    def num_gpus = task.accelerator ? "--num-gpus $task.accelerator.request" : ''
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    cp $fasta \$INDEX

    pbrun \\
        fq2bam \\
        --ref \$INDEX \\
        $in_fq_command \\
        --out-bam ${prefix}.bam \\
        $known_sites_command \\
        $known_sites_output \\
        $interval_file_command \\
        $num_gpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
            pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //' )
    END_VERSIONS
    """
}
