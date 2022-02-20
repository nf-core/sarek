process FREEBAYES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py38ha193a2f_3' :
        'quay.io/biocontainers/freebayes:1.3.5--py38ha193a2f_3' }"

    input:
    tuple val(meta), path(input_1), path(input_1_index), path(input_2), path(input_2_index), path(target_bed)
    path fasta
    path fasta_fai
    path samples
    path populations
    path cnv

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input            = input_2        ? "${input_1} ${input_2}"        : "${input_1}"
    def targets_file     = target_bed     ? "--target ${target_bed}"       : ""
    def samples_file     = samples        ? "--samples ${samples}"         : ""
    def populations_file = populations    ? "--populations ${populations}" : ""
    def cnv_file         = cnv            ? "--cnv-map ${cnv}"             : ""

    if (task.cpus > 1) {
        """
        freebayes-parallel \\
            <(fasta_generate_regions.py $fasta_fai 10000) $task.cpus \\
            -f $fasta \\
            $targets_file \\
            $samples_file \\
            $populations_file \\
            $cnv_file \\
            $args \\
            $input > ${prefix}.vcf

        bgzip ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
        END_VERSIONS
        """

    } else {
        """
        freebayes \\
            -f $fasta \\
            $targets_file \\
            $samples_file \\
            $populations_file \\
            $cnv_file \\
            $args \\
            $input > ${prefix}.vcf

        bgzip ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            freebayes: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
        END_VERSIONS
        """
    }
}
