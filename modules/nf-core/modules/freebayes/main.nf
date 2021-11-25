// Import generic module functions
include { initOptions; saveFiles; getProcessName; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FREEBAYES {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::freebayes=1.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/freebayes:1.3.5--py38ha193a2f_3"
    } else {
        container "quay.io/biocontainers/freebayes:1.3.5--py38ha193a2f_3"
    }

    input:
    tuple val(meta), path(input_1), path(input_1_index), path(input_2), path(input_2_index)
    path fasta
    path fai
    path targets
    path samples
    path populations
    path cnv

    output:
    tuple val(meta), path("*.vcf.gz")   , emit: vcf
    path  "versions.yml"                , emit: versions

    script:
    def prefix           = options.suffix ? "${meta.id}${options.suffix}"  : "${meta.id}"
    def input            = input_2        ? "${input_1} ${input_2}"        : "${input_1}"
    def targets_file     = targets        ? "--target ${targets}"          : ""
    def samples_file     = samples        ? "--samples ${samples}"         : ""
    def populations_file = populations    ? "--populations ${populations}" : ""
    def cnv_file         = cnv            ? "--cnv-map ${cnv}"             : ""

    if (task.cpus > 1) {
        """
        freebayes-parallel \\
            <(fasta_generate_regions.py ${fasta}.fai 10000) ${task.cpus} \\
            -f $fasta \\
            $targets_file \\
            $samples_file \\
            $populations_file \\
            $cnv_file \\
            $options.args \\
            $input > ${prefix}.vcf

        gzip --no-name ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
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
            $options.args \\
            $input > ${prefix}.vcf

        gzip --no-name ${prefix}.vcf

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(freebayes --version 2>&1) | sed 's/version:\s*v//g' )
        END_VERSIONS
        """
    }
}
