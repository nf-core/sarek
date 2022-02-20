process GATK4_MUTECT2 {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta) , path(input) , path(input_index) , path(intervals), val(which_norm)
    val  run_single
    val  run_pon
    val  run_mito
    path fasta
    path fai
    path dict
    path germline_resource
    path germline_resource_tbi
    path panel_of_normals
    path panel_of_normals_tbi

    output:
    tuple val(meta), path("*.vcf.gz")     , emit: vcf
    tuple val(meta), path("*.tbi")        , emit: tbi
    tuple val(meta), path("*.stats")      , emit: stats
    tuple val(meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def panels_command = ''
    def normals_command = ''

    def inputs_command = '-I ' + input.join( ' -I ')
    def interval = intervals ? "-L ${intervals}" : ""

    if(run_pon) {
        panels_command = ''
        normals_command = ''

    } else if(run_single) {
        panels_command = " --germline-resource $germline_resource --panel-of-normals $panel_of_normals"
        normals_command = ''

    } else if(run_mito){
        panels_command = "-L ${intervals} --mitochondria-mode"
        normals_command = ''

    } else {
        panels_command = " --germline-resource $germline_resource --panel-of-normals $panel_of_normals --f1r2-tar-gz ${prefix}.f1r2.tar.gz"
        normals_command = '-normal ' + which_norm.join( ' -normal ')
    }

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" Mutect2 \\
        -R ${fasta} \\
        ${inputs_command} \\
        ${normals_command} \\
        ${panels_command} \\
        ${interval} \\
        -O ${prefix}.vcf.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
