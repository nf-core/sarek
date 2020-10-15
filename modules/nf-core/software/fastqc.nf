include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

environment = params.enable_conda ? "bioconda::fastqc=0.11.9" : null
container = "quay.io/biocontainers/fastqc:0.11.9--0"
if (workflow.containerEngine == 'singularity' && !params.pull_docker_container) container = "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"

process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    label 'cpus_2'

    publishDir "${params.outdir}/${options.publish_dir}/${meta.sample}/${meta.id}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (options.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    conda environment
    container container

    input:
        tuple val(meta), path(reads)

    output:
        path "*.html",        emit: html
        path "*.version.txt", emit: version
        path "*.zip",         emit: zip

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    prefix = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz

    fastqc ${options.args} --threads ${task.cpus} ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz

    fastqc --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > fastqc.version.txt
    """
}