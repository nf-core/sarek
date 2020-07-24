process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    label 'cpus_2'

    publishDir "${params.outdir}/${options.publish_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (options.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/fastqc:0.11.9--0"

    input:
        tuple val(meta), path(reads)
        val options

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