process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    label 'cpus_2'

    publishDir "${params.outdir}/Reports/${meta.sample}/FastQC/${meta.id}",
        mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith('.version.txt')) null
                    else filename }

    container "quay.io/biocontainers/fastqc:0.11.9--0"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*.html"), emit: html
        tuple val(meta), path("*.zip"),  emit: zip
        path "*.version.txt",            emit: version

    script:
    prefix = "${meta.id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    fastqc --threads ${task.cpus} ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz

    fastqc --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > fastqc.version.txt
    """
}
