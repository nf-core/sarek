/*
 * FastQC
 */
process FASTQC {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    tuple val(name), val(single_end), path(reads)

    output:
    path "*.{zip,html}"

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    if (single_end) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        fastqc --quiet --threads $task.cpus ${name}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        fastqc --quiet --threads $task.cpus ${name}_1.fastq.gz ${name}_2.fastq.gz
        """
    }
}
