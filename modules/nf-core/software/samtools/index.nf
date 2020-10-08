include { initOptions; saveFiles; getSoftwareName } from './../functions'

process SAMTOOLS_INDEX {
   label 'cpus_8'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (options.publish_results == "none") null
                    else if (filename.endsWith('.version.txt')) null
                    else "${options.publish_dir_up}/${meta.sample}/${options.publish_dir_down}/${filename}" }

    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
        tuple val(meta), path(bam)
        val options

    output:
        tuple val(meta), path("${name}.bam"), path("*.bai")

    script:
    name = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    [ ! -f  ${name}.bam ] && ln -s ${bam} ${name}.bam

    samtools index ${name}.bam
    """
}