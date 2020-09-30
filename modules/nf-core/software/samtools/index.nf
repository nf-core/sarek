process SAMTOOLS_INDEX {
   label 'cpus_8'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (params.save_bam_mapped) "Preprocessing/${meta.sample}/Mapped/${it}"
            else null
        }

    container "quay.io/biocontainers/samtools:1.10--h2e538c0_3"

    conda (params.conda ? "bioconda::samtools=1.10" : null)

    input:
        tuple val(meta), path(bam)
        val options

    output:
        tuple val(meta), path("${prefix}.bam"), path("*.bai")

    script:
    prefix = options.suffix ? "${meta.id}.${options.suffix}" : "${meta.id}"
    """
    [ ! -f  ${prefix}.bam ] && ln -s ${bam} ${prefix}.bam

    samtools index ${prefix}.bam
    """
}