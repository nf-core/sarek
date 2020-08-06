process SAMTOOLS_INDEX {
   label 'cpus_8'

    tag "${meta.id}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (save_bam_mapped) "Preprocessing/${meta.sample}/Mapped/${it}"
            else null
        }

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