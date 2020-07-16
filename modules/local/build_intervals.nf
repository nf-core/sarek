process BUILD_INTERVALS {
    tag "${fai}"

    publishDir params.outdir, mode: params.publish_dir_mode,
    saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path fai

    output:
        path "${fai.baseName}.bed"

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fai} > ${fai.baseName}.bed
    """
}