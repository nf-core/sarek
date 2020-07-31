process SAMTOOLS_INDEX {
   label 'cpus_8'

    tag "${meta.id}"

//     publishDir params.outdir, mode: params.publish_dir_mode,
//         saveAs: {
//             if (save_bam_mapped) "Preprocessing/${idSample}/Mapped/${it}"
//             else null
//         }

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path(bam), path("*.bai")
        //set idPatient, idSample into tsv_bam_indexed

    script:
    """
    samtools index $bam
    """
    //     samtools index ${idSample}.recal.bam TODO: is the naming here relevant?
}