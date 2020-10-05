process QUALIMAP_BAMQC {
    label 'memory_max'
    label 'cpus_16'

    tag "${meta.id}"

    publishDir "${params.outdir}/Reports/${meta.id}/bamQC", mode: params.publish_dir_mode

    container "quay.io/biocontainers/qualimap:2.2.2d--1"

    conda (params.conda ? "bioconda::qualimap=2.2.2d" : null)

    input:
         tuple val(meta), path(bam)
         path(target_bed)

     output:
         path("${bam.baseName}") 

     script:
     use_bed = params.target_bed ? "-gff ${target_bed}" : ''
     """
     qualimap --java-mem-size=${task.memory.toGiga()}G \
        bamqc \
        -bam ${bam} \
        --paint-chromosome-limits \
        --genome-gc-distr HUMAN \
        ${use_bed} \
        -nt ${task.cpus} \
        -skip-duplicated \
        --skip-dup-mode 0 \
        -outdir ${bam.baseName} \
        -outformat HTML
     """
}
