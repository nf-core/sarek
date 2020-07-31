process QUALIMAP_BAMQC {
    // label 'memory_max'
    // label 'cpus_16'

    // tag "${meta.id}"

    // publishDir "${params.outdir}/Reports/${meta.id}/bamQC", mode: params.publish_dir_mode

    input:
         tuple val(meta), path(bam)
         path(targetBED) 

     output:
         path("${bam.baseName}") 

    // //when: !('bamqc' in skip_qc)

     script:
    // use_bed = params.target_bed ? "-gff ${targetBED}" : ''
     """
   # // qualimap --java-mem-size=${task.memory.toGiga()}G \
    #//     bamqc \
    #//     -bam ${bam} \
    #//     --paint-chromosome-limits \
    #//     --genome-gc-distr HUMAN \
    #//     $use_bed \
    #//     -nt ${task.cpus} \
    #//     -skip-duplicated \
    #//     --skip-dup-mode 0 \
    #//     -outdir ${bam.baseName} \
    #//     -outformat HTML
     """
}
