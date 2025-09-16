include { MINIMAP2_ALIGN } from '../../../modules/local/minimap2/align/main'

workflow BAM_ALIGN_MINIMAP2 {
    take:
    reads     // channel: [mandatory] meta, reads
    index     // channel: [mandatory] meta, index
    sort      // boolean: [mandatory] true -> sort, false -> don't sort
    fasta     // channel: fasta
    fasta_fai // channel: fai

    main:
    versions = Channel.empty()

    // Align with minimap2
    // The index channel already has the meta information from PREPARE_GENOME
    MINIMAP2_ALIGN(
        reads,
        index,
        sort,
        fasta,
        fasta_fai
    )
    
    versions = versions.mix(MINIMAP2_ALIGN.out.versions)

    emit:
    bam      = MINIMAP2_ALIGN.out.bam
    bai      = MINIMAP2_ALIGN.out.bai
    versions = versions
}