//
// BAM/CRAM to FASTQ conversion, paired end only
//

include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_MAP_MAP     } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_UNMAP_UNMAP } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_UNMAP_MAP   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_VIEW         as SAMTOOLS_VIEW_MAP_UNMAP   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_MERGE        as SAMTOOLS_MERGE_UNMAP      } from '../../../modules/nf-core/samtools/merge'
include { SAMTOOLS_COLLATEFASTQ as COLLATE_FASTQ_UNMAP       } from '../../../modules/nf-core/samtools/collatefastq'
include { SAMTOOLS_COLLATEFASTQ as COLLATE_FASTQ_MAP         } from '../../../modules/nf-core/samtools/collatefastq'
include { CAT_FASTQ                                          } from '../../../modules/nf-core/cat/fastq'

workflow BAM_CONVERT_SAMTOOLS {
    take:
    input       // channel: [meta, alignment (BAM or CRAM), index (optional)]
    fasta       // channel: [meta, fasta] — reference, staged in to decode CRAM input
    fasta_fai   // channel: [meta, fai]
    interleaved // value: true/false

    main:

    // Combined [ meta, fasta, fai ] reference tuple for the updated samtools modules
    fasta_and_fai = fasta.combine(fasta_fai).map { meta, fasta_, _fai_meta, fai -> [ meta, fasta_, fai ] }

    // MAP - MAP
    SAMTOOLS_VIEW_MAP_MAP(input, fasta_and_fai, [[], []], [[], []], [])

    // UNMAP - UNMAP
    SAMTOOLS_VIEW_UNMAP_UNMAP(input, fasta_and_fai, [[], []], [[], []], [])

    // UNMAP - MAP
    SAMTOOLS_VIEW_UNMAP_MAP(input, fasta_and_fai, [[], []], [[], []], [])

    // MAP - UNMAP
    SAMTOOLS_VIEW_MAP_UNMAP(input, fasta_and_fai, [[], []], [[], []], [])

    // Merge UNMAP
    all_unmapped_bam = SAMTOOLS_VIEW_UNMAP_UNMAP.out.bam
        .join(SAMTOOLS_VIEW_UNMAP_MAP.out.bam, failOnDuplicate: true, remainder: true)
        .join(SAMTOOLS_VIEW_MAP_UNMAP.out.bam, failOnDuplicate: true, remainder: true)
        .map{ meta, unmap_unmap, unmap_map, map_unmap -> [ meta, [ unmap_unmap, unmap_map, map_unmap ] ] }

    SAMTOOLS_MERGE_UNMAP(all_unmapped_bam.map { meta, bams -> [ meta, bams, [] ] }, fasta_and_fai.map { meta, fasta_, fai -> [ meta, fasta_, fai, [] ] })

    // Collate & convert unmapped
    COLLATE_FASTQ_UNMAP(SAMTOOLS_MERGE_UNMAP.out.bam, fasta_and_fai, interleaved)

    // Collate & convert mapped
    COLLATE_FASTQ_MAP(SAMTOOLS_VIEW_MAP_MAP.out.bam, fasta_and_fai, interleaved)

    // join Mapped & unmapped fastq

    reads_to_concat = COLLATE_FASTQ_MAP.out.fastq
        .join(COLLATE_FASTQ_UNMAP.out.fastq, failOnDuplicate: true, failOnMismatch: true)
        .map{ meta, mapped_reads, unmapped_reads -> [ meta, [ mapped_reads[0], mapped_reads[1], unmapped_reads[0], unmapped_reads[1] ] ] }

    // Concatenate Mapped_R1 with Unmapped_R1 and Mapped_R2 with Unmapped_R2
    CAT_FASTQ(reads_to_concat)
    reads = CAT_FASTQ.out.reads

    emit:
    reads

}
