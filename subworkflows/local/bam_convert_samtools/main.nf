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
    input       // channel: [meta, alignment (BAM), index (optional)]
    interleaved // value: true/false

    main:

    // BAM -> FASTQ needs no reference, but the updated samtools modules still
    // require the reference tuple, so pass empty placeholders.
    reference = [[id: 'fasta'], [], []] // [ meta, fasta, fai ]

    // MAP - MAP
    SAMTOOLS_VIEW_MAP_MAP(input, reference, [[], []], [[], []], [])

    // UNMAP - UNMAP
    SAMTOOLS_VIEW_UNMAP_UNMAP(input, reference, [[], []], [[], []], [])

    // UNMAP - MAP
    SAMTOOLS_VIEW_UNMAP_MAP(input, reference, [[], []], [[], []], [])

    // MAP - UNMAP
    SAMTOOLS_VIEW_MAP_UNMAP(input, reference, [[], []], [[], []], [])

    // Merge UNMAP
    all_unmapped_bam = SAMTOOLS_VIEW_UNMAP_UNMAP.out.bam
        .join(SAMTOOLS_VIEW_UNMAP_MAP.out.bam, failOnDuplicate: true, remainder: true)
        .join(SAMTOOLS_VIEW_MAP_UNMAP.out.bam, failOnDuplicate: true, remainder: true)
        .map{ meta, unmap_unmap, unmap_map, map_unmap -> [ meta, [ unmap_unmap, unmap_map, map_unmap ] ] }

    SAMTOOLS_MERGE_UNMAP(all_unmapped_bam.map { meta, bams -> [ meta, bams, [] ] }, reference + [[]])

    // Collate & convert unmapped
    COLLATE_FASTQ_UNMAP(SAMTOOLS_MERGE_UNMAP.out.bam, reference, interleaved)

    // Collate & convert mapped
    COLLATE_FASTQ_MAP(SAMTOOLS_VIEW_MAP_MAP.out.bam, reference, interleaved)

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
