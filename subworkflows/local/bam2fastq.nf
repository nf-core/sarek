//
// BAM/CRAM to FASTQ conversion, paired end only
//


//include { FASTQC                 } from '../../modules/nf-core/modules/fastqc/main'
include { SAMTOOLS_INDEX                                    } from '../../modules/nf-core/modules/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MAP_MAP            } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_UNMAP_UNMAP        } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_UNMAP_MAP          } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MAP_UNMAP          } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_UNMAPPED         } from '../../modules/nf-core/modules/samtools/merge/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_UNMAPPED         } from '../../modules/local/samtools/fastq/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED           } from '../../modules/local/samtools/fastq/main'
include { CAT_FASTQ as CAT_FASTQ_R1                         } from '../../modules/nf-core/modules/cat/fastq/main'
include { CAT_FASTQ as CAT_FASTQ_R2                         } from '../../modules/nf-core/modules/cat/fastq/main'

workflow ALIGNMENT_TO_FASTQ {
    take:
    input // channel: [meta, alignment (BAM or CRAM), index (optional)]
    fasta // optional: reference file if CRAM format and reference not in header


    main:

    //Index File if not PROVIDED -> this also requires updates to samtools view possibly URGH

    //QC input BAM? -> needs another FASTQC module implementation

    //MAP - MAP
    SAMTOOLS_VIEW_MAP_MAP(input, fasta)

    // UNMAP - UNMAP
    SAMTOOLS_VIEW_UNMAP_UNMAP(input, fasta)

    // UNMAP - MAP
    SAMTOOLS_VIEW_UNMAP_MAP(input, fasta)

    //MAP - UNMAP
    SAMTOOLS_VIEW_MAP_UNMAP(input, fasta)

    // Merge UNMAP
    SAMTOOLS_VIEW_UNMAP_UNMAP.bam.join(SAMTOOLS_VIEW_UNMAP_MAP.bam, remainder: true)
                .join(SAMTOOLS_VIEW_MAP_UNMAP.bam, remainder: true)
                .set{ all_unmapped_bam }

    SAMTOOLS_MERGE(all_unmapped_bam, fasta)

    // Collate & convert unmapped
    SAMTOOLS_FASTQ_UNMAPPED(all_unmapped_bam)

    // Collate & convert mapped
    SAMTOOLS_FASTQ_MAPPED(SAMTOOLS_VIEW_MAP_MAP.out.bam)

    // join Mapped & unmapped fastq
    SAMTOOLS_FASTQ_UNMAPPED.out.reads.map{ meta, read1, read2, other, singleton ->
                                           [meta, read1]
                                        }.set(unmapped_r1)
    SAMTOOLS_FASTQ_UNMAPPED.out.reads.map{ meta, read1, read2, other, singleton ->
                                           [meta, read2]
                                        }.set(unmapped_r2)

    SAMTOOLS_FASTQ_MAPPED.out.reads.map{ meta, read1, read2, other, singleton ->
                                           [meta, read1]
                                        }.set(mapped_r1)
    SAMTOOLS_FASTQ_MAPPED.out.reads.map{ meta, read1, read2, other, singleton ->
                                           [meta, read2]
                                        }.set(mapped_r2)

    //TODO test if it is not actually possible to use the module with paired reads and concatenate pairs
    CAT_FASTQ_R1()
    CAT_FASTQ_R2()

    emit:
    reads

}
