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
include { CAT_FASTQ                                         } from '../../modules/nf-core/modules/cat/fastq/main'

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
    SAMTOOLS_VIEW_UNMAP_UNMAP.out.bam.join(SAMTOOLS_VIEW_UNMAP_MAP.out.bam, remainder: true)
                .join(SAMTOOLS_VIEW_MAP_UNMAP.out.bam, remainder: true)
                .map{ meta, unmap_unmap, unmap_map, map_unmap ->
                    [meta, [unmap_unmap, unmap_map, map_unmap]]
                }.set{ all_unmapped_bam }

    SAMTOOLS_MERGE_UNMAPPED(all_unmapped_bam, fasta)

    // Collate & convert unmapped
    SAMTOOLS_FASTQ_UNMAPPED(SAMTOOLS_MERGE_UNMAPPED.out.bam)

    // Collate & convert mapped
    SAMTOOLS_FASTQ_MAPPED(SAMTOOLS_VIEW_MAP_MAP.out.bam)

    // join Mapped & unmapped fastq
    SAMTOOLS_FASTQ_UNMAPPED.out.reads.map{ meta, reads ->
                fq_1 = reads.findAll{ it.toString().endsWith("_1.fq.gz") }.get(0)
                fq_2 = reads.findAll{ it.toString().endsWith("_2.fq.gz") }.get(0)
                [meta, [ fq_1, fq_2]]
                }.set{unmapped_reads}

    SAMTOOLS_FASTQ_MAPPED.out.reads.map{ meta, reads ->
                fq_1 = reads.findAll{ it.toString().endsWith("_1.fq.gz") }.get(0)
                fq_2 = reads.findAll{ it.toString().endsWith("_2.fq.gz") }.get(0)
                [meta, [ fq_1, fq_2]]
                }.set{mapped_reads}

    mapped_reads.join(unmapped_reads).map{ meta, mapped_reads, unmapped_reads ->
                                [meta, [mapped_reads[0], mapped_reads[1], unmapped_reads[0], unmapped_reads[1]]]
                                }.set{ reads_to_concat }

    // Concatenate Mapped_R1 with Unmapped_R1 and Mapped_R2 with Unmapped_R2
    CAT_FASTQ(reads_to_concat)

    emit:
    reads = CAT_FASTQ.out.reads

}
