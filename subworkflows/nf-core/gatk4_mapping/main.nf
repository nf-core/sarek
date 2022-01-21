//
// MAPPING
//

include { BWAMEM2_MEM                     } from '../../../modules/nf-core/modules/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM          } from '../../../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_INDEX as INDEX_MAPPING } from '../../../modules/local/samtools/index/main'
include { SAMTOOLS_MERGE                  } from '../../../modules/nf-core/modules/samtools/merge/main'
include { SEQKIT_SPLIT2                   } from '../../../modules/nf-core/modules/seqkit/split2/main'

workflow GATK4_MAPPING {
    take:
        aligner             // string:  [mandatory] "bwa-mem" or "bwa-mem2"
        bwa                 // channel: [mandatory] bwa
        fasta               // channel: [mandatory] fasta
        fasta_fai           // channel: [mandatory] fasta_fai
        reads_input         // channel: [mandatory] meta, reads_input
        skip_markduplicates // boolean: true/false
        save_bam_mapped     // boolean: true/false

    main:

    bam_indexed = Channel.empty()

    if (params.split_fastq > 1) {
        reads_input_split = SEQKIT_SPLIT2(reads_input).reads.map{ key, reads ->
            //TODO maybe this can be replaced by a regex to include part_001 etc.

            //sorts list of split fq files by :
            //[R1.part_001, R2.part_001, R1.part_002, R2.part_002,R1.part_003, R2.part_003,...]
            //TODO: determine whether it is possible to have an uneven number of parts, so remainder: true woud need to be used, I guess this could be possible for unfiltered reads, reads that don't have pairs etc.
            read_files = reads.sort{ a,b -> a.getName().tokenize('.')[ a.getName().tokenize('.').size() - 3] <=> b.getName().tokenize('.')[ b.getName().tokenize('.').size() - 3]}.collate(2)
            key.size = read_files.size()
            [key, read_files]
        }.transpose()
    } else reads_input_split = reads_input



    bam_bwamem1      = Channel.empty()
    bam_bwamem2      = Channel.empty()
    bam_from_aligner = Channel.empty()
    tool_versions    = Channel.empty()

    if (aligner == "bwa-mem") {
        BWAMEM1_MEM(reads_input_split, bwa, true)

        bam_bwamem1 = BWAMEM1_MEM.out.bam

        bwamem1_version = BWAMEM1_MEM.out.versions.first()

        tool_versions = tool_versions.mix(bwamem1_version)
    } else {
        BWAMEM2_MEM(reads_input_split, bwa, true)

        bam_bwamem2 = BWAMEM2_MEM.out.bam

        bwamem2_version = BWAMEM2_MEM.out.versions.first()

        tool_versions = tool_versions.mix(bwamem2_version)
    }

    bam_from_aligner = bam_from_aligner.mix(bam_bwamem1)
    bam_from_aligner = bam_from_aligner.mix(bam_bwamem2)

    bam_from_aligner.map{ meta, bam ->
        new_meta = meta.clone()
        new_meta.remove('read_group')
        new_meta.id = meta.sample

        // groupKey is to makes sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all pieces are mapped
        def groupKey = groupKey(meta, meta.numLanes * meta.size)
        //Returns the values we need
        tuple(groupKey, new_meta, bam)
    }.groupTuple(by:[0,1]).map{ groupKey, new_meta, bam ->
        println new_meta.getClass()
        println bam.getClass()
        [new_meta, bam]}.set{bam_mapped}

    bam_mapped.view()
    // GATK markduplicates can handle multiple BAMS as input
    // So no merging/indexing at this step
    // Except if and only if skipping markduplicates
    // Or saving mapped BAMs

    if (save_bam_mapped || skip_markduplicates) {
        bam_mapped.branch{
            single:   it[1].size() == 1
            multiple: it[1].size() > 1
        }.set{bam_to_merge}

        SAMTOOLS_MERGE(bam_to_merge.multiple, [])
        bam_merged = bam_to_merge.single.mix(SAMTOOLS_MERGE.out.bam)

        INDEX_MAPPING(bam_merged)
        bam_indexed = INDEX_MAPPING.out.bam_bai
    }

    emit:
        bam         = bam_mapped
        bam_indexed = bam_indexed
        versions    = tool_versions
}
