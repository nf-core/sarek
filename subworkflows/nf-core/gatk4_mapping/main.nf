//
// MAPPING
//

include { BWAMEM2_MEM                     } from '../../../modules/nf-core/modules/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM          } from '../../../modules/nf-core/modules/bwa/mem/main'
include { SAMTOOLS_INDEX as INDEX_MAPPING } from '../../../modules/local/samtools/index/main'
include { SAMTOOLS_MERGE as MERGE_MAPPING } from '../../../modules/nf-core/modules/samtools/merge/main'
include { SEQKIT_SPLIT2                   } from '../../../modules/nf-core/modules/seqkit/split2/main'

workflow GATK4_MAPPING {
    take:
        reads_input // channel: [mandatory] meta, reads_input
        bwa         // channel: [mandatory] bwa
        fasta       // channel: [mandatory] fasta
        fasta_fai   // channel: [mandatory] fasta_fai

    main:

    ch_versions = Channel.empty()

    reads_input = reads_input.map{ meta, reads ->
        meta.size = 1
        [meta, reads]
    }

    SEQKIT_SPLIT2(reads_input)

    reads_input_split = SEQKIT_SPLIT2.out.reads.map{ key, reads ->
        //TODO maybe this can be replaced by a regex to include part_001 etc.

        //sorts list of split fq files by :
        //[R1.part_001, R2.part_001, R1.part_002, R2.part_002,R1.part_003, R2.part_003,...]
        //TODO: determine whether it is possible to have an uneven number of parts, so remainder: true woud need to be used, I guess this could be possible for unfiltered reads, reads that don't have pairs etc.
        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[ a.getName().tokenize('.').size() - 3] <=> b.getName().tokenize('.')[ b.getName().tokenize('.').size() - 3]}.collate(2)
        key.size = read_files.size()
        [key, read_files]
    }.transpose()

    reads_input_all = reads_input.mix(reads_input_split)

    BWAMEM1_MEM(reads_input_all, bwa, true)
    BWAMEM2_MEM(reads_input_all, bwa, true)

    ch_versions = ch_versions.mix(BWAMEM1_MEM.out.versions.first())
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    bam_mapped = BWAMEM1_MEM.out.bam.mix(BWAMEM2_MEM.out.bam).map{ meta, bam ->
        new_meta = meta.clone()
        new_meta.remove('read_group')
        new_meta.remove('size')
        new_meta.id = meta.sample

        // groupKey is to makes sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all pieces are mapped
        def groupKey = groupKey(meta, meta.numLanes * meta.size)

        //Returns the values we need
        tuple(groupKey, new_meta, bam)
    }.groupTuple(by:[0,1]).map{
        groupKey, new_meta, bam ->
        [new_meta, bam]
    }

    // GATK markduplicates can handle multiple BAMS as input
    // So no merging/indexing at this step
    // Except if and only if skipping markduplicates
    // Or saving mapped BAMs

    bam_mapped.branch{
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set{bam_to_merge}

    MERGE_MAPPING(bam_to_merge.multiple, [])

    INDEX_MAPPING(bam_to_merge.single.mix(MERGE_MAPPING.out.bam))

    emit:
        bam         = bam_mapped
        bam_indexed = INDEX_MAPPING.out.bam_bai
        versions    = ch_versions
}
