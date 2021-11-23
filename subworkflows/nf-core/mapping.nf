//
// MAPPING
//

params.bwamem1_mem_options       = [:]
params.bwamem1_mem_tumor_options = [:]
params.bwamem2_mem_options       = [:]
params.bwamem2_mem_tumor_options = [:]
params.merge_bam_options         = [:]
params.samtools_index_options    = [:]
params.seqkit_split2_options     = [:]

include { BWAMEM2_MEM as BWAMEM2_MEM_T } from '../../modules/local/bwamem2/mem/main'                addParams(options: params.bwamem2_mem_tumor_options)
include { BWAMEM2_MEM }                  from '../../modules/local/bwamem2/mem/main'                addParams(options: params.bwamem2_mem_options)
include { BWA_MEM as BWAMEM1_MEM }       from '../../modules/local/bwa/mem/main'                    addParams(options: params.bwamem1_mem_options)
include { BWA_MEM as BWAMEM1_MEM_T }     from '../../modules/local/bwa/mem/main'                    addParams(options: params.bwamem1_mem_tumor_options)
include { SAMTOOLS_INDEX }               from '../../modules/local/samtools/index/main'             addParams(options: params.samtools_index_options)
include { SAMTOOLS_MERGE }               from '../../modules/nf-core/modules/samtools/merge/main'   addParams(options: params.merge_bam_options)
include { SEQKIT_SPLIT2 }                from '../../modules/nf-core/modules/seqkit/split2/main.nf' addParams(options: params.seqkit_split2_options)

workflow MAPPING {
    take:
        aligner             // string:  [mandatory] "bwa-mem" or "bwa-mem2"
        bwa                 // channel: [mandatory] bwa
        fai                 // channel: [mandatory] fai
        fasta               // channel: [mandatory] fasta
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
            [key, reads.sort{ a,b -> a.getName().tokenize('.')[ a.getName().tokenize('.').size() - 3] <=> b.getName().tokenize('.')[ b.getName().tokenize('.').size() - 3]}.collate(2)]
        }.transpose()
    } else reads_input_split = reads_input

    // If meta.status is 1, then sample is tumor
    // else, (even is no meta.status exist) sample is normal
    reads_input_split.branch{ meta, reads ->
        tumor:  meta.status == 1
        normal: true
    }.set{reads_input_status}

    bam_bwamem      = Channel.empty()
    bam_bwamem1     = Channel.empty()
    bam_bwamem2     = Channel.empty()
    tool_versions   = Channel.empty()

    if (aligner == "bwa-mem") {
        BWAMEM1_MEM(reads_input_status.normal, bwa)
        BWAMEM1_MEM_T(reads_input_status.tumor, bwa)

        bam_bwamem1_n = BWAMEM1_MEM.out.bam
        bam_bwamem1_t = BWAMEM1_MEM_T.out.bam
        bam_bwamem1   = bam_bwamem1_n.mix(bam_bwamem1_t)

        bwamem1_n_version = BWAMEM1_MEM.out.version
        bwamem1_t_version = BWAMEM1_MEM_T.out.version

        bwamem1_version = bwamem1_n_version.mix(bwamem1_t_version).first()

        tool_versions = tool_versions.mix(bwamem1_version)
    } else {
        BWAMEM2_MEM(reads_input_status.normal, bwa)
        BWAMEM2_MEM_T(reads_input_status.tumor, bwa)

        bam_bwamem2_n = BWAMEM2_MEM.out.bam
        bam_bwamem2_t = BWAMEM2_MEM_T.out.bam
        bam_bwamem2   = bam_bwamem2_n.mix(bam_bwamem2_t)

        bwamem2_n_version = BWAMEM2_MEM.out.version
        bwamem2_t_version = BWAMEM2_MEM_T.out.version

        bwamem2_version = bwamem2_n_version.mix(bwamem2_t_version).first()
        tool_versions = tool_versions.mix(bwamem2_version)
    }

    bam_bwamem = bam_bwamem.mix(bam_bwamem1)
    bam_bwamem = bam_bwamem.mix(bam_bwamem2)

    bam_bwamem.map{ meta, bam ->
        new_meta = meta.clone()
        new_meta.remove('read_group')
        new_meta.id = meta.sample

        // groupKey is to makes sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all pieces are mapped
        def groupKey = groupKey(meta, meta.numLanes * params.split_fastq)
        tuple(groupKey, bam)

        [new_meta, bam]
    }.groupTuple().set{bam_mapped}

    // GATK markduplicates can handle multiple BAMS as input
    // So no merging/indexing at this step
    // Except if and only if skipping markduplicates
    // Or saving mapped BAMs

    if (save_bam_mapped || skip_markduplicates) {
        bam_mapped.branch{
            single:   it[1].size() == 1
            multiple: it[1].size() > 1
        }.set{bam_to_merge}

        SAMTOOLS_MERGE(bam_to_merge.multiple)
        bam_merged = bam_to_merge.single.mix(SAMTOOLS_MERGE.out.bam)

        SAMTOOLS_INDEX(bam_merged)
        bam_indexed = SAMTOOLS_INDEX.out.bam_bai
    }

    emit:
        bam         = bam_mapped
        bam_indexed = bam_indexed
        versions    = tool_versions
}
