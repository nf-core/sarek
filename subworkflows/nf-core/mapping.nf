/*
========================================================================================
    MAPPING
========================================================================================
*/

params.seqkit_split2_options     = [:]
params.bwamem1_mem_options       = [:]
params.bwamem1_mem_tumor_options = [:]
params.bwamem2_mem_options       = [:]
params.bwamem2_mem_tumor_options = [:]

include { SEQKIT_SPLIT2 }                from '../../modules/nf-core/modules/seqkit/split2/main.nf' addParams(options: params.seqkit_split2_options)
include { BWA_MEM as BWAMEM1_MEM }       from '../../modules/local/bwa/mem/main'                    addParams(options: params.bwamem1_mem_options)
include { BWA_MEM as BWAMEM1_MEM_T }     from '../../modules/local/bwa/mem/main'                    addParams(options: params.bwamem1_mem_tumor_options)
include { BWAMEM2_MEM }                  from '../../modules/local/bwamem2/mem/main'                addParams(options: params.bwamem2_mem_options)
include { BWAMEM2_MEM as BWAMEM2_MEM_T } from '../../modules/local/bwamem2/mem/main'                addParams(options: params.bwamem2_mem_tumor_options)

workflow MAPPING {
    take:
        aligner         // string:  [mandatory] "bwa-mem" or "bwa-mem2"
        bwa             // channel: [mandatory] bwa
        fai             // channel: [mandatory] fai
        fasta           // channel: [mandatory] fasta
        reads_input     // channel: [mandatory] meta, reads_input

    main:

    bam_mapped_index = Channel.empty()
    bam_reports      = Channel.empty()


    if(params.split_fastq > 1){
        reads_input_split = SEQKIT_SPLIT2(reads_input).reads.map{
                key, reads ->
                    //TODO maybe this can be replaced by a regex to include part_001 etc.

                    //sorts list of split fq files by :
                    //[R1.part_001, R2.part_001, R1.part_002, R2.part_002,R1.part_003, R2.part_003,...]
                    //TODO: determine whether it is possible to have an uneven number of parts, so remainder: true woud need to be used, I guess this could be possible for unfiltered reads, reads that don't have pairs etc.
                    return [key, reads.sort{ a,b -> a.getName().tokenize('.')[ a.getName().tokenize('.').size() - 3] <=> b.getName().tokenize('.')[ b.getName().tokenize('.').size() - 3]}
                                            .collate(2)]
            }.transpose()
    }else{
        reads_input_split = reads_input
    }

    // If meta.status is 1, then sample is tumor
    // else, (even is no meta.status exist) sample is normal
    reads_input_split.branch{
            tumor:  it[0].status == 1
            normal: true
        }.set{ reads_input_status }

    bam_bwamem1 = Channel.empty()
    bam_bwamem2 = Channel.empty()

    if (aligner == "bwa-mem") {
        BWAMEM1_MEM(reads_input_status.normal, bwa)
        bam_bwamem1_n = BWAMEM1_MEM.out.bam

        BWAMEM1_MEM_T(reads_input_status.tumor, bwa)
        bam_bwamem1_t = BWAMEM1_MEM_T.out.bam

        bam_bwamem1 = bam_bwamem1_n.mix(bam_bwamem1_t)
    } else {
        BWAMEM2_MEM(reads_input_status.normal, bwa)
        bam_bwamem2_n = BWAMEM2_MEM.out.bam

        BWAMEM2_MEM_T(reads_input_status.tumor, bwa)
        bam_bwamem2_t = BWAMEM2_MEM_T.out.bam

        bam_bwamem2 = bam_bwamem2_n.mix(bam_bwamem2_t)
    }
    bam_bwa = bam_bwamem1.mix(bam_bwamem2)

    bam_bwa.map{ meta, bam ->
        meta.remove('read_group')
        meta.id = meta.sample
        // groupKey is to makes sure that the correct group can advance as soon as it is complete
        // and not stall the workflow until all pieces are mapped
        def groupKey = groupKey(meta, meta.numLanes * params.split_fastq)
        tuple(groupKey, bam)
        [meta, bam]
    }.groupTuple()
    .set{bam_mapped}

    // MarkDuplicates can handles multiple BAMS as input, so no merging/indexing at this step

    emit:
        bam = bam_mapped
}
