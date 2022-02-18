//
// SPLIT_FASTQ
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SEQKIT_SPLIT2 } from '../../modules/nf-core/modules/seqkit/split2/main'

workflow SPLIT_FASTQ {
    take:
        reads_input // channel: [mandatory] meta, reads_input

    main:

    ch_versions = Channel.empty()

    reads_input = reads_input.map{ meta, reads ->
        meta.size = 1
        [meta, reads]
    }

    // Only if we want to split fastq files
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

    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    emit:
        reads       = reads_input_all
        versions    = ch_versions
}

