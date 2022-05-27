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

    // Only if we want to split fastq files
    SEQKIT_SPLIT2(reads_input)

    // Remapping the channel
    reads = SEQKIT_SPLIT2.out.reads.map{ key, reads ->
        //TODO maybe this can be replaced by a regex to include part_001 etc.

        //sorts list of split fq files by :
        //[R1.part_001, R2.part_001, R1.part_002, R2.part_002,R1.part_003, R2.part_003,...]
        //TODO: determine whether it is possible to have an uneven number of parts, so remainder: true woud need to be used, I guess this could be possible for unfiltered reads, reads that don't have pairs etc.
        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[ a.getName().tokenize('.').size() - 3] <=> b.getName().tokenize('.')[ b.getName().tokenize('.').size() - 3]}.collate(2)

        [[patient: key.patient, sample:key.sample, gender:key.gender, status:key.status, id:key.id, numLanes:key.numLanes, read_group:key.read_group, data_type:key.data_type, size:read_files.size()],
         read_files]
    }.transpose()

    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    emit:
        reads    = reads
        versions = ch_versions
}

