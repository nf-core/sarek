//
// SPLIT_FASTQ
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

workflow SPLIT_FASTQ {
    take:
        reads_input // channel: [mandatory] meta, reads_input

    main:

    // Remapping the channel
    reads = reads_input.map{ key, reads ->

        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
        [[patient: key.patient, sample:key.sample, gender:key.gender, status:key.status, id:key.id, numLanes:key.numLanes, read_group:key.read_group, data_type:key.data_type, size:read_files.size()],
        read_files]
    }.transpose()

    emit:
        reads    = reads
}

