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
    //SEQKIT_SPLIT2(reads_input)

    // Remapping the channel
    reads = reads_input.map{ key, reads ->
        //sorts list of split fq files by :
        //[0001.Read_1.trim.fastq.gz,0001.Read_2.trim.fastq.gz,0002.Read_1.trim.fastq.gz,0002.Read_2.trim.fastq.gz]
        // 1. Get element until first . -> 0001
        // 2. sort by it: 0001_Read1, 0001_Read2, 0002_Read2, 0002_Read1 etc
        // 3. Collate two elements together: [0001_Read1, 0001_Read2], [0002_Read2, 0002_Read1]
        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0]}.collate(2)
        [[patient: key.patient, sample:key.sample, gender:key.gender, status:key.status, id:key.id, numLanes:key.numLanes, read_group:key.read_group, data_type:key.data_type, size:read_files.size()],
        read_files]
    }.transpose()
    //ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    emit:
        reads    = reads
        versions = ch_versions
}

