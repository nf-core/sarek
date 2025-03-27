//
// Uncompress and prepare reference genome files
//
include { UNTAR as UNTAR_BBSPLIT_INDEX      } from '../../../modules/nf-core/untar'

include { BBMAP_BBSPLIT                     } from '../../../modules/nf-core/bbmap/bbsplit'

workflow PREPARE_BBSPLIT {
    take:
    bbsplit_fasta_list   //      file: /path/to/bbsplit_fasta_list.txt
    bbsplit_index        // directory: /path/to/rsem/index/
    prepare_tool_indices //      list: tools to prepare indices for

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress BBSplit index or generate from scratch if required
    //
    ch_bbsplit_index = Channel.empty()
    if ('bbsplit' in prepare_tool_indices) {
        if (bbsplit_index) {
            if (bbsplit_index.endsWith('.tar.gz')) {
                ch_bbsplit_index = UNTAR_BBSPLIT_INDEX ( [ [:], bbsplit_index ] ).untar.map { it[1] }
                ch_versions      = ch_versions.mix(UNTAR_BBSPLIT_INDEX.out.versions)
            } else {
                ch_bbsplit_index = Channel.value(file(bbsplit_index))
            }
        } else {
            Channel
                .from(file(bbsplit_fasta_list))
                .splitCsv() // Read in 2 column csv file: short_name,path_to_fasta
                .flatMap { id, fasta -> [ [ 'id', id ], [ 'fasta', file(fasta, checkIfExists: true) ] ] } // Flatten entries to be able to groupTuple by a common key
                .groupTuple()
                .map { it -> it[1] } // Get rid of keys and keep grouped values
                .collect { [ it ] } // Collect entries as a list to pass as "tuple val(short_names), path(path_to_fasta)" to module
                .set { ch_bbsplit_fasta_list }

            ch_bbsplit_index = BBMAP_BBSPLIT ( [ [:], [] ], [], ch_fasta, ch_bbsplit_fasta_list, true ).index
            ch_versions      = ch_versions.mix(BBMAP_BBSPLIT.out.versions)
        }
    }

    emit:
    bbsplit_index    = ch_bbsplit_index          // channel: path(bbsplit/index/)
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
