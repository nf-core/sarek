#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { MARKDUPLICATES } from '../../../../subworkflows/nf-core/markduplicates' addParams(
    markduplicates_options:      modules['markduplicates'],
    markduplicatesspark_options: modules['markduplicatesspark']
)

workflow test_markduplicates {
    input = [[id: 'test'],
            [file(params.test_data['nf-core']['test_paired_end_sorted_bam'], checkIfExists: true)],
            [file(params.test_data['nf-core']['test_paired_end_sorted_bai'], checkIfExists: true)]]

    MARKDUPLICATES ( input, false, true )
}
