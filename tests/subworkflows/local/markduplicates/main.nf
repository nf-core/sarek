#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { MARKDUPLICATES } from '../../../../subworkflow/local/markduplicates' addParams(
    markduplicates_options: modules['markduplicates']
)

workflow test_markduplicates {
    input = [[id: 'test'],
            [file(params.test_data['nf-core']['test_paired_end_sorted_bam'], checkIfExists: true)]]

    MARKDUPLICATES ( input, true )
}
