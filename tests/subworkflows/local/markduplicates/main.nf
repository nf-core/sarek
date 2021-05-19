#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MARKDUPLICATES } from '../../../../subworkflow/local/markduplicates' addParams(
    markduplicates_options:          modules['markduplicates']
)

workflow test_markduplicates {
    input = [[ id: 'test' ],
             [ file(params.test_data['nf-core']['test_paired_end_sorted_bam'], checkIfExists: true)]]

    step = 'preparerecalibration'

    MARKDUPLICATES ( input, step )
}
