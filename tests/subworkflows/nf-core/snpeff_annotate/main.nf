#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { SNPEFF_ANNOTATE } from '../../../../subworkflows/nf-core/snpeff_annotate' addParams(
    bgziptabix_snpeff_options: modules['bgziptabix_snpeff'],
    snpeff_options:            modules['snpeff'],
    snpeff_tag:                "${modules['snpeff'].tag_base}.WBcel235",
    use_cache:                 false
)

workflow test_snpeff_annotate {
    input = [[id: 'test'],
            [file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)]]

    SNPEFF_ANNOTATE (
        input,
        "WBcel235.99",
        [])
}
