#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { VEP_ANNOTATE } from '../../../../subworkflows/nf-core/vep_annotate' addParams(
    bgziptabix_vep_options: modules['bgziptabix_vep'],
    use_cache:              false,
    vep_options:            modules['vep'],
    vep_tag:                "${modules['vep'].tag_base}.WBcel235"
)

workflow test_vep_annotate {
    input = [[id: 'test'],
            [file(params.test_data['sarscov2']['illumina']['test_vcf'], checkIfExists: true)]]

    VEP_ANNOTATE (
        input,
        "WBcel235",
        "caenorhabditis_elegans",
        "104",
        [])
}
