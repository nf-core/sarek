#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { ANNOTATE } from '../../../../subworkflows/local/annotate' addParams(
    annotation_cache:               false,
    bgziptabix_merge_vep_options:   modules['bgziptabix_merge_vep'],
    bgziptabix_snpeff_options:      modules['bgziptabix_snpeff'],
    bgziptabix_vep_options:         modules['bgziptabix_vep'],
    merge_vep_options:              modules['merge_vep'],
    snpeff_options:                 modules['snpeff'],
    snpeff_tag:                     "${modules['snpeff'].tag_base}.WBcel235",
    vep_options:                    modules['vep'],
    vep_tag:                        "${modules['vep'].tag_base}.WBcel235"
)

workflow test_annotate {
    input = [[id: 'test'],
            [file(params.test_data['nf-core']['test_vcf'], checkIfExists: true)]]

    ANNOTATE(
        input,
        ["snpeff","vep","merge"],
        "WBcel235.99",
        [],
        "WBcel235",
        "caenorhabditis_elegans",
        "104",
        [])
}
