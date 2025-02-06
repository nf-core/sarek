// DISCLAIMER:
// This subworkflow is just to test the functions and the schema
// It should not be used in any pipeline

// This include statement can also be deleted
include { samplesheetToList } from 'plugin/nf-schema'

workflow UTILS_REFERENCES {
    take:
    yaml_reference
    param_file
    param_value
    attribute_file
    attribute_value
    basepath

    main:
    references = Channel.fromList(samplesheetToList(yaml_reference, "${projectDir}/subworkflows/nf-core/utils_references/schema_references.json"))

    // GIVING up writing a test for the functions, so writing a subworkflow to test it
    references_file = get_references_file(references, param_file, attribute_file, basepath)
    references_value = get_references_value(references, param_value, attribute_value)

    emit:
    references_file
    references_value
}
// You can delete everything before this line (including this line)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS TO EXTRACT REFERENCES FILES OR VALUES FROM THE REFERENCES YAML OR PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def get_references_file(references, param, attribute, basepath) {
    return references
        .map { meta, _readme ->
            if (param || meta[attribute]) {
                [meta.subMap(['id']), file(param ?: meta[attribute].replace('${params.igenomes_base}', basepath), checkIfExists: true)]
            }
            else {
                null
            }
        }
        .collect()
}

def get_references_value(references, param, attribute) {
    return references
        .map { meta, _readme ->
            if (param || meta[attribute]) {
                [meta.subMap(['id']), param ?: meta[attribute]]
            }
            else {
                null
            }
        }
        .collect()
}
