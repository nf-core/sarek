// JUST TO TEST THE SUBWORKFLOW

workflow UTILS_REFERENCES {
    take:
    references
    param_file
    param_value
    attribute_file
    attribute_value
    basepath

    main:
    // GIVING up writing a test for the functions, so writing a subworkflow to test it
    references_file = extract_references_file(references, param_file, attribute_file, basepath)
    references_value = extract_references_value(references, param_value, attribute_value)

    emit:
    references_file
    references_value
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS TO EXTRACT REFERENCES FILES OR VALUES FROM THE REFERENCES YAML OR PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def extract_references_file(references, param, attribute, basepath) {
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

def extract_references_value(references, param, attribute) {
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
