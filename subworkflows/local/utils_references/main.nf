/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS TO EXTRACT REFERENCES FILES OR VALUES FROM THE REFERENCES YAML OR PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def extract_references_file(references, param, attribute) {
    def item = references.map { meta, _readme ->
        if (param || meta[attribute]) {
            [meta.subMap(['id']), file(param ?: meta[attribute].replace('${params.igenomes_base}', params.igenomes_base))]
        }
        else {
            null
        }
    }
    return item.collect()
}

def extract_references_value(references, param, attribute) {
    def item = references.map { meta, _readme ->
        if (param || meta[attribute]) {
            [meta.subMap(['id']), param ?: meta[attribute]]
        }
        else {
            null
        }
    }
    return item.collect()
}
