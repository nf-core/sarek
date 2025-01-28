/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS TO EXTRACT REFERENCES FILES OR VALUES FROM THE REFERENCES YAML OR PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def extract_references_file(references, param, attribute) {
    references.map { meta, _readme ->
        if (param || meta[attribute]) {
            return [meta.subMap(['id']), file(param ?: meta[attribute] ?: [])]
        }
        else {
            return null
        }
    }
}

def extract_references_vcf(references, param, vcf_attribute, attribute) {
    references.map { meta, _readme ->
        if (param || meta.vcf[vcf_attribute][attribute]) {
            return [meta.subMap(['id']), file(param ?: meta.vcf[vcf_attribute][attribute] ?: [])]
        }
        else {
            return null
        }
    }
}

def extract_references_value(references, param, attribute) {
    references.map { meta, _readme ->
        if (param || meta[attribute]) {
            return [meta.subMap(['id']), param ?: meta[attribute] ?: []]
        }
        else {
            return null
        }
    }
}
