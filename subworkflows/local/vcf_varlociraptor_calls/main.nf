//
// VARLOCIRAPTOR CALLS
//
include { VARLOCIRAPTOR_CALLVARIANTS                } from '../../../modules/nf-core/varlociraptor/callvariants/main.nf'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES } from '../../../modules/nf-core/varlociraptor/estimatealignmentproperties/main'
include { VARLOCIRAPTOR_PREPROCESS                  } from '../../../modules/nf-core/varlociraptor/preprocess/main'

workflow VARLOCIRAPTOR_CALLS {

    take:
    vcfs
    cram
    fasta
    fasta_fai
    scenario

    main:
    versions = Channel.empty()

    // Possibly pull this out to a different preparation workflow
    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES( cram, fasta.map{ it -> [[id:it.baseName], it]}, fasta_fai.map{ it -> [[id:it.baseName], it]})

    //TODO Scatter candidates with rbt vcfsplit

    input_preprocess = vcfs.map{ meta, vcf -> [meta - meta.subMap('variantcaller'), meta, vcf] }
                            .combine(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json, by: 0)
                            .combine(cram, by: 0)
                            .map{meta_simple, meta, vcf, json, cram, crai -> [meta, cram, crai, vcf, json]}

    VARLOCIRAPTOR_PREPROCESS(input_preprocess, fasta.map{ it -> [[id:it.baseName], it]}, fasta_fai.map{ it -> [[id:it.baseName], it]})

    VARLOCIRAPTOR_CALLVARIANTS(VARLOCIRAPTOR_PREPROCESS.out.vcf_gz.map{meta, vcf -> [meta, vcf, []]}, scenario, [])

    versions = versions.mix(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.versions)
    versions = versions.mix(VARLOCIRAPTOR_PREPROCESS.out.versions)
    versions = versions.mix(VARLOCIRAPTOR_CALLVARIANTS.out.versions)

    emit:
    vcfs = Channel.empty() //VARLOCIRAPTOR_PREPROCESS.out.vcf_gz // post processed vcfs
    versions // channel: [ versions.yml ]
}
