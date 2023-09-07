//
// VARLOCIRAPTOR CALLS
//
include { VARLOCIRAPTOR_CALLVARIANTS                } from '../modules/nf-core/varlociraptor/callvariants/main.nf'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES } from '../modules/nf-core/varlociraptor/estimatealignmentproperties/main'
include { VARLOCIRAPTOR_PREPROCESS                  } from '../modules/nf-core/varlociraptor/preprocess/main'

workflow VARLOCIRAPTOR_CALLS {

    take:
    //vcfs
    cram
    fasta
    fasta_fai
    //scenario

    main:
    versions = Channel.empty()

    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES( cram, fasta, fai)

    input_preprocess = vcfs.join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json)
    input_preprocess.view()

    //VARLOCIRAPTOR_PREPROCESS(vcfs.join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json), fasta, fai)

    //VARLOCIRAPTOR_CALLVARIANTS(VARLOCIRAPTOR_PREPROCESS.out.vcf_gz.map{ meta, vcf -> [meta, vcf, []]}, scenario, VARLOCIRAPTOR_PREPROCESS.out.vcf_gz.map{meta, vcf -> meta.id }  )

    emit:
    vcfs = // post processed vcfs
    versions // channel: [ versions.yml ]


}
