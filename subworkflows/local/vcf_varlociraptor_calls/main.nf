//
// VARLOCIRAPTOR CALLS
//
include { VARLOCIRAPTOR_CALLVARIANTS                } from '../modules/nf-core/varlociraptor/callvariants/main'
include { VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES } from '../modules/nf-core/varlociraptor/estimatealignmentproperties/main'
include { VARLOCIRAPTOR_PREPROCESS                  } from '../modules/nf-core/varlociraptor/preprocess/main'

workflow VARLOCIRAPTOR_CALLS {

    take:
    //vcfs
    fasta
    fasta_fai
    cram
    //scenario

    main:
    versions = Channel.empty()

    emit:
    vcfs =  // post processed vcfs

    VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES( cram, fasta, fai)

    input_preprocess = vcfs.join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json)
    input_preprocess.view()

    //VARLOCIRAPTOR_PREPROCESS(vcfs.join(VARLOCIRAPTOR_ESTIMATEALIGNMENTPROPERTIES.out.alignment_properties_json), fasta, fai)

    //VARLOCIRAPTOR_CALLVARIANTS(VARLOCIRAPTOR_PREPROCESS.out.vcf_gz.map{ meta, vcf -> [meta, vcf, []]}, scenario, VARLOCIRAPTOR_PREPROCESS.out.vcf_gz.map{meta, vcf -> meta.id }  )

    versions // channel: [ versions.yml ]
}
